# kate: syntax Python;
"""
Required modules:

module use /proj/b2014342/sw/modules
module load igypipe

These are the main files that the pipeline creates, in the order in which they
are created:

merged.fastq.gz -- merged reads
trimmed.fastq.gz -- primers removed from merged reads
filtered.fasta  -- too short sequences removed, converted to FASTA
unique.fasta -- collapsed sequences
unique.igblast.txt -- raw IgBLAST output
unique.table.tab -- result of parsing IgBLAST output
groups.tab -- sequences grouped by barcode
consensus.fasta -- contains one consensus sequence for each group
consensus.igblast.txt -- consensus sequences sent through IgBLAST
consensus.table.tab -- result of parsing IgBLAST output
"""
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import FigureCanvasPdf
import seaborn as sns
from sqt.dna import reverse_complement
import igypipe

# Set defaults for some configurable values.

MERGE_PROGRAM = 'pear'

# Limit processing to this number of reads. Set to None to process all.
LIMIT = None

# Which program to use for clustering. This can be either a path to the
# usearch binary or vsearch.
CLUSTER_PROGRAM = 'vsearch'

MULTIALIGN_PROGRAM = 'muscle-medium'

# Filter out reads that have more than this number of expected errors.
# Set to None to disable.
MAXIMUM_EXPECTED_ERRORS = None

# Whether to trim primers. Can be set to True or False.
TRIM_PRIMERS = False

# IgBLAST mismatch penalty (can be -1, -2, -3, -4 or None)
MISMATCH_PENALTY = None

BARCODE_LENGTH = 0

try:
	include: "pipeline.conf"
except WorkflowError:
	sys.exit("Pipeline configuration file 'pipeline.conf' not found. Please create it!")

if not REVERSE_PRIMERS:
	sys.exit("The list of REVERSE_PRIMERS is empty. This will currently not work correctly.")

if not FORWARD_PRIMERS:
	sys.exit("The list of FORWARD_PRIMERS is empty. This will currently not work correctly.")

# This command is run before every shell command and helps to catch errors early
shell.prefix("set -euo pipefail;")

if LIMIT:
	READS1 = 'reads-{}.1.fastq'.format(LIMIT)
	READS2 = 'reads-{}.2.fastq'.format(LIMIT)
else:
	READS1 = 'reads.1.fastq'
	READS2 = 'reads.2.fastq'


rule all:
	input:
		READS1, READS2,
		#expand("fastqc/reads.{r}.zip", r=(1, 2)),
		"stats/readlengthhisto.pdf",
		"stats/barcodes.txt",
		"stats/counts.txt",
		"stats/unique.correlationVJ.pdf",
		"stats/consensus.correlationVJ.pdf",
		'stats/unique.shmhistograms.pdf',
		'stats/consensus.shmhistograms.pdf',
		"unique.table.tab",
		"consensus.table.tab",
		"unique.v_usage.tab",
		"unique.v_usage.pdf",
		"consensus.v_usage.tab",
		"consensus.v_usage.pdf"


if LIMIT:
	rule limit_reads:
		output: 'reads-{}.{{nr}}.fastq'.format(LIMIT)
		input: 'reads.{nr}.fastq.gz'
		shell:
			'sqt-fastqmod --limit {LIMIT} {input} > {output}'
else:
	rule ungzip:
		output: "{file}.fastq"
		input: "{file}.fastq.gz"
		shell: "zcat {input} > {output}"


if MERGE_PROGRAM == 'flash':
	rule flash_merge:
		"""Use FLASH to merge paired-end reads"""
		output: "merged.fastq.gz"
		input: READS1, READS2
		resources: time=60
		threads: 8
		log: 'flash.log'
		shell:
			# -M: maximal overlap (2x300, 420-450bp expected fragment size)
			"flash -t {threads} -c -M {FLASH_MAXIMUM_OVERLAP} {input} 2> >(tee {log} >&2) | pigz > {output}"
elif MERGE_PROGRAM == 'pear':
	rule pear_merge:
		"""Use pear to merge paired-end reads"""
		output: 'pear.unassembled.forward.fastq', 'pear.unassembled.reverse.fastq', 'pear.discarded.fastq', fastq="merged.fastq.gz"
		input: READS1, READS2
		resources: time=60
		threads: 8
		shell:
			r"""
			pear -j {threads} -f {input[0]} -r {input[1]} -o pear && \
			gzip < pear.assembled.fastq > {output.fastq} && rm pear.assembled.fastq
			"""
else:
	sys.exit("MERGE_PROGRAM not recognized")


rule read_length_histogram:
	# TODO on which data should this be created?
	output:
		txt="stats/readlengthhisto.txt",
		pdf="stats/readlengthhisto.pdf"
	input:
		fastq="merged.fastq.gz"
	shell:
		"sqt-readlenhisto --plot {output.pdf} {input}  > {output.txt}"


rule barcodes:
	"""Print out number of random barcodes in the library
	TODO
	- run this only if random barcodes are actually used
	- make sure that a stranded protocol is used
	"""
	output: txt="stats/barcodes.txt"
	input: fastq="merged.fastq.gz"
	shell:
		"""
		zcat {input} | awk 'NR%4==2 {{print substr($1,1,12)}}' | grep -v N | sort -u | wc -l > {output}
		"""

rule stats_numbers:
	output: txt="stats/counts.txt"
	input:
		reads=READS1,
		merged="merged.fastq.gz",
		unique="unique.fasta",
		unique_table="unique.table.tab",
		consensus_table="consensus.table.tab",
	shell:
		"""
		echo -n "Number of paired-end reads: " > {output}
		awk 'END {{ print NR/4 }}' {input.reads} >> {output}
		echo -n "Number of barcodes (looking at 1st read in pair): " >> {output}
		awk 'NR % 4 == 2 {{ print substr($1, 1, 12) }}' {input.reads} | sort -u | grep -v N | wc -l >> {output}
		echo -n "Number of merged sequences: " >> {output}
		zcat {input.merged} | awk 'END {{ print NR/4 }}' >> {output}
		echo -n "Number of barcodes in merged sequences: " >> {output}
		zcat {input.merged} | awk 'NR % 4 == 2 {{ print substr($1, 1, 12) }}' | sort -u | grep -v N | wc -l >> {output}
		echo -n "Number of unique sequences: " >> {output}
		grep -c '^>' unique.fasta >> {output}
		echo -n "Number of barcodes in unique sequences: " >> {output}
		grep -A 1 '^>' {input.unique} | awk '!/^>/ && $1 != "--" {{ print substr($1,1,12) }}' | sort -u | grep -v N | wc -l >> {output}
		echo -n "Number of barcodes in unique table: " >> {output}
		cut -f35 {input.unique_table} | sed 1d | awk '{{print substr($1, 1, 12) }}' | sort -u | grep -v N | wc -l >>  {output}
		echo -n "Number of sequences in final consensus table: " >> {output}
		sed 1d {input.consensus_table} | wc -l >> {output}
		"""


rule stats_correlation_V_J:
	output:
		pdf="stats/{base}.correlationVJ.pdf"
	input:
		table="{base}.table.tab"
	run:
		table = igypipe.read_table(input.table)
		fig = Figure(figsize=(30/2.54, 21/2.54))
		ax = fig.gca()
		ax.set_xlabel('V%SHM')
		ax.set_ylabel('J%SHM')
		ax.scatter(table['V_SHM'], table['J_SHM'])
		ax.set_xlim(left=-0.1)
		ax.set_ylim(bottom=-0.1)
		ax.set_title('Correlation between V%SHM and J%SHM')
		FigureCanvasPdf(fig).print_figure(output.pdf)


rule stats_plot_shmhistograms:
	output:
		pdf='stats/{base}.shmhistograms.pdf',
	input:
		table='{base}.table.tab'
	shell:
		'igypipe singledisco --plot {output.pdf} {input.table}'


# Adjust the primer sequences so they are correctly reverse-complemented.
# When a forward primer fwd and a reverse primer rev are given, then we need to
# search for:
# * fwd in the beginning, revcomp(rev) in the end
# * rev in the beginning, revcomp(fwd) in the end
PRIMERS = [
	FORWARD_PRIMERS[:], [reverse_complement(seq) for seq in REVERSE_PRIMERS]
]
PRIMERS[0].extend(REVERSE_PRIMERS)
PRIMERS[1].extend(reverse_complement(seq) for seq in FORWARD_PRIMERS)

rule trim_primers_fivep:
	output: fastq="trimmed-fivep.fastq.gz"
	input: fastq="merged.fastq.gz"
	resources: time=60
	log: 'cutadapt-fiveprime.log'
	params:
		five_p=" ".join("-g ^{}".format(seq) for seq in PRIMERS[0])
	shell:
		r"""
		cutadapt --discard-untrimmed {params.five_p} -o {output.fastq} {input.fastq} | tee {log}
		"""


rule trim_primers_threep:
	output: fastq="trimmed.fastq.gz"
	input: fastq="trimmed-fivep.fastq.gz"
	resources: time=60
	log: 'cutadapt-threeprime.log'
	params:
		three_p=" ".join("-a {}$".format(seq) for seq in PRIMERS[1])
	shell:
		r"""
		cutadapt --discard-untrimmed {params.three_p} -o {output.fastq} {input.fastq} | tee {log}
		"""


rule fastqc:
	output:
		zip='fastqc/{file}.zip',
		png='fastqc/{file}/Images/per_base_quality.png',
		html='fastqc/{file}_fastqc.html'
	input: fastq='{file}.fastq'
	shell:
		r"""
		rm -rf fastqc/{wildcards.file}/ fastqc/{wildcards.file}_fastqc/ && \
		fastqc -o fastqc {input} && \
		mv fastqc/{wildcards.file}_fastqc.zip {output.zip} && \
		unzip -o -d fastqc/ {output.zip} && \
		mv fastqc/{wildcards.file}_fastqc/ fastqc/{wildcards.file}
		"""


rule fastq_to_fasta:
	"""
	* Convert from FASTQ to FASTA
	* Remove low-quality sequences
	* Discard too short sequences
	"""
	output: fasta="filtered.fasta"
	input: fastq="trimmed.fastq.gz" if TRIM_PRIMERS else "merged.fastq.gz"
	params: max_errors="--max-errors {MAXIMUM_EXPECTED_ERRORS}" if MAXIMUM_EXPECTED_ERRORS is not None else ""
	shell:
		"sqt-fastqmod {params.max_errors} --minimum-length {MINIMUM_MERGED_READ_LENGTH} --fasta {input.fastq} > {output.fasta}"


rule dereplicate:
	"""Collapse identical sequences with VSEARCH"""
	output: fasta="unique.fasta"
	input: fasta="filtered.fasta"
	shell:
		"""vsearch --derep_fulllength {input.fasta} --strand both --output {output.fasta} --sizeout"""


rule cluster:
	"""
	TODO Daniel ran this three times (at 99%, 98% and 97% identity) in order to
	avoid a specific type of misclustering.
	"""
	output:
		fasta="clustered.fasta",  # centroids
		uc="clustered.uc"
	input: fasta="unique.fasta"
	resources: time=36*60, mem=32000
	threads: 4
	shell:
		# TODO -idprefix 5?
		r"""
		{CLUSTER_PROGRAM} -threads {threads} -cluster_fast {input.fasta} -id 0.97 -uc {output.uc} \
			-idprefix 5 -sizeout --centroids {output.fasta}
		"""


rule makeblastdb:
	input: fasta="database/{organism}_{gene}.fasta"
	output: "database/{organism}_{gene}.nhr" # and nin nog nsd nsi nsq
	params: dbname="database/{organism}_{gene}"
	log: 'database/{organism}_{gene}.log'
	threads: 100  # force to run as only job
	shell:
		r"""
		makeblastdb -parse_seqids -dbtype nucl -in {input.fasta} -out {params.dbname} >& {log}
		grep '^Error:' {log} && {{ echo "makeblastdb failed when creating {params.dbname}"; false; }} || true
		"""


rule igypipe_igblast:
	output:
		txt="{base}.igblast.txt"
	input:
		fasta="{base}.fasta",
		db_v="database/{species}_V.nhr".format(species=SPECIES),
		db_d="database/{species}_D.nhr".format(species=SPECIES),
		db_j="database/{species}_J.nhr".format(species=SPECIES)
	params:
		penalty='--penalty {}'.format(MISMATCH_PENALTY) if MISMATCH_PENALTY is not None else ''
	threads: 16
	shell:
		#-auxiliary_data $IGDATA/optional_file/{SPECIES}_gl.aux
		r"""
		igypipe igblast --threads {threads} {params.penalty} --species {SPECIES} database/ {input.fasta} > {output.txt}
		"""


rule igypipe_parse:
	output:
		tab="{base}.table.tab"
	input:
		txt="{base}.igblast.txt",
		fasta="{base}.fasta"
	params:
		dirname=os.path.basename(os.getcwd())
	shell:
		"igypipe parse --rename {params.dirname}_ {input.txt} {input.fasta} > {output.tab}"


rule igypipe_group:
	"""Group by barcode"""
	output:
		pdf="stats/groupsizes.pdf",
		tab="groups.tab",
		fasta="consensus.fasta"
	input:
		tab="unique.table.tab"
	shell:
		"igypipe group --program={MULTIALIGN_PROGRAM} --barcode-length {BARCODE_LENGTH} --plot-sizes {output.pdf} --groups-output {output.tab} {input.tab} > {output.fasta}"


rule count_and_plot:
	output:
		plot="{base}.v_usage.pdf",
		counts="{base}.v_usage.tab"
	input:
		v_reference="database/{SPECIES}_V.fasta".format(SPECIES=SPECIES),
		tab="{base}.table.tab"
	shell:
		"igypipe count --reference {input.v_reference} {input.tab} {output.plot} > {output.counts}"

