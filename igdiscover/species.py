"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from sqt.dna import amino_acid_regex

# This is a slightly improved version of the regular expression by
# D’Angelo et al.: The antibody mining toolbox.
# http://dx.doi.org/10.4161/mabs.27105
# The amino-acid version of the expression is:
# [FY][FWVHY]C[ETNGASDRIKVM]X{5,32}W[GAV]
CDR3REGEX = re.compile("""
	(TT[TC] | TA[CT])                            # F or Y
	(TT[CT] | TA[TC] | CA[TC] | GT[AGCT] | TGG)  # any of F, W, V, H, Y
	(TG[TC])                                     # C
	(?P<cdr3>                                    # actual CDR3 starts here
		(([GA][AGCT]) | TC | CG) [AGCT]          # any of ETNGASDRIKVM
		([ACGT]{3}){5,31}                        # between five and 31 codons
	)                                            # end of CDR3
	TGG                                          # W
	G[GCT][GCTA]                                 # G, A or V
	""", re.VERBOSE)


def _build_V_gene_regex():
	r = '('
	r += '|'.join(amino_acid_regex(aa) for aa in 'DV EA EM EV LV QV QL QE VE'.split())
	r += ')' + amino_acid_regex('Q')
	r += '([ACGT]{3})*'

	#r += '(' + amino_acid_regex('F') + '|' + amino_acid_regex('Y') + ')'
	#r += '([ACGT]{3})*'  # any codon
	# beginning of the CDR3 expression
	#r += '(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC|CG)[AGCT]'
	#r += '[ACGT]{4,6}$'
	return re.compile(r)


# TODO this is unused
_V_GENE_REGEX = _build_V_gene_regex()


# TODO this is work in progress
def looks_like_V_gene(s):
	"""
	Check whether the given sequence matches our expectations of how a V gene
	should look like.
	"""
	s = s.upper()
	for start in 'CAGGT CAGCT CAGGA GAGGT GAAGT GACGT GAAAT GTGGA GAGGC GAGAT CTGGT'.split():
		if s.startswith(start):
			return True #break
	else:
		return False
	for end in 'TATTACTGT TTTTACTGT TATTACTGC TATTACTGC TATTGTGCA TATTACTGC TATTATTGT'.split():
		if s[-len(end)-13:].find(end) != -1:
			return True
	return False

	#return bool(_V_GENE_REGEX.match(s))
