"""
For each V gene, plot a histogram of error rates for those sequences that were
assigned to that gene.
"""
import logging
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from .table import read_table

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('errorplot', help=__doc__)
	subparser.set_defaults(func=errorplot_command)
	subparser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene (default: %(default)s)')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	subparser.add_argument('pdf', help='Plot error frequency histograms to this PDF file', default=None)
	return subparser


def plot_error_histogram(group, v_gene, bins=np.arange(20.1)):
	"""
	Plot a histogram of error rates for a specific V gene.

	v_gene -- name of the gene
	"""
	error_rates = list(group.V_SHM)
	z = error_rates.count(0)
	fig = Figure(figsize=(297/25.4, 210/25.4))
	ax = fig.gca()
	ax.set_xlabel('Error rate (%)')
	ax.set_ylabel('Frequency')
	ax.set_title('Gene ' + v_gene, fontsize=18)
	ax.text(0.95, 0.95, '{} sequences with zero differences'.format(z), transform=ax.transAxes, fontsize=15, ha='right', va='top')
	ax.text(0.95, 0.90, '{} different J genes used'.format(len(set(group.J_gene))), transform=ax.transAxes, fontsize=15, ha='right', va='top')

	_ = ax.hist(error_rates, bins=bins)
	return fig


def errorplot_command(args):
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	with PdfPages(args.pdf) as pages:
		for gene, group in table.groupby('V_gene'):
			if len(group) < args.minimum_group_size:
				continue
			fig = plot_error_histogram(group, gene)
			n += 1
			FigureCanvasPdf(fig).print_figure(pages)
	logger.info('%s plots created (rest had too few sequences)', n)