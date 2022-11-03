"""
Group assigned sequences by clonotype

Two sequences have the same clonotype if
- their V and J assignments are the same
- the length of their CDR3 is identical
- the difference between their CDR3s (in terms of mismatches)
  is not higher than a given threshold (by default 1)

The output is a table with one row per clonotype, written to
standard output.

Optionally, a full table of all members (sequences belonging to a clonotype)
can be created with one row per input sequence, sorted by
clonotype, plus an empty line between each group of sequences
that have the same clonotype.

The tables are by default sorted by clonotype, but can instead be sorted
by the group size (number of members of a clonotype).
"""
import itertools
import logging
import time
from contextlib import ExitStack
from collections import Counter

import pandas as pd
from xopen import xopen
from tinyalign import hamming_distance, edit_distance

from ..table import read_table
from ..cluster import hamming_single_linkage
from ..utils import slice_arg

CLONOTYPE_COLUMNS = ['sequence_id', 'count', 'v_call', 'd_call', 'j_call', 'cdr3', 'cdr3_aa',
    'FR1_SHM', 'CDR1_SHM', 'FR2_SHM', 'CDR2_SHM', 'FR3_SHM', 'FR4_SHM',
    'FR1_aa_mut', 'CDR1_aa_mut', 'FR2_aa_mut', 'CDR2_aa_mut', 'FR3_aa_mut', 'V_aa_mut', 'J_aa_mut',
    'V_errors', 'J_errors', 'V_SHM', 'J_SHM', 'barcode', 'VDJ_nt', 'VDJ_aa',
]


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--sort', action='store_true', default=False,
        help='Sort by group size (largest first). Default: Sort by V/D/J gene names')
    arg('--limit', metavar='N', type=int, default=None,
        help='Print out only the first N groups')
    arg('--v-shm-threshold', default=5, type=float,
        help='V SHM threshold for _mindiffrate computations')
    add_clonotyping_cdr3_arguments(arg)
    arg('--no-mindiffrate', dest='mindiffrate', action='store_false', default=True,
        help='Do not add _mindiffrate columns')
    arg('--members', metavar='FILE',
        help='Write member table to FILE')
    arg('--clustered', metavar='FILE',
        help='Fast mode: Only write a table with added clonotype ids to FILE, then exit. '
             'This table contains the same information as the --members table, '
             'but is easier to parse (the clonotype_id column denotes the clonotype a '
             'sequence belongs to). No members or clonotypes tables are created if you use this.')
    arg('table', help='Table with parsed and filtered IgBLAST results')


def add_clonotyping_cdr3_arguments(arg):
    arg('--cdr3-core', default=None,
        type=slice_arg, metavar='START:END',
        help='START:END defines the non-junction region of CDR3 '
            'sequences. Use negative numbers for END to count '
            'from the end. Regions before and after are considered to '
            'be junction sequence, and for two CDR3s to be considered '
            'similar, at least one of the junctions must be identical. '
            'Default: no junction region.')
    arg('--mismatches', default=1, type=float,
        help='No. of allowed mismatches between CDR3 sequences. '
            'Can also be a fraction between 0 and 1 (such as 0.15), '
            'interpreted relative to the length of the CDR3 (minus the front non-core). '
            'Default: %(default)s')
    arg('--aa', default=False, action='store_true',
        help='Count CDR3 mismatches on amino-acid level. Default: Compare nucleotides.')


def main(args):
    run_clonotypes(**vars(args))


def run_clonotypes(
    table,
    sort=False,
    limit=None,
    v_shm_threshold=5,
    aa=False,
    mismatches=1,
    members=None,
    clustered=None,
    cdr3_core=None,
    mindiffrate=True,
):
    logger.info('Reading input table ...')
    usecols = CLONOTYPE_COLUMNS
    # TODO backwards compatibility
    if 'FR1_aa_mut' not in pd.read_table(table, nrows=0).columns:
        usecols = [col for col in usecols if not col.endswith('_aa_mut')]

    table = read_table(table, usecols=usecols)
    logger.info('Read table with %s rows', len(table))
    table.insert(5, 'CDR3_length', table['cdr3'].apply(len))
    table = table[table['CDR3_length'] > 0]
    table = table[table['cdr3_aa'].map(lambda s: '*' not in s)]
    logger.info('After discarding rows with unusable CDR3, %s remain', len(table))
    with ExitStack() as stack:
        if members:
            members_file = stack.enter_context(xopen(members, 'w'))
        else:
            members_file = None

        columns = usecols[:]
        columns.remove('barcode')
        columns.remove('count')
        columns.insert(0, 'count')
        columns.insert(columns.index('cdr3'), 'CDR3_length')
        if not clustered:
            print(*columns, sep='\t')
        members_header = True
        output_header = True
        cdr3_column = 'cdr3_aa' if aa else 'cdr3'

        if mindiffrate:
            barcode_col_index = table.columns.get_loc('barcode')  # insert before this column
            for column in ['cdr3', 'cdr3_aa', 'VDJ_nt', 'VDJ_aa'][::-1]:
                table.insert(barcode_col_index, column + '_mindiffrate', None)

        logger.info("Computing clonotypes ...")
        table = add_clonotype_id(table, mismatches, cdr3_column, cdr3_core)
        if clustered:
            table.to_csv(clustered, sep="\t", index=False)
            logger.info('Found %d clonotypes', table["clonotype_id"].max() + 1)
            return

        grouped = group_by_clonotype(table, mismatches, sort, cdr3_core, cdr3_column)
        started = time.time()
        n = k = 0
        progress_updated = 0
        chunk = []
        for group in itertools.islice(grouped, 0, limit):
            if mindiffrate:
                group = augment_group(group, v_shm_threshold=v_shm_threshold)
            if members_file:
                # We get an intentional empty line between groups since
                # to_csv() already includes a line break
                print(
                    group.drop("clonotype_id", axis=1).to_csv(sep='\t', header=members_header, index=False),
                    file=members_file
                )
                members_header = False

            rep = representative(group)
            print(*[rep[col] for col in columns], sep='\t')
            n += 1
            k += len(group)
            if n % 1000 == 0:
                elapsed = time.time() - started
                if elapsed >= progress_updated + 60:
                    hours = int(elapsed / 3600)
                    minutes = int(elapsed) % 3600 // 60
                    seconds = int(elapsed % 60)
                    logger.info(
                        f"{hours:3d}:{minutes:02d}:{seconds:02d} h:"
                        f" {n} clonotypes and {k} sequences written"
                    )
                    progress_updated = elapsed

    logger.info('%d clonotypes and %d sequences written', n, k)


def group_by_clonotype(table, mismatches, sort, cdr3_core, cdr3_column):
    """
    Yield clonotype groups. Each item is a DataFrame with all the members of the
    clonotype. The input table must have "clonotype_id" column.
    """
    prev_v = None
    groups = []
    for _, group in table.groupby("clonotype_id"):
        v_gene = group["v_call"].iloc[0]
        if prev_v != v_gene:
            logger.info('Processing %s', v_gene)
        prev_v = v_gene
        if sort:
            # When sorting by group size is requested, we need to buffer results
            groups.append(group)
        else:
            yield group
    if sort:
        logger.info("Sorting by size ...")
        yield from sorted(groups, key=len, reverse=True)


def add_clonotype_id(table: pd.DataFrame, mismatches: int, cdr3_column: str, cdr3_core) -> pd.DataFrame:
    """
    Cluster sequences into clonotypes and add a 'clonotype_id' column
    """
    table.insert(
        0, "vjlen_id", table.groupby(["v_call", "j_call", "CDR3_length"]).ngroup()
    )
    table.insert(1, "cdr3_cluster_id", 0)
    table = table.groupby("vjlen_id", group_keys=True).apply(
        assign_cdr3_cluster_id,
        mismatches=mismatches,
        cdr3_core=cdr3_core,
        cdr3_column=cdr3_column,
    )
    table.insert(
        3, "clonotype_id", table.groupby(["vjlen_id", "cdr3_cluster_id"]).ngroup()
    )
    table.drop(["vjlen_id", "cdr3_cluster_id"], axis=1, inplace=True)
    return table


def assign_cdr3_cluster_id(table, mismatches, cdr3_core, cdr3_column):
    """
    Cluster the rows of the table by Hamming distance between
    their CDR3 sequences and use it to add a "cdr3_cluster_id" column.
    """
    # Cluster all unique CDR3s by Hamming distance
    sequences = sorted(set(table[cdr3_column]))
    if len(sequences) == 1:
        return table

    def linked(s, t):
        return is_similar_with_junction(s, t, mismatches, cdr3_core)

    clusters = hamming_single_linkage(sequences, mismatches, linked=linked)

    # Create dict that maps CDR3 sequences to a numeric cluster id
    cluster_ids = dict()
    for cluster_id, cdr3s in enumerate(clusters):
        for cdr3 in cdr3s:
            cluster_ids[cdr3] = cluster_id

    # Assign cluster id to each row
    table["cdr3_cluster_id"] = table[cdr3_column].apply(lambda cdr3: cluster_ids[cdr3])
    return table


def is_similar_with_junction(s, t, mismatches, cdr3_core):
    """
    Return whether strings s and t have at most the given number of mismatches
    *and* have at least one identical junction.
    """
    # TODO see issue #81
    if len(s) != len(t):
        return False
    if 0 < mismatches < 1:
        delta = cdr3_core.start if cdr3_core is not None else 0
        distance_ok = hamming_distance(s, t) <= (len(s) - delta) * mismatches
    else:
        distance_ok = hamming_distance(s, t) <= mismatches
    if cdr3_core is None:
        return distance_ok
    return distance_ok and (
            (s[:cdr3_core.start] == t[:cdr3_core.start]) or
            (s[cdr3_core.stop:] == t[cdr3_core.stop:]))


def representative(table):
    """
    Given a table with members of the same clonotype, return a representative
    as a dict.
    """
    n = len(table)
    if n == 1:
        return table.iloc[0]
    elif n == 2:
        result = table.iloc[0].copy()
    else:
        c = Counter()
        for row in table.itertuples():
            c[row.VDJ_nt] += row.count
        most_common_vdj_nt = c.most_common(1)[0][0]
        result = table[table['VDJ_nt'] == most_common_vdj_nt].iloc[0]
    result.at['count'] = table['count'].sum()
    return result


def augment_group(table, v_shm_threshold=5, suffix='_mindiffrate'):
    """
    Fill in the _mindiffrate columns. These contain percentage difference of VDJ_nt, VDJ_aa, cdr3,
    cdr3_aa to the least mutated (in terms of V_SHM) sequence in this group.
    """
    if table.empty:
        return table

    # Find row whose V is least mutated
    root = table.loc[table['V_SHM'].idxmin()]
    if root['V_SHM'] > v_shm_threshold:
        return table

    for column in ['cdr3', 'cdr3_aa', 'VDJ_nt', 'VDJ_aa']:
        root_seq = root[column]
        table[column + suffix] = table[column].apply(lambda s:
            round(edit_distance(root_seq, s, maxdiff=int(0.2 * len(root_seq))) / len(root_seq) * 100., 1)
        )

    return table
