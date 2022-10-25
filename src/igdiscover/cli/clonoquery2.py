"""
Query one or more tables of assigned sequences by clonotype

"""
import itertools
import logging
import time
from contextlib import ExitStack

import pandas as pd
from xopen import xopen

from ..table import read_table
from .clonotypes import (
    CLONOTYPE_COLUMNS,
    add_clonotyping_cdr3_arguments,
    group_by_clonotype,
    representative, add_clonotype_id,
)
from .clonoquery import collect, read_and_preprocess_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    add_clonotyping_cdr3_arguments(arg)
    arg(
        "querytable",
        help="Query table with IgBLAST results (assigned.tab or filtered.tab)",
    )
    arg(
        "reftables",
        nargs="+",
        help="Reference tables with parsed and filtered "
        "IgBLAST results (filtered.tab)",
    )


def main(args):
    run_clonotypes2(**vars(args))


def run_clonotypes2(
    querytable,
    reftables,
    clustered: str=None,
    cdr3_core=None,
    mismatches=1,
    aa=False,
):
    logger.info("Reading tables ...")
    reftable_paths = reftables
    usecols = CLONOTYPE_COLUMNS
    reftables = [read_table(path, usecols=usecols) for path in reftable_paths]
    querytable = read_and_preprocess_table(querytable, filter_empty=True)

    if len(querytable) > len(reftables[-1]):
        logger.warning(
            "The last reference table is smaller than the "
            "query table! Did you swap query and reference? (Query needs to be given first.)"
        )

    reftable = pd.concat(reftables, keys=range(len(reftables)), names=["library"])
    logger.info("Read %d tables with %s rows in total", len(reftables), len(reftable))
    del reftables
    cdr3_column = "cdr3_aa" if aa else "cdr3"

    reftable.insert(5, "CDR3_length", reftable["cdr3"].apply(len))
    reftable = reftable[reftable["CDR3_length"] > 0]
    reftable = reftable[reftable["cdr3_aa"].map(lambda s: "*" not in s)]
    logger.info("After discarding rows with unusable CDR3, %s remain", len(reftable))

    logger.info("Clustering into clonotypes ...")
    reftable = add_clonotype_id(reftable, mismatches, cdr3_column, cdr3_core)

    if clustered:
        reftable.to_csv(clustered, sep="\t", index=False)
        logger.info('Found %d clonotypes', reftable["clonotype_id"].max() + 1)

    with ExitStack() as stack:

        started = time.time()
        n = k = 0
        progress_updated = 0
        for group in itertools.islice(grouped, 0, limit):
            rep = representative(group)
            print(*[rep[col] for col in columns], sep="\t")
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
    logger.info("%d clonotypes and %d sequences written", n, k)


def query(querytable, reftable, mismatches, cdr3_core_slice, cdr3_column):
    """
    Find the clonotypes of the querytable in the reftable.

    Find all queries from the querytable in the reftable.

    Yield tuples (query_rows, similar_rows) where the query_rows is a list
    with all the rows that have the same result. similar_rows is a DataFrame
    whose rows are the ones matching the query.
    """

    # The vjlentype is a "clonotype without CDR3 sequence" (only V, J, CDR3 length)
    # Determine set of vjlentypes to query
    query_vjlentypes = defaultdict(list)
    for row in querytable.itertuples():
        vjlentype = (row.v_call, row.j_call, len(row.cdr3))
        query_vjlentypes[vjlentype].append(row)

    groupby = ['v_call', 'j_call', 'CDR3_length']
    for vjlentype, vjlen_group in reftable.groupby(groupby):
        # (v_gene, j_gene, cdr3_length) = vjlentype
        if vjlentype not in query_vjlentypes:
            continue

        # Collect results for this vjlentype. The results dict
        # maps row indices (into the vjlen_group) to each query_row,
        # allowing us to group identical results together.
        results = defaultdict(list)
        for query_row in query_vjlentypes.pop(vjlentype):
            cdr3 = getattr(query_row, cdr3_column)
            # Save indices of the rows that are similar to this query
            indices = tuple(index for index, r in enumerate(vjlen_group.itertuples())
                if is_similar_with_junction(cdr3, getattr(r, cdr3_column), mismatches, cdr3_core_slice))
            results[indices].append(query_row)

        # Yield results, grouping queries that lead to the same result
        for indices, query_rows in results.items():
            if not indices:
                for query_row in query_rows:
                    yield ([query_row], reftable.head(0))
                continue

            similar_group = vjlen_group.iloc[list(indices), :].copy()
            yield (query_rows, similar_group)

    # Yield result tuples for all the queries that have not been found
    for queries in query_vjlentypes.values():
        for query_row in queries:
            yield ([query_row], reftable.head(0))
