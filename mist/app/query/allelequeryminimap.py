import json
from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd
import polars as pl
import orjson

from mist.app import model
from mist.app.loggers.logger import logger
from mist.app.query.bestmatching import ImperfectMatchDetector, InvalidLengthException
from mist.app.query.seqholder import SeqHolder
from mist.app.utils import (
    minimap2utils,
    sequenceutils,
    unique_preserve_order,
)


class MultiStrategy(Enum):
    """
    Strategy to handle multiple perfect hits.
    """

    ALL = 'all'
    FIRST = 'first'
    LONGEST = 'longest'


def merge_results(results: list[model.AlleleResult], multi: str) -> str:
    """
    Merges allele results into a single allele string.
    :param results: Results
    :param multi: Multi-hit strategy
    :return: Allele string
    """
    if len(results) == 0:
        return '-'
    if len(results) == 1:
        return results[0].allele
    multi_hit_strategy = MultiStrategy(multi)
    if multi_hit_strategy == MultiStrategy.ALL:
        unique_alleles = unique_preserve_order([r.allele for r in results])
        return '__'.join(unique_alleles)
    elif multi_hit_strategy == MultiStrategy.LONGEST:
        max_length = max(res.alignment.length for res in results)
        unique_alleles = unique_preserve_order([r.allele for r in results if r.alignment.length == max_length])
        return '__'.join(unique_alleles)
    else:
        raise ValueError(f'Invalid multi-hit strategy: {multi}')


class AlleleQueryMinimap2:
    """
    Queries alleles using Minimap2.
    """

    def __init__(
        self,
        dir_db: Path,
        dir_out: Optional[Path] = None,
        multi_strategy: MultiStrategy = MultiStrategy.LONGEST,
        min_id_novel: int = 99,
        save_minimap2: bool = False,
    ) -> None:
        """
        Initializes this class.
        :param dir_db: Path to the database
        :param dir_out: Output directory to store BLAST output (optional)
        :param multi_strategy: Strategy to handle multiple perfect hits
        :param min_id_novel: Minimum % identity for novel alleles
        :param save_minimap2: If True, the Minimap2 output is stored
        :return: None
        """
        self._dir_db = dir_db
        self._dir_out = dir_out
        self._multi_strategy = multi_strategy
        self._min_id_novel = min_id_novel
        self._seq_holder: SeqHolder | None = None
        self._save_minimap2 = save_minimap2

    @staticmethod
    def __get_coord(x: pd.Series) -> tuple[str, int, int]:
        """
        Constructs the coordinate string to retrieve the target sequence from the input assembly.
        :param x: BLAST output records
        :return: Coordinate string
        """
        if x['sstrand'] in ('plus', '+'):
            overhang_left = x['sstart']
            overhang_right = x['slen'] - x['send']
            return x['qseqid'], x['qstart'] + 1 - overhang_left, x['qend'] + overhang_right
        else:
            overhang_left = x['sstart']
            overhang_right = x['slen'] - x['send']
            return x['qseqid'], x['qstart'] + 1 - overhang_right, x['qend'] + overhang_left

    def _extract_exact_match(self, seq: str, row: pd.Series, data_locus: pd.DataFrame) -> model.AlleleResult | None:
        """
        Checks for an exact match.
        :param seq: Allele sequence
        :param row: Alignment record
        :param data_locus: Locus data
        :return: Match
        """
        allele = data_locus['hashes'].get(sequenceutils.hash_sequence(seq))
        if allele is None:
            allele = data_locus['hashes'].get(sequenceutils.hash_sequence(sequenceutils.rev_complement(seq)))
        if allele is None:
            return None
        return model.AlleleResult(
            allele=allele,
            alignment=model.Alignment(
                seq_id=row['query_name'], start=row['query_start'], end=row['query_end'], strand=row['sstrand']
            ),
            length=len(seq),
        )

    def _extract_partial_match(self, df_alignment: pd.DataFrame, locus_name: str) -> model.QueryResult | None:
        """
        Checks for an exact match.
        :param df_alignment: Minimap2 alignment data
        :param locus_name: Locus name
        :return: Match
        """
        divergence = df_alignment['tag_dv'].str.rsplit(':', n=1).str[-1].astype(float)
        best_idx = divergence.idxmin()
        seq = self._seq_holder.get_seq(
            df_alignment.loc[best_idx, 'query_name'],
            df_alignment.loc[best_idx, 'query_start'],
            df_alignment.loc[best_idx, 'query_end'],
            df_alignment.loc[best_idx, 'sstrand'],
        )
        if len(seq) == 0:
            # Allele coordinates extend beyond assembly boundaries —> likely contig edge
            return model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.EDGE])

        # Screen for imperfect matches
        best_matching = ImperfectMatchDetector(self._dir_db / locus_name)
        try:
            seq_ids_closest = best_matching.retrieve_best_matching(seq, self._min_id_novel)
        except InvalidLengthException:
            return model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.INDEL])

        # No imperfect hits to existing alleles found
        if len(seq_ids_closest) == 0:
            return model.QueryResult(model.ALLELE_MISSING, [], [])

        # Potential novel allele
        allele_hash = sequenceutils.hash_sequence(seq, rev_comp=False)
        allele_hash_shown = f'*{allele_hash[:4]}'
        return model.QueryResult(
            allele_str=allele_hash_shown,
            allele_results=[
                model.AlleleResult(
                    allele=allele_hash_shown,
                    alignment=model.Alignment(
                        seq_id=df_alignment.loc[best_idx, 'query_name'],
                        start=int(df_alignment.loc[best_idx, 'query_start']),
                        end=int(df_alignment.loc[best_idx, 'query_end']),
                        strand=df_alignment.loc[best_idx, 'sstrand'],
                    ),
                    length=len(seq),
                    sequence=seq,
                    closest_alleles=seq_ids_closest,
                )
            ],
            tags=[model.Tag.NOVEL],
        )

    def _process_locus_polars(self, locus_name: str, df_alignment: pl.DataFrame) -> model.QueryResult:
        """
        Types the input locus.
        :param locus_name: Locus name
        :param df_alignment: Alignment data
        :return: Results for the locus
        """
        # Parse database information
        dir_locus = self._dir_db / locus_name
        with open(dir_locus / 'mist_db.json', 'rb') as handle: # Open in binary
            data_locus = orjson.loads(handle.read())

        # Process seed alignments
        matches = []
        for row in df_alignment.iter_rows(named=True):
            # row is now a dictionary: {'query_name': '...', 'query_start': 123, ...}
            seq = self._seq_holder.get_seq(row['query_name'], row['query_start'], row['query_end'])

            if not seq:  # More concise than (seq is None) or (len(seq) == 0)
                continue

            match_perfect = self._extract_exact_match(seq, row, data_locus)
            if match_perfect is not None:
                matches.append(match_perfect)

        # Check if perfect matches have been found
        matches.sort(key=lambda res: data_locus['alleles'][res.allele]['idx'])
        if len(matches) > 0:
            return model.QueryResult(
                merge_results(matches, self._multi_strategy.value),
                allele_results=matches,
                tags=[model.Tag.MULTI] if len(matches) > 1 else [model.Tag.EXACT],
            )

        # Check imperfect matches
        logger.debug(f'Screening for imperfect hits for: {locus_name}')
        return self._extract_partial_match(df_alignment, locus_name)

    def _process_locus(self, locus_name: str, df_alignment: pd.DataFrame) -> model.QueryResult:
        """
        Types the input locus.
        :param locus_name: Locus name
        :param df_alignment: Alignment data
        :return: Results for the locus
        """
        # Parse database information
        dir_locus = self._dir_db / locus_name
        with open(dir_locus / 'mist_db.json') as handle:
            data_locus = json.load(handle)

        # Process seed alignments
        matches = []
        for _, row in df_alignment.iterrows():
            # Retrieve the full sequence
            seq = self._seq_holder.get_seq(row['query_name'], row['query_start'], row['query_end'])
            if (seq is None) or (len(seq) == 0):
                continue

            # Check for an exact match
            match_perfect = self._extract_exact_match(seq, row, data_locus)
            if match_perfect is not None:
                matches.append(match_perfect)

        # Check if perfect matches have been found
        matches.sort(key=lambda res: data_locus['alleles'][res.allele]['idx'])
        if len(matches) > 0:
            return model.QueryResult(
                merge_results(matches, self._multi_strategy.value),
                allele_results=matches,
                tags=[model.Tag.MULTI] if len(matches) > 1 else [model.Tag.EXACT],
            )

        # Check imperfect matches
        logger.debug(f'Screening for imperfect hits for: {locus_name}')
        return self._extract_partial_match(df_alignment, locus_name)

    def query(self, path_fasta: Path, loci: list[str] | None = None, threads: int = 1) -> dict[str, model.QueryResult]:
        """
        Queries the database with the given FASTA file.
        :param path_fasta: Input FASTA file
        :param loci: List of target loci
        :param threads: Threads (number of threads to use)
        :return: Matching allele(s)
        """
        # Retrieve all loci
        with open(self._dir_db / 'loci.txt') as handle:
            all_loci = [l.strip() for l in handle]

        # Seed alignment
        logger.info('Performing seed alignment with Minimap2')
        data_mm2 = minimap2utils.align(
            path_fasta, self._dir_db / 'loci_repr.fasta', include_cigar=False, threads=threads
        )
        data_mm2_pl = pl.from_pandas(data_mm2)

        logger.info(f'{len(data_mm2_pl):,} seed alignments')
        if self._save_minimap2:
            path_out = self._dir_out / 'minimap2_parsed.tsv'
            data_mm2.to_csv(path_out, sep='\t', index=False)
            logger.info(f'Saved minimap2 output to: {path_out}')

        # Check for empty results
        if len(data_mm2) == 0:
            logger.warning('No seed alignments found. Please verify that the correct scheme has been specified.')
            return {locus: model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.ABSENT]) for locus in all_loci}

        data_mm2_pl = data_mm2_pl.with_columns(
            locus=pl.col("sseqid")
            .str.split("_")  # Split into a list: ["ABC", "DEF", "GHI"]
            .list.slice(0, pl.col("sseqid").str.count_matches("_"))  # Keep all but last
            .list.join("_")  # Join back: "ABC_DEF"
        )
        nb_loci = data_mm2_pl['locus'].n_unique()
        logger.info(f"{nb_loci:,}/{len(all_loci):,} loci aligned ({100 * nb_loci / len(all_loci):.2f}%)")

        overhang_left = pl.col("sstart")
        overhang_right = pl.col("slen") - pl.col("send")

        # 2. Apply the conditional logic using expressions
        data_mm2_pl = data_mm2_pl.with_columns(
            [
                pl.col("qseqid").alias("query_name"),
                pl.when(pl.col("sstrand").is_in(["plus", "+"]))
                .then(pl.col("qstart") + 1 - overhang_left)
                .otherwise(pl.col("qstart") + 1 - overhang_right)
                .alias("query_start"),
                pl.when(pl.col("sstrand").is_in(["plus", "+"]))
                .then(pl.col("qend") + overhang_right)
                .otherwise(pl.col("qend") + overhang_left)
                .alias("query_end"),
            ]
        )
        data_mm2_pl = data_mm2_pl.unique(subset=['locus', 'query_name', 'query_start', 'query_end'], keep='first')
        logger.info(f'{len(data_mm2_pl):,} seed alignments (without duplicates)')

        # Save the input sequence in memory
        self._seq_holder = SeqHolder(path_fasta)

        # Perform the querying
        results_by_locus: dict[str, model.QueryResult] = {}
        if loci is not None:
            data_mm2_pl = data_mm2_pl.filter(pl.col("locus").is_in(loci))
        groups = data_mm2_pl.partition_by("locus", as_dict=True)

        for (locus_tuple,), df_locus in groups.items():
            locus_name = str(locus_tuple)
            results_by_locus[locus_name] = self._process_locus_polars(locus_name, df_locus)

        return {
            locus: results_by_locus.get(locus, model.QueryResult(model.ALLELE_MISSING, [], [model.Tag.ABSENT]))
            for locus in all_loci
        }
