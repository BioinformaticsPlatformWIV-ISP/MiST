import json
from collections import Counter
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from mist.app import NAME_REPR_INFO, model
from mist.app.loggers.logger import logger
from mist.app.query.bestmatching import ImperfectMatchDetector, InvalidLengthException
from mist.app.query.seqholder import SeqHolder
from mist.app.utils import (
    dbutils,
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

    # Minimap2 default lower bound for the minimizer occurrence cutoff
    MIN_MID_OCC_DEFAULT = 10

    # Factor applied to the largest nb. of representative alleles per locus to obtain the occurrence cutoff
    MID_OCC_FACTOR = 2

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
    def _add_query_coords(data_mm2: pd.DataFrame) -> None:
        """
        Adds the coordinates of the full-length target sequence in the input assembly, by extending the alignment with
        the unaligned parts of the representative allele.
        :param data_mm2: Minimap2 alignment data (updated in place)
        :return: None
        """
        is_plus = data_mm2['sstrand'].isin(['plus', '+'])
        overhang_left = data_mm2['sstart']
        overhang_right = data_mm2['slen'] - data_mm2['send']
        data_mm2['query_name'] = data_mm2['qseqid']
        data_mm2['query_start'] = data_mm2['qstart'] + 1 - overhang_left.where(is_plus, overhang_right)
        data_mm2['query_end'] = data_mm2['qend'] + overhang_right.where(is_plus, overhang_left)

    @staticmethod
    def _extract_exact_match(
        seq: str, alignment: model.Alignment, data_locus: dict[str, Any]
    ) -> model.AlleleResult | None:
        """
        Checks for an exact match.
        :param seq: Allele sequence
        :param alignment: Alignment in the input assembly
        :param data_locus: Locus data
        :return: Match
        """
        allele = data_locus['hashes'].get(sequenceutils.hash_sequence(seq))
        if allele is None:
            allele = data_locus['hashes'].get(sequenceutils.hash_sequence(sequenceutils.rev_complement(seq)))
        if allele is None:
            return None
        return model.AlleleResult(allele=allele, alignment=alignment, length=len(seq))

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
        best_matching = ImperfectMatchDetector(self._dir_db / locus_name, len(seq))
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
        cols = ['query_name', 'query_start', 'query_end', 'sstrand']
        for seq_id, start, end, strand in df_alignment[cols].itertuples(index=False, name=None):
            # Retrieve the full sequence
            seq = self._seq_holder.get_seq(seq_id, start, end)
            if (seq is None) or (len(seq) == 0):
                continue

            # Check for an exact match
            alignment = model.Alignment(seq_id=seq_id, start=start, end=end, strand=strand)
            match_perfect = self._extract_exact_match(seq, alignment, data_locus)
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

    def _get_min_mid_occ(self) -> int:
        """
        Returns the minimum minimizer occurrence cutoff for Minimap2. Minimizers occurring more often than the cutoff
        are ignored. Loci with many representatives (e.g., variable-length repeats) exceed the default and are missed.
        :return: Minimum occurrence cutoff
        """
        path_repr_info = self._dir_db / NAME_REPR_INFO
        if path_repr_info.exists():
            with path_repr_info.open() as handle:
                counts = Counter(json.load(handle)['nb_repr_by_locus'])
        else:
            logger.debug(f'{NAME_REPR_INFO} not found, counting the representative alleles')
            counts = dbutils.count_alleles_by_locus(self._dir_db / 'loci_repr.fasta')
        if len(counts) == 0:
            return AlleleQueryMinimap2.MIN_MID_OCC_DEFAULT
        locus, nb_repr = counts.most_common(1)[0]
        logger.debug(f'Largest nb. of representative alleles: {nb_repr:,} ({locus})')
        return max(AlleleQueryMinimap2.MIN_MID_OCC_DEFAULT, AlleleQueryMinimap2.MID_OCC_FACTOR * nb_repr)

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

        # Seed alignment (use the pre-built index if available)
        path_db = self._dir_db / 'loci_repr.fasta.mni'
        if not path_db.exists():
            logger.warning(f'Minimap2 index not found ({path_db.name}), indexing the representative alleles')
            path_db = self._dir_db / 'loci_repr.fasta'
        min_mid_occ = self._get_min_mid_occ()
        logger.info(f'Performing seed alignment with Minimap2 (min. occurrence cutoff: {min_mid_occ:,})')
        data_mm2 = minimap2utils.align(
            path_fasta, path_db, include_cigar=False, threads=threads, min_mid_occ=min_mid_occ
        )
        logger.info(f'{len(data_mm2):,} seed alignments')
        if self._save_minimap2:
            path_out = self._dir_out / 'minimap2_parsed.tsv'
            data_mm2.to_csv(path_out, sep='\t', index=False)
            logger.info(f'Saved minimap2 output to: {path_out}')

        # Check for empty results
        if len(data_mm2) == 0:
            logger.warning('No seed alignments found. Please verify that the correct scheme has been specified.')
            return {locus: model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.ABSENT]) for locus in all_loci}

        # Extract loci
        data_mm2['locus'] = data_mm2['sseqid'].map(dbutils.get_locus_from_id)
        nb_loci = len(data_mm2['locus'].unique())
        logger.info(f"{nb_loci:,}/{len(all_loci):,} loci aligned ({100 * nb_loci / len(all_loci):.2f}%)")

        # Calculate query string and remove duplicates
        AlleleQueryMinimap2._add_query_coords(data_mm2)
        data_mm2.drop_duplicates(['locus', 'query_name', 'query_start', 'query_end'], keep='first', inplace=True)
        logger.info(f'{len(data_mm2):,} seed alignments (without duplicates)')

        # Save the input sequence in memory
        self._seq_holder = SeqHolder(path_fasta)

        # Perform the querying
        results_by_locus: dict[str, model.QueryResult] = {}
        for locus, data in data_mm2.groupby('locus'):
            if (loci is not None) and (locus not in loci):
                continue
            results_by_locus[str(locus)] = self._process_locus(str(locus), data)
        return {
            locus: results_by_locus.get(locus, model.QueryResult(model.ALLELE_MISSING, [], [model.Tag.ABSENT]))
            for locus in all_loci
        }
