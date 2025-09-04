import json
from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd

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
            save_minimap2: bool = False
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

            # Check for exact match
            allele = data_locus['hashes'].get(sequenceutils.hash_sequence(seq))
            if allele is None:
                allele = data_locus['hashes'].get(sequenceutils.hash_sequence(sequenceutils.rev_complement(seq)))
            if allele is None:
                continue
            matches.append(model.AlleleResult(
                allele=allele,
                alignment=model.Alignment(
                    seq_id=row['query_name'],
                    start=row['query_start'],
                    end=row['query_end'],
                    strand=row['sstrand']),
                )
            )
        matches.sort(key=lambda res: data_locus['alleles'][res.allele]['idx'])
        if len(matches) > 0:
            return model.QueryResult(merge_results(matches, self._multi_strategy.value), allele_results=matches, tags=[])

        # Check imperfect matches
        logger.debug(f'Screening for imperfect hits for: {locus_name}')
        divergence = df_alignment['tag_dv'].str.rsplit(':', n=1).str[-1].astype(float)
        best_idx = divergence.idxmin()
        seq = self._seq_holder.get_seq(
            df_alignment.loc[best_idx, 'query_name'],
            df_alignment.loc[best_idx, 'query_start'],
            df_alignment.loc[best_idx, 'query_end'],
            df_alignment.loc[best_idx, 'sstrand'],
        )
        if len(seq) == 0:
            return model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.EDGE])
        best_matching = ImperfectMatchDetector(self._dir_db / locus_name)
        try:
            bms = best_matching.retrieve_best_matching(seq, self._min_id_novel)
        except InvalidLengthException:
            return model.QueryResult(model.ALLELE_MISSING, [], tags=[model.Tag.INDEL])

        if len(bms) == 0:
            # No imperfect hits found
            return model.QueryResult(model.ALLELE_MISSING, [], [])
        allele_hash = sequenceutils.hash_sequence(seq)
        return model.QueryResult(
            allele_str=f'n{allele_hash[:6]}',
            allele_results=[model.AlleleResult(
                allele=f'n{allele_hash[:6]}',
                alignment=model.Alignment(
                    seq_id=df_alignment.loc[best_idx, 'query_name'],
                    start=int(df_alignment.loc[best_idx, 'query_start']),
                    end=int(df_alignment.loc[best_idx, 'query_end']),
                    strand=df_alignment.loc[best_idx, 'sstrand'],
                ),
                sequence=seq,
                closest_alleles=bms,
            )],
            tags=[model.Tag.NOVEL]
        )

    def query(
            self,
            path_fasta: Path,
            loci: list[str] | None = None,
            threads: int = 1) -> dict[str, model.QueryResult]:
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
            path_fasta, self._dir_db / 'loci_repr.fasta', include_cigar=False, threads=threads)
        logger.info(f'{len(data_mm2):,} seed alignments')
        if self._save_minimap2:
            data_mm2.to_csv(self._dir_out / 'minimap2_parsed.tsv', sep='\t', index=False)

        # Check for empty results
        if len(data_mm2) == 0:
            logger.warning('No seed alignments found. Please verify that the correct scheme has been specified.')
            return {locus: model.QueryResult(model.ALLELE_MISSING, [], []) for locus in all_loci}

        # Extract loci
        data_mm2['locus'] = data_mm2['sseqid'].str.rsplit('_', n=1).str[0]
        nb_loci = len(data_mm2['locus'].unique())
        logger.info(f"{nb_loci:,}/{len(all_loci)} loci aligned ({100*nb_loci/len(all_loci):.2f}%)")

        # Calculate query string and remove duplicates
        data_mm2[['query_name', 'query_start', 'query_end']] = (
            data_mm2.apply(lambda x: AlleleQueryMinimap2.__get_coord(x), axis=1, result_type='expand'))
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
        return {locus: results_by_locus.get(locus, model.QueryResult(model.ALLELE_MISSING, [], [model.Tag.ABSENT])) for locus in all_loci}
