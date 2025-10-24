import dataclasses
import json
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mist.app import model, NAME_DB_INFO
from mist.app.loggers.logger import logger
from mist.app.model import CustomEncoder
from mist.app.query.allelequeryminimap import AlleleQueryMinimap2, MultiStrategy
from mist.app.query.profilequery import ProfileQuery
from mist.app.utils.dependencies import check_dependencies
from mist.version import __version__


class MistCaller:
    """
    Main class to call alleles.
    """

    def __init__(self, dir_db: Path, multi: str, loci: list[str] | None = None, keep_minimap2: bool = False,
                 export_novel: bool = False, min_id_novel: int = 99) -> None:
        """
        Initializes a query.
        :param dir_db: Database directory
        :param multi: Multi-hit strategy
        :param loci: Restrict typing to these loci
        :param keep_minimap2: Keep Minimap2 output
        :param export_novel: If True, novel alleles are exported in FASTA format
        :param min_id_novel: Minimum sequence identity for novel alleles
        :return: None
        """
        check_dependencies(['minimap2'])
        self._dir_db = dir_db
        self._multi = multi
        self._loci = loci
        self._keep_minimap2 = keep_minimap2
        self._export_novel = export_novel
        self._min_id_novel = min_id_novel

    def _results_to_df(self, result_by_locus: dict[str, model.QueryResult]) -> pd.DataFrame:
        """
         Converts the typing results into a dataframe.
        :param result_by_locus: Results by locus.
        :return: DataFrame
        """
        data_out = pd.DataFrame([{
            # Main result
            'locus': locus,
            'allele': res.allele_str,
            'length': ';'.join(str(x.length) for x in res.allele_results) if len(res.allele_results) > 0 else '-',

            # Location in genome
            'contig': ';'.join(x.alignment.seq_id for x in res.allele_results) if len(res.allele_results) > 0 else '-',
            'start': ';'.join(str(x.alignment.start) for x in res.allele_results) if len(res.allele_results) > 0 else '-',
            'end': ';'.join(str(x.alignment.end) for x in res.allele_results) if len(res.allele_results) > 0 else '-',
            'strand': ';'.join(str(x.alignment.strand) for x in res.allele_results) if len(res.allele_results) > 0 else '-',

            # Novel alleles
            'is_novel': model.Tag.NOVEL in res.tags,
            'closest_alleles': ';'.join(res.allele_results[0].closest_alleles) if model.Tag.NOVEL in res.tags else '-',

            # Additional tags
            'tags': ';'.join(t.value for t in res.tags) if len(res.tags) > 0 else '-'}
            for locus, res in result_by_locus.items()
        ])

        # noinspection PyTypeChecker,PyUnresolvedReferences
        nb_detected = (data_out['allele'] != '-').sum()
        nb_novel = data_out['is_novel'].sum()
        logger.info(
            f"Detected {nb_detected}/{len(data_out)} loci ({100 * nb_detected / len(data_out):.2f}%), "
            f"including {nb_novel:,} (potential) novel alleles"
        )
        return data_out

    def _export_novel_allele_seqs(
            self, data_results: pd.DataFrame, res_by_locus: dict[str, model.QueryResult], dir_out: Path) -> None:
        """
        Exports the novel allele sequences.
        :param data_results: DataFrame with typing results
        :param res_by_locus: Results by locus
        :param dir_out: Output directory
        :return: None
        """
        for locus in data_results[data_results['is_novel']]['locus']:
            dir_novel_alleles = dir_out / 'novel_alleles'
            dir_novel_alleles.mkdir(parents=True, exist_ok=True)
            allele_id = res_by_locus[locus].allele_str.replace('*', '')
            with open(dir_novel_alleles / f'{locus}-{allele_id}.fasta', 'w') as handle:
                SeqIO.write(SeqRecord(
                    id=f'{locus}_{allele_id}',
                    description=f"closest_match={','.join(res_by_locus[locus].allele_results[0].closest_alleles)}",
                    seq=Seq(res_by_locus[locus].allele_results[0].sequence),
                ), handle, 'fasta')

    def call_alleles(self, path_fasta: Path, out_json: Path, out_dir: Path, out_tsv: Path, threads: int) -> None:
        """
        Calls the alleles from the input FASTA file.
        :param path_fasta: Input FASTA path
        :param out_json: JSON output path
        :param out_dir: Output directory
        :param out_tsv: Output TSV path
        :param threads: Number of threads
        :return: None
        """
        allele_query = AlleleQueryMinimap2(
            dir_db=self._dir_db,
            dir_out=out_dir,
            multi_strategy=MultiStrategy(self._multi),
            min_id_novel=self._min_id_novel,
            save_minimap2=self._keep_minimap2,
        )
        result_by_locus = allele_query.query(path_fasta, loci=self._loci, threads=threads)

        # Create DataFrame to calculate statistics
        data_results = self._results_to_df(result_by_locus)

        # Query the profiles
        if (self._dir_db / 'profiles.tsv').exists():
            profile_query = ProfileQuery(self._dir_db / 'profiles.tsv')
            profile, pct_match = profile_query.query(result_by_locus)
            logger.info(f'Matching ST: {profile.name} ({pct_match:.2f}% match)')
        else:
            profile, pct_match = None, None

        # Create output files
        self._export_json(result_by_locus, path_out=out_json, profile=profile, pct_match=pct_match)
        if self._export_novel and data_results['is_novel'].sum() > 0:
            logger.info("Exporting novel allele sequences")
            self._export_novel_allele_seqs(data_results, result_by_locus, out_dir)

        # Export TSV output
        if out_tsv is not None:
            data_results.to_csv(out_tsv, sep='\t', index=False)

    def _export_json(
            self, results_by_locus: dict[str, model.QueryResult], path_out: Path, profile: model.Profile,
            pct_match: float) -> None:
        """
        Exports the results in JSON format.
        :param results_by_locus: Result(s) by locus
        :param path_out: Output path
        :param profile: Detected profile
        :param pct_match: Percent match for the profile
        :return: None
        """
        # Parse the database information (if available)
        path_db_info = self._dir_db / NAME_DB_INFO
        if path_db_info.exists():
            with path_db_info.open() as handle:
                data_db = json.load(handle)
        else:
            data_db = None

        with open(path_out, 'w') as handle:
            json.dump({
                'alleles': {locus: dataclasses.asdict(res) for locus, res in results_by_locus.items()},
                'profile': {**dataclasses.asdict(profile), 'pct_match': pct_match} if profile is not None else None,
                'metadata': {
                    'timestamp': datetime.now().isoformat(),
                    'tool_version': __version__,
                    'db': data_db
                }
            }, handle, indent=2, cls=CustomEncoder)
        logger.info(f'Output stored in: {path_out}')
