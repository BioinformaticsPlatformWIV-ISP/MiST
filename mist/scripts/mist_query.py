#!/usr/bin/env python
import argparse
import dataclasses
import json
from collections.abc import Sequence
from datetime import datetime
from importlib.resources import files
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mist.app import model
from mist.app.loggers.logger import initialize_logging, logger
from mist.app.model import CustomEncoder
from mist.app.query.allelequeryminimap import AlleleQueryMinimap2, MultiStrategy
from mist.app.query.profilequery import ProfileQuery
from mist.app.utils.dependencies import check_dependencies
from mist.version import __version__


class MistQuery:
    """
    Class to query a scheme.
    """

    def __init__(self, dir_db: Path, multi: str, loci: list[str] | None = None, keep_minimap2: bool = False,
                 export_novel: bool = False) -> None:
        """
        Initializes a query.
        :param dir_db: Database directory
        :param multi: Multi-hit strategy
        :param loci: Restricts typing to these loci
        :param keep_minimap2: Keep Minimap2 output
        :param export_novel: If True, novel alleles are exported in FASTA format
        :return: None
        """
        self._dir_db = dir_db
        self._multi = multi
        self._loci = loci
        self._keep_minimap2 = keep_minimap2
        self._export_novel = export_novel

    def _results_to_df(self, result_by_locus: dict[str, model.QueryResult]) -> pd.DataFrame:
        """
         Converts the typing results into a dataframe.
        :param result_by_locus: Results by locus.
        :return: DataFrame
        """
        data_out = pd.DataFrame(
            [{
                'locus': locus,
                'allele': res.allele_str,
            } for locus, res in result_by_locus.items()]
        )
        data_out['is_novel'] = data_out['allele'].str.startswith('n')
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
            allele_id = res_by_locus[locus].allele_str
            with open(dir_novel_alleles / f'{locus}-{allele_id}.fasta', 'w') as handle:
                SeqIO.write(SeqRecord(
                    id=f'{locus}_{allele_id}',
                    description=f"closest_match={','.join(res_by_locus[locus].allele_results[0].closest_alleles)}",
                    seq=Seq(res_by_locus[locus].allele_results[0].sequence),
                ), handle, 'fasta')

    def query(self, path_fasta: Path, out_json: Path, out_dir: Path, out_tsv: Path, threads: int) -> None:
        """
        Queries a scheme.
        :param path_fasta: Input FASTA path
        :param out_json: JSON output path
        :param out_dir: Output directory
        :param out_tsv: Output TSV path
        :param threads: Number of threads
        :return: None
        """
        # Query alleles
        allele_query = AlleleQueryMinimap2(
            dir_db=self._dir_db,
            dir_out=out_dir,
            multi_strategy=MultiStrategy(self._multi),
            min_id_novel=99,
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
        export_json(result_by_locus, out_json, profile, pct_match)
        if self._export_novel and data_results['is_novel'].sum() > 0:
            self._export_novel_allele_seqs(data_results, result_by_locus, out_dir)

        # Export TSV output
        if out_tsv is not None:
            data_results.to_csv(out_tsv, sep='\t', index=False)

def _parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parses the command line arguments.
    :param args: Optional arguments
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument('-f', '--fasta', type=Path, required=True, help='Input FASTA path')
    parser.add_argument('-d', '--db', type=Path, required=True, help='Database path')

    # Output
    parser.add_argument('-o', '--out-json', type=Path, help='JSON output file', default=Path('mist.json'))
    parser.add_argument('--out-tsv', type=Path, help='TSV output file')
    parser.add_argument('--out-dir', type=Path, help='Output directory')
    parser.add_argument('--export-novel', action='store_true', help='Create FASTA files for (potential) novel alleles')
    parser.add_argument('--keep-minimap2', action='store_true', help='Store the minimap2 output')

    # Various
    parser.add_argument('-t', '--threads', type=int, default=1, help='Nb. of threads to use')
    parser.add_argument('--min-id-novel', type=int, default=99, help='Minimum %% identity for novel alleles')
    parser.add_argument(
        '-m', '--multi', choices=[s.value for s in MultiStrategy], default=MultiStrategy.ALL.value,
        help='Strategy to handle multiple perfect hits')
    parser.add_argument('--loci', help='Limit to these loci')

    # Version
    parser.add_argument(
        '--version', help='Print version and exit', action='version', version=f'MiST {__version__}',)
    return parser.parse_args(args)


def export_json(results_by_locus: dict[str, model.QueryResult], path_out: Path, profile: model.Profile, pct_match: float) -> None:
    """
    Exports the results in JSON format.
    :param results_by_locus: Result(s) by locus
    :param path_out: Output path
    :param profile: Detected profile
    :param pct_match: Percent match for the profile
    :return: None
    """
    with open(path_out, 'w') as handle:
        json.dump({
            'alleles': {locus: dataclasses.asdict(res) for locus, res in results_by_locus.items()},
            'profile': {**dataclasses.asdict(profile), 'pct_match': pct_match} if profile is not None else None,
            'metadata': {
                'timestamp': datetime.now().isoformat(),
                'tool_version': __version__,
            }
        }, handle, indent=2, cls=CustomEncoder)
    logger.info(f'Output stored in: {path_out}')


def main(args_str: Optional[Sequence[str]] = None) -> None:
    """
    Runs the main script.
    :param args_str: Command line arguments (optional)
    :return: None
    """
    # Setup
    t0 = datetime.now()
    check_dependencies(['minimap2'])
    args = _parse_args(args_str)
    # Create the output directory before initializing logging
    if args.out_dir is not None:
        args.out_dir.mkdir(parents=True, exist_ok=True)
    initialize_logging(dir_logs=args.out_dir if args.out_dir else None)

    # Validate args
    if (args.export_novel or args.keep_minimap2) and (args.out_dir is None):
        raise ValueError("Output directory ('--out-dir') must be specified when exporting FASTA / Minimap2 files")

    # Perform the query
    query = MistQuery(
        dir_db=args.db,
        multi=args.multi,
        keep_minimap2=args.keep_minimap2,
        loci=args.loci,
        export_novel=args.export_novel,
    )
    query.query(
        path_fasta=args.fasta,
        out_json=args.out_json,
        out_dir=args.out_dir,
        out_tsv=args.out_tsv,
        threads=args.threads,
    )

    # Typing completed
    logger.info("Make sure to cite the corresponding database when using this in your research")
    path_citation = files('mist').joinpath('resources/citation.txt')
    logger.info(f'Please cite: {path_citation.read_text()}')
    logger.info(f"Processing time: {(datetime.now() - t0).total_seconds():.2f} seconds")

if __name__ == '__main__':
    main()
