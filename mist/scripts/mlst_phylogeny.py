#!/usr/bin/env python
import argparse
import itertools
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import pandas as pd

from mist.app.loggers.logger import initialize_logging, logger


class MLSTPhylogeny:
    """
    Class to generate phylogenies from sequence typing outputs.

    The expected input is a list of TSV with the following columns:
    locus (string), allele (string / int).

    For example:
    locus	allele	is_novel
    SAUR0001        1       False
    SAUR0002        2       False
    SAUR0003        1       False
    SAUR0004        n815d8c       True
    SAUR0005        1       False
    SAUR0006        16      False
    SAUR0007        1       False

    Missing alleles should be indicated with allele '-'.
    """

    MIN_NB_DATASETS = 3

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :return: None
        """
        self._args = MLSTPhylogeny.parse_args(args)

    @staticmethod
    def parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (optional) arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--tsv', nargs='+', type=Path)
        parser.add_argument('-j', '--json', nargs='+', type=Path)
        parser.add_argument('-o', '--out-matrix', type=Path, help='Filtered allele matrix (TSV)', default=Path('allele_matrix.tsv'))
        parser.add_argument(
            '-d', '--out-dists', type=Path, help='Pairwise distance matrix (TSV)', default=Path('distances.tsv'))
        parser.add_argument(
            '-l', '--min-perc-loci', type=int, default=90,
            help='Minimum percentage of loci that should be present in a dataset')
        parser.add_argument(
            '-s', '--min-perc-samples', type=int, default=90,
            help='Minimum percentage of datasets where loci should be present')
        return parser.parse_args(args)

    @staticmethod
    def calc_distance(alleles_a: pd.Series, alleles_b: pd.Series) -> int:
        """
        Calculate the pairwise allele distance.
        :param alleles_a: Alleles a
        :param alleles_b: Alleles b
        :return: Number of allelic differences
        """
        total_dist = 0
        for allele_a, allele_b in zip(alleles_a, alleles_b):
            if allele_a == '-' and allele_b == '-':
                continue
            else:
                total_dist += 0 if allele_a == allele_b else 1
        return total_dist

    @staticmethod
    def parse_tsv(path_in: Path) -> tuple[str, pd.Series]:
        """
        Parses the input TSV file.
        :param path_in: Input path
        :return: Parsed data
        """
        logger.debug(f'Parsing file: {path_in}')
        data_tsv = pd.read_table(path_in, dtype=str, index_col='locus')
        return path_in.name.replace('.tsv', ''), data_tsv['allele']

    @staticmethod
    def parse_json(path_in: Path) -> tuple[str, pd.Series]:
        """
        Parses the input JSON file.
        :param path_in: Input path
        :return: Parsed data
        """
        logger.debug(f'Parsing file: {path_in}')
        with path_in.open() as handle:
            data = json.load(handle)
        alleles = pd.Series({locus: row['allele_str'] for locus, row in data['alleles'].items()})
        return path_in.name.replace('.json', ''), alleles

    def _log_nb_perfect_hits(self, name, alleles: pd.Series) -> None:
        """
        Logs the number and percentage of perfect hits for the input dataset.
        :param name: Dataset name
        :param alleles: Allele calls
        :return: None
        """
        nb_perfect = len(alleles) - alleles.eq('-').sum()
        perc_perfect = 100 * nb_perfect / len(alleles)
        logger.info(f"{name}: {nb_perfect}/{len(alleles):,} ({perc_perfect:.2f}%) perfect hits")

    def parse_input_files(self) -> pd.DataFrame:
        """
        Parse the input files.
        :return: DataFrame with detected alleles (index = sample name)
        """
        # Check if there are enough datasets
        paths_tsv = self._args.tsv if self._args.tsv is not None else []
        paths_json = self._args.json if self._args.json is not None else []
        if len(paths_tsv) + len(paths_json) < MLSTPhylogeny.MIN_NB_DATASETS:
            raise ValueError(
                f'At least 3 input files are required (found: {len(paths_tsv) + len(paths_json)})')

        # Parse the input
        datasets = {}
        for path_tsv in paths_tsv:
            name, alleles = MLSTPhylogeny.parse_tsv(path_tsv)
            datasets[name] = alleles
            self._log_nb_perfect_hits(name, alleles)
        for path_json in paths_json:
            name, alleles = MLSTPhylogeny.parse_json(path_json)
            datasets[name] = alleles
            self._log_nb_perfect_hits(name, alleles)

        # Create merged DataFrame
        return pd.DataFrame(data=list(datasets.values()), index=list(datasets.keys()), dtype=str)

    def calculate_distance_matrix(self, allele_data_filtered: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates the pairwise distance matrix
        :param allele_data_filtered: Filtered allele matrix
        :return: Distance matrix
        """
        # Calculate pair-wise distances
        distance_by_dataset_pair = {}
        for dataset_a, dataset_b in itertools.combinations(allele_data_filtered.index, r=2):
            key = tuple(sorted([dataset_a, dataset_b]))
            dist = MLSTPhylogeny.calc_distance(allele_data_filtered.loc[dataset_a], allele_data_filtered.loc[dataset_b])
            distance_by_dataset_pair[key] = dist

        # Create data frame with pairwise distances
        records_out = []
        for dataset_a in allele_data_filtered.index:
            records_out.append({
                dataset_b: distance_by_dataset_pair.get(tuple(sorted([dataset_a, dataset_b])), 0) for
                dataset_b in allele_data_filtered.index
            })
        return pd.DataFrame(records_out, index=allele_data_filtered.index)

    def filter_allele_matrix(self, allele_data: pd.DataFrame) -> pd.DataFrame:
        """
        Filters the allele matrix by removing:
        - Datasets with less than x% of loci detected
        - Loci present in less than x% of datasets
        :param allele_data: Allele data
        :return: Filtered allele data, loci cutoff, datasets cutoff
        """
        # Filter allele matrix (nb. of loci detected per dataset)
        nb_loci_detected = allele_data.apply(lambda x: len(x) - list(x).count('-'), axis=1)
        cutoff_loci = int(self._args.min_perc_loci * len(allele_data.columns) / 100)
        logger.info(f"Removing datasets with <{cutoff_loci} ({self._args.min_perc_loci}%) loci detected")
        allele_data_filt = allele_data[nb_loci_detected > cutoff_loci]
        logger.info(f"{len(allele_data_filt)}/{len(allele_data)} datasets passed filtering")

        # Filter allele matrix (loci detected in nb. of datasets)
        cutoff_datasets = int(self._args.min_perc_samples * len(allele_data_filt) / 100)
        logger.info(f"Removing loci detected in <{cutoff_datasets} ({self._args.min_perc_samples}%) samples")
        # noinspection PyUnresolvedReferences
        locus_present_in_datasets = (allele_data_filt != '-').sum(axis=0)
        allele_data_filt = allele_data_filt.loc[:, locus_present_in_datasets > cutoff_datasets]
        logger.info(f"{len(allele_data_filt.columns)}/{len(allele_data.columns)} loci passed filtering")
        return allele_data_filt

    def run(self) -> None:
        """
        The main method to construct the MLST phylogeny.
        :return: None
        """
        # Parse the input files
        df_alleles = self.parse_input_files()

        # Filter the allele matrix
        df_alleles_filt = self.filter_allele_matrix(df_alleles)
        df_alleles_filt.to_csv(self._args.out_matrix, index_label='ID', sep='\t')
        logger.info(f'Allele matrix exported to: {self._args.out_matrix}')

        # Calculate the distance matrix
        df_dists = self.calculate_distance_matrix(df_alleles_filt)

        # Check the distance matrix
        if all(max(values) == 0 for _, values in df_dists.iterrows()):
            raise ValueError('Empty distance matrix')

        # Export the distance matrix
        df_dists.to_csv(self._args.out_dists, index_label='ID', sep='\t')
        logger.info(f'Distance matrix exported to: {self._args.out_dists}')

        # Log the command for GrapeTree
        logger.info(f'You can construct a phylogeny using: grapetree --profile {self._args.out_matrix.name} --method MSTreeV2')


def main(args_str: Optional[Sequence[str]] = None) -> None:
    """
    Main script to run the MLST phylogeny.
    :param args_str: Command line arguments (optional)
    :return: None
    """
    initialize_logging()
    mlst_phylo = MLSTPhylogeny(args_str)
    mlst_phylo.run()


if __name__ == '__main__':
    main()
