import concurrent.futures
import json
import logging
import re
import shutil
from pathlib import Path
from typing import Any

from Bio import SeqIO

from mist.app import NAME_DB_INFO
from mist.app.loggers.logger import logger
from mist.app.utils import (
    clusterutils,
    dbutils,
    minimap2utils,
    sequenceutils,
)
from mist.app.utils.clustersplit import ClusterSplit
from mist.app.utils.dependencies import check_dependencies


class MistIndex:
    """
    Main class to index schemes using MiST.
    """

    def __init__(self, paths_fasta: list[Path], path_profiles: Path | None, cutoff: int = 80, debug: bool = False) -> None:
        """
        Initializes the index script.
        :param paths_fasta: List of input FASTA files
        :param path_profiles: Path to profiles file (optional)
        :param cutoff: Clustering % identity cutoff
        :param debug: Enable debug mode
        :return: None
        """
        check_dependencies(['CD-HIT', 'minimap2', 'nucmer'])
        self._paths_fasta = paths_fasta
        self._path_profiles = path_profiles
        self._loci = [p.name.replace('.fasta', '') for p in paths_fasta]
        self._cutoff = cutoff
        self._debug = debug

    @staticmethod
    def process_locus(path_fasta: Path, dir_out: Path, cutoff: int, debug: bool) -> dict[str, Any]:
        """
        Processes a single locus
        :param path_fasta: Input FASTA path
        :param dir_out: Output directory
        :param cutoff: Clustering cutoff
        :param debug: Enable debug mode
        :return: Database information
        """
        dir_out.mkdir(parents=True, exist_ok=True)

        # Convert the input FASTA file to standard format
        path_fasta_full = dir_out / path_fasta.name
        dbutils.reformat_fasta(path_fasta, path_fasta_full)

        # Parse the FASTA file
        seq_by_id = {}
        with open(path_fasta_full) as handle:
            for i, seq in enumerate(SeqIO.parse(handle, 'fasta')):
                seq_by_id[seq.id] = {'seq': seq, 'idx': i}
        logger.info(f'Parsed: {len(seq_by_id):,} sequences from {path_fasta.name}')

        # Create clustered FASTA
        path_fasta_clustered = clusterutils.cluster_fasta(path_fasta_full, dir_out, cutoff, threads=1, debug=debug)
        clusters = clusterutils.parse_cluster_from_file(path_fasta_clustered)

        # Calculate hashes
        orientation_by_seq_id = {
            member['seq_id']: member['ori'] if member['ori'] is not None else '+'
            for cluster in clusters for member in cluster['members']
        }
        data_hashes = sequenceutils.extract_hashes(path_fasta_full, orientation_by_seq_id)
        logging.info(f'Parsed {len(data_hashes)} hashes from {path_fasta.name}')

        # Output data
        data_out = {
            'allele_id_regex': r'_(\d+)$',
            'alleles': {
                re.search(sequenceutils.REGEX_ALLELE, id_).group(1): {
                    'idx': seq_info['idx'],
                    'length': len(seq_info['seq']),
                }
                for id_, seq_info in seq_by_id.items()
            },
            'clusters': {c['members'][0]['seq_id']: c for c in clusters},
            'cutoff': cutoff,
            'fasta_clustered': path_fasta_clustered.name,
            'fasta_full': path_fasta_full.name,
            'name': path_fasta.name.replace('.fasta', ''),
            'hashes': {
                h: allele_id
                for allele_id, h in data_hashes[['allele_id', 'hash']].itertuples(
                    index=False, name=None
                )
            },
            'nb_seqs': len(data_hashes)
        }

        with open(dir_out / 'mist_db.json', 'w') as handle:
            json.dump(data_out, handle, indent=2)

        # Split up misaligned clusters
        splitter = ClusterSplit(dir_out, debug=debug)
        splitter.run()

        # Return the DB information
        return data_out

    def _extract_repr_alleles(self, dir_out: Path, threads: int) -> None:
        """
        Extracts the representative alleles from all FASTA files.
        :param dir_out: Output directory
        :param threads: Nb. of threads
        :return: None
        """
        res_by_locus = {}
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                future_to_locus = {
                    executor.submit(
                        MistIndex.process_locus, **{
                            'path_fasta': path_fasta,
                            'dir_out': dir_out / locus,
                            'cutoff': self._cutoff,
                            'debug': self._debug
                        },
                    ): locus for locus, path_fasta in zip(self._loci, self._paths_fasta)
                }

                for future in concurrent.futures.as_completed(future_to_locus):
                    locus = future_to_locus[future]
                    try:
                        res_by_locus[locus] = future.result()
                    except BaseException as err:
                        logging.error(f"Locus '{locus}' generated exception: {err}")
                        raise RuntimeError(err)
        except KeyboardInterrupt as err:
            logger.warning("Keyboard interrupt received. Cancelling all pending tasks.")
            for future in future_to_locus:
                future.cancel()
            raise err

    def create_index(self, dir_out: Path, threads: int = 4) -> None:
        """
        Creates the index in the target output directory.
        :param dir_out: Output directory
        :param threads: Number of threads
        :return: None
        """
        logger.info(f'Creating index for {len(self._paths_fasta):,} FASTA files')

        # Create the output directory
        dir_out = dir_out.expanduser().resolve()
        dir_out.mkdir(parents=True, exist_ok=True)

        # Extract representative alleles by clustering
        self._extract_repr_alleles(dir_out, threads)
        logger.info("Clustering completed successfully")

        # Create a FASTA file with all loci combined
        path_fasta_out = dir_out / 'loci_repr.fasta'
        nb_seqs = sequenceutils.merge_fasta_files(
            paths_fasta=[dir_out / locus / f'{locus}-clustered.fasta' for locus in self._loci],
            path_out=path_fasta_out,
        )
        minimap2utils.create_index(path_fasta_out)
        logger.info(f'Combined FASTA file created: {path_fasta_out} ({nb_seqs:,} sequences)')

        # Create a TXT file with all loci
        self._copy_db_info(dir_out)
        with open(dir_out / 'loci.txt', 'w') as handle:
            for locus in self._loci:
                handle.write(locus + '\n')

        # Copy profiles file
        if self._path_profiles:
            shutil.copyfile(self._path_profiles, dir_out / 'profiles.tsv')
        logger.info('Indexing finished successfully')

    def _copy_db_info(self, dir_out: Path) -> None:
        """
        Copies the database metadata information (if available).
        :param dir_out: Output directory
        :return: None
        """
        path_metadata = self._paths_fasta[0].parent / NAME_DB_INFO
        if not path_metadata.exists():
            logger.debug('No db_info.json file found, not copying DB metadata')
            return
        shutil.copyfile(path_metadata, dir_out / path_metadata.name)
