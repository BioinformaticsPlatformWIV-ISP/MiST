from pathlib import Path
from typing import Any

from furl import furl

from mist.app.dbs.basedownloader import BaseDownloader
from mist.app.loggers.logger import logger
from mist.app.utils import restutils, sequenceutils


class BIGSDbDownloader(BaseDownloader):
    """
    Downloader for BIGSdb databases.
    Note: this downloader does **NOT** handle authorization and will not download all available data.
    """

    DOWNLOADER_KEY = 'bigsdb'

    def _download_fasta_files(self, metadata: dict[str, Any]) -> list[Path]:
        """
        Downloads the scheme FASTA files.
        :param metadata: Scheme metadata
        :return: List of FASTA paths
        """
        paths_fasta = []
        for href_locus in metadata['loci']:
            locus = furl(href_locus).path.segments[-1]
            url_fasta = str(furl(href_locus).add(path='alleles_fasta'))
            path_fasta = restutils.download_fasta(locus_name=locus, url=url_fasta, dir_out=self.dir_out)
            paths_fasta.append(path_fasta)
            logger.info(
                f'Downloaded {path_fasta.name} ({sequenceutils.count_sequences(path_fasta):,} sequences)'
            )
        return paths_fasta

    def _download(self, url: str, dir_out: Path, include_profiles: bool = False) -> None:
        """
        Downloads the target scheme.
        :param url: Scheme URL
        :param dir_out: Output directory
        :param include_profiles: Include profiles
        :return: None
        """
        logger.warning('Downloading from PubMLST **without** authentication, not all data will be retrieved!')

        # Download the scheme metadata
        metadata = restutils.retrieve_page_data(url).json()
        logger.info(f"Metadata received: {metadata['locus_count']:,} loci")

        # Download FASTA files
        paths_fasta = self._download_fasta_files(metadata)

        # Download profiles
        if include_profiles:
            href_profiles = str(furl(url).add(path='profiles_csv'))
            response = restutils.retrieve_page_data(href_profiles)
            with open(self.dir_out / 'profiles.tsv', 'w') as handle:
                handle.write(response.text)
            logger.info('Profiles downloaded')

        # Create a TXT file with all FASTA files
        self.create_fasta_list(paths_fasta)
        logger.info('Scheme downloaded successfully')

    def get_available_schemes(self, base_url: furl, **kwargs: Any) -> list[tuple[str, str]]:
        """
        Retrieve a list of available cgMLST schemes and their URL.
        :param base_url: Base URL
        :param kwargs: Keyword arguments
        :return: List of schemes / species and URL
        """
        data = restutils.retrieve_page_data(str(base_url))

        # Species
        schemes = []

        if kwargs.get('species', False) is True:
            for row in data.json():
                # Retrieve the available species
                try:
                    db_info = next(d for d in row['databases'] if d['name'].endswith('_seqdef'))
                except StopIteration:
                    logger.debug(f"No SeqDef database found for: {row['description']}")
                    continue
                schemes.append((row['description'].strip(), db_info['href'].split('/')[-1]))
        else:
            for data_scheme in data.json()['schemes']:
                schemes.append((data_scheme['description'], data_scheme['scheme']))
        return schemes
