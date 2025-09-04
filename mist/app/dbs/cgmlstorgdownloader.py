import time
import zipfile
from pathlib import Path
from typing import Any

import requests
from furl import furl

from mist.app.dbs.basedownloader import BaseDownloader
from mist.app.loggers.logger import logger
from mist.app.utils.restutils import retrieve_page_data


class CgMLSTOrgDownloader(BaseDownloader):
    """
    Downloader for cgMLST.org databases.
    """

    @staticmethod
    def retrieve_page_data(
        url: str, retries: int = 3, timeout: int = 60 * 5
    ) -> requests.Response:
        """
        Retrieves data from the given URL.
        :param url: URL
        :param retries: Number of retries
        :param timeout: Timeout (in seconds)
        :return: URL data
        """
        error = None
        for retry in range(retries):
            try:
                response = requests.get(url, timeout=timeout)
                response.raise_for_status()
                return response
            except BaseException as err:
                logger.warning(
                    f'Error retrieving page data ({err}), retrying (attempt {retry + 1})'
                )
                error = err
                time.sleep(1)
        raise RuntimeError(f'Error retrieving page data: {url}: {error}')

    @staticmethod
    def retrieve_fasta_archive(url: str, output_dir: Path) -> Path:
        """
        Retrieves the FASTA archive from the cgmlst.org website
        :param url: URL
        :param output_dir: Output directory
        :return: Path to the output FASTA
        """
        output_path = output_dir / 'fasta.zip'
        with open(output_path, 'wb') as handle:
            handle.write(CgMLSTOrgDownloader.retrieve_page_data(url).content)
        logger.info(f'Archive downloaded to: {output_path}')
        return output_path

    @staticmethod
    def extract_fasta_files(input_archive: Path, output_dir: Path) -> list[Path]:
        """
        Extracts the FASTA files from the zip archive
        :param input_archive: Downloaded zip archive
        :param output_dir: Output directory
        :return: List of paths to the FASTA files
        """
        with zipfile.ZipFile(input_archive, 'r') as handle:
            if handle.testzip() is not None:
                raise zipfile.BadZipFile(f"Corrupt zip file: {input_archive}")
            handle.extractall(output_dir)
        return [f for f in output_dir.iterdir() if f.name.endswith('.fasta')]

    def _download(self, url: str, dir_out: Path, include_profiles: bool = False, **kwargs: Any) -> None:
        """
        Downloads the target scheme.
        :param url: Scheme URL
        :param dir_out: Output directory
        :param include_profiles: Include profiles
        :param kwargs: Additional arguments
        :return: None
        """
        # Download the archive containing the FASTA files
        href_alleles = str(furl(url).add(path='alleles'))
        path_fasta_arch = CgMLSTOrgDownloader.retrieve_fasta_archive(href_alleles, dir_out)
        paths_fasta = CgMLSTOrgDownloader.extract_fasta_files(path_fasta_arch, dir_out)
        logger.info(f"Found {len(paths_fasta):,} FASTA files")

        # Download the profiles (if specified)
        if include_profiles:
            href_profiles = str(furl(url).add(path='locus/?content-type=csv'))
            with open(self.dir_out / 'profiles.tsv', 'w') as handle:
                handle.write(retrieve_page_data(href_profiles).text)
            logger.info('Profiles downloaded')

        # Create a TXT file with all FASTA files
        self.create_fasta_list(paths_fasta)
        logger.info('Scheme downloaded successfully')
