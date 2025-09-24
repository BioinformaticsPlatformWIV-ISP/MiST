import gzip
import re
import shutil
from pathlib import Path
from typing import Any

import bs4
import requests
from furl import furl

from mist.app.dbs.basedownloader import BaseDownloader
from mist.app.loggers.logger import logger


class EnteroBaseDownloader(BaseDownloader):
    """
    Downloads (cg)MLST schemes from EnteroBase.
    """

    DOWNLOADER_KEY = 'enterobase'
    STATUS_CODE_OK = 200

    @staticmethod
    def extract_fasta_gz_urls(url: furl) -> list[furl]:
        """
        Extracts all `.fasta.gz` file URLs from an HTML directory listing.
        :param url: Base directory URL (as `furl`)
        :return: List of full URLs (as `furl`) to `.fasta.gz` files
        """
        response = requests.get(url.tostr())
        if response.status_code != EnteroBaseDownloader.STATUS_CODE_OK:
            raise RuntimeError(f'Error accessing: {url} (status code: {response.status_code})')
        content = bs4.BeautifulSoup(response.text, 'html.parser')
        links = content.find_all('a')
        return [
            url.copy().add(path=a['href'])
            for a in links
            if a['href'].endswith('.fasta.gz')
        ]

    @staticmethod
    def download_fasta_file(dir_out: Path, url: furl) -> Path:
        """
        Downloads the given FASTA file.
        :param dir_out: Directory to store the FASTA file
        :param url: URL to the FASTA file
        :return: Path to FASTA file
        """
        logger.debug(f"Downloading FASTA file: {url.path.segments[-1]}")
        locus_name = re.search(r'(.+)\.fasta\.gz$', url.path.segments[-1]).group(1)
        path_out = dir_out / f'{locus_name}.fasta'

        # Download and decompress
        response = requests.get(url.tostr(), stream=True)
        if response.status_code != EnteroBaseDownloader.STATUS_CODE_OK:
            raise RuntimeError(f'Error accessing: {url} (status code: {response.status_code})')
        with (
            gzip.GzipFile(fileobj=response.raw) as gzip_file,
            open(path_out, 'wb') as handle_out,
        ):
            shutil.copyfileobj(gzip_file, handle_out)
        logger.debug(f'FASTA file downloaded: {path_out.name}')
        return path_out

    @staticmethod
    def download_profiles(dir_out: Path, url: furl) -> Path:
        """
        Downloads the profiles file.
        :param dir_out: Output directory
        :param url: Scheme url
        :return: Path to profiles file
        """
        href_profiles = url.copy().add(path='profiles.list.gz')
        path_out = dir_out / 'profiles.tsv'

        # Download and decompress
        response = requests.get(href_profiles.tostr(), stream=True)
        with (
            gzip.GzipFile(fileobj=response.raw) as gzip_file,
            open(path_out, 'wb') as handle_out,
        ):
            shutil.copyfileobj(gzip_file, handle_out)
        logger.info(f'Profiles file downloaded: {path_out.name}')
        return path_out

    def _download(self, url: str, dir_out: Path, include_profiles: bool = False) -> None:
        """
        Downloads the target scheme.
        :param url: Scheme URL
        :param dir_out: Output directory
        :param include_profiles: Include profiles
        :return: None
        """
        # Get the URLs to the FASTA files
        urls_fasta = EnteroBaseDownloader.extract_fasta_gz_urls(furl(url))
        logger.info(f'Found {len(urls_fasta):,} FASTA files')
        if len(urls_fasta) == 0:
            raise RuntimeError('No FASTA file URLs extracted')

        # Download the FASTA files
        paths_fasta = []
        for url_fasta in urls_fasta:
            paths_fasta.append(EnteroBaseDownloader.download_fasta_file(dir_out, url_fasta))
        self.create_fasta_list(paths_fasta)

        # Download the profiles file
        if include_profiles:
            EnteroBaseDownloader.download_profiles(dir_out, furl(url))
        logger.info(f'Scheme downloaded to: {dir_out}')

    def get_available_schemes(self, base_url: furl, **kwargs: Any) -> list[tuple[str, str]]:
        """
        Retrieve a list of available cgMLST schemes and their URL.
        :param base_url: Base URL
        :param kwargs: Keyword arguments
        :return: None
        """
        response = requests.get(str(base_url))
        soup = bs4.BeautifulSoup(response.text, 'html.parser')
        schemes = []
        for a in soup.select('pre a'):
            href = a.get('href')
            name = a.get_text(strip=True).rstrip('/')

            # Skip parent directory entry
            if href == '../':
                continue
            url = base_url.copy().add(path=href)
            schemes.append((name, str(url)))
        return schemes
