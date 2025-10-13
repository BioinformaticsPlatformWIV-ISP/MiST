import abc
import datetime
import json
from pathlib import Path
from typing import Any, ClassVar, Optional

from furl import furl

from mist.app import NAME_DB_INFO
from mist.app.loggers.logger import logger


class BaseDownloader(metaclass=abc.ABCMeta):
    """
    Baseclass for scheme downloaders.
    """
    DOWNLOADER_KEY: ClassVar[str]

    def __init__(self, **kwargs: Any) -> None:
        """
        Initializes the downloader.
        :param kwargs: Keyword arguments
        :return: None
        """
        self.dir_out: Optional[Path] = None

    @abc.abstractmethod
    def _download(self, url: str, dir_out: Path, include_profiles: bool = False) -> None:
        """
        Downloads the target scheme.
        :param url: Scheme URL
        :param dir_out: Output directory
        :param include_profiles: Include profiles
        :return: None
        """
        pass

    def download_scheme(self, url: str, dir_out: Path, include_profiles: bool = False) -> None:
        """
        Downloads the target scheme.
        :param url: Scheme URL
        :param dir_out: Output directory
        :param include_profiles: Include profiles
        :return: None
        """
        self.dir_out = dir_out
        self._download(url, dir_out, include_profiles)
        self.export_metadata_file(url)
        logger.info(
            f'You can create the index using:\n'
            f"mist index --fasta-list {dir_out / 'fasta_list.txt'} --output DB_NAME --threads 4")

    def create_fasta_list(self, paths_fasta: list[Path]) -> None:
        """
        Creates a TXT file with the FASTA file paths.
        :param paths_fasta: List of FASTA paths
        :return: None
        """
        path_fasta_list = self.dir_out / 'fasta_list.txt'
        with open(path_fasta_list, 'w') as handle:
            for path_fasta in paths_fasta:
                handle.write(str(path_fasta.absolute()))
                handle.write('\n')
        logger.info(f'FASTA list created: {path_fasta_list}')

    def export_metadata_file(self, url: str) -> None:
        """
        Exports the database metadata file.
        :param url: Scheme url
        :return: None
        """
        path_out = self.dir_out / NAME_DB_INFO
        with path_out.open('w') as handle:
            json.dump({
                'url': url,
                'downloader': self.DOWNLOADER_KEY,
                'download_date': datetime.datetime.now().isoformat(),
            }, handle, indent=2)
        logger.debug(f'DB info exported to: {path_out}')

    def get_available_schemes(self, base_url: furl, **kwargs: Any) -> list[tuple[str, str]]:
        """
        Retrieves the available schemes that can be downloaded.
        :param base_url: Base URL
        :param kwargs: Keyword arguments
        :return: List of available schemes
        """
        return []
