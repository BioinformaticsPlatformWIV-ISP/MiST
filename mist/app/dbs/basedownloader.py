import abc
from pathlib import Path
from typing import Any, Optional

from mist.app.loggers.logger import logger


class BaseDownloader(metaclass=abc.ABCMeta):
    """
    Baseclass for scheme downloaders.
    """

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
