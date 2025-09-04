import shutil
from importlib.resources import files
from pathlib import Path
from typing import Any

from furl import furl

from mist.app.dbs.basedownloader import BaseDownloader
from mist.app.loggers.logger import logger
from mist.app.utils import restutils, sequenceutils
from mist.app.utils.command import Command


class BIGSDbAuthDownloader(BaseDownloader):
    """
    Downloader for BIGSdb databases with authentication.
    """

    def __init__(self, **kwargs: Any) -> None:
        """
        Initializes the downloader.
        :param kwargs: Keyword arguments
        :return: None
        """
        super().__init__()
        try:
            self._dir_tokens: Path = Path(kwargs['dir_tokens'])
            self._key_name: str = kwargs['key_name']
            self._site: str = kwargs['site']
        except KeyError as err:
            raise RuntimeError(f"Required argument {err} missing")
        if self._dir_tokens is None:
            raise ValueError("Token directory should be set")

    @property
    def path_script(self) -> str:
        """
        Returns the path to the downloader script
        :return: Path to the downloader script
        """
        return str(files('mist').joinpath('resources/pubmlst/download_pubmlst.py'))

    def _check_credentials(self, dir_keys: Path) -> None:
        """
        Checks if the required credential are present.
        :param dir_keys: Path to the directory containing the credentials.
        :return: None
        """
        if not dir_keys.exists():
            raise FileNotFoundError(f'Directory does not exist: {dir_keys}')

    def _retrieve_page(self, url: str, path_out: Path) -> Path:
        """
        Downloads a FASTA file using OAuth.
        :param url: URL to the FASTA file
        :param path_out: Output path
        :return: None
        """
        command = Command(' '.join([
            'python', self.path_script,
            '--key_name', self._key_name,
            '--site', self._site,
            '--url', url,
            '--token_dir', str(self._dir_tokens.absolute()),
            '--cron',
            f'> {path_out.absolute()}'
        ]))
        command.run(self.dir_out, timeout=60)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error retrieving: {url}\n{command.stderr}')
        return path_out

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
            path_fasta = self._retrieve_page(url_fasta, self.dir_out / f'{locus}.fasta')
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
        # Download the scheme metadata
        metadata = restutils.retrieve_page_data(url).json()
        logger.info(f"Metadata received: {metadata['locus_count']:,} loci")

        # Create a temporary directory to store the tokens (to avoid overwriting the original files)
        dir_tokens_temp = dir_out / 'tokens'
        dir_tokens_temp.mkdir(exist_ok=True)
        for p_file in self._dir_tokens.iterdir():
            if p_file.is_file():
                shutil.copyfile(p_file, dir_tokens_temp / p_file.name)
        logger.info(f"Tokens copied to: {dir_tokens_temp}")
        self._dir_tokens = dir_tokens_temp

        # Download FASTA files
        paths_fasta = self._download_fasta_files(metadata)

        # Download profiles
        if include_profiles:
            href_profiles = str(furl(url).add(path='profiles_csv'))
            self._retrieve_page(href_profiles, self.dir_out / 'profiles.tsv')
            logger.info('Profiles downloaded')

        # Create a TXT file with all FASTA files
        self.create_fasta_list(paths_fasta)
        logger.info('Scheme downloaded successfully')
