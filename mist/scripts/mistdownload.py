from pathlib import Path

from mist.app.dbs import DOWNLOADERS


class MistDownload:
    """
    Wrapper for downloading schemes.
    """

    def __init__(self, url: str, output: Path, include_profiles: bool, downloader: str, dir_tokens: Path, key_name: str,
                 site: str) -> None:
        """
        Initialize the scheme download class.
        :param url: Scheme URL
        :param output: Output directory
        :param include_profiles: If True, profiles are kept
        :param downloader: Key of the downloader class
        :param dir_tokens: Token directory
        :param key_name: Key name
        :param site: Site name
        :return: None
        """
        self._url = url
        self._output = output
        self._include_profiles = include_profiles
        self._downloader = downloader
        self._dir_tokens = dir_tokens
        self._key_name = key_name
        self._site = site

    def run(self) -> None:
        """
        Runs the downloader.
        :return: None
        """
        self._output.mkdir(exist_ok=True, parents=True)
        downloader = DOWNLOADERS[self._downloader](
            dir_tokens=self._dir_tokens,
            key_name=self._key_name,
            site=self._site,
        )
        downloader.download_scheme(self._url, self._output, self._include_profiles)
