#!/usr/bin/env python
import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

from mist.app.dbs import DOWNLOADERS
from mist.app.loggers.logger import initialize_logging


class DownloadScheme:
    """
    Wrapper for downloading schemes.
    """

    def __init__(self, args: argparse.Namespace) -> None:
        """
        Initializes the downloader.
        :param args: Command line arguments
        :return: None
        """
        self._args = args

    @staticmethod
    def parse_args(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Command line arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('-u', '--url', help='URL of the scheme', required=True)
        parser.add_argument('-o', '--output', type=Path, default=Path('mist_download'), help='Output directory')
        parser.add_argument('-p', '--include-profiles', action='store_true', help='Download the profiles')
        parser.add_argument('-d', '--downloader', choices=list(DOWNLOADERS.keys()), required=True, help='Downloader')
        # Authorization for BIGSdb
        parser.add_argument('--dir-tokens', type=Path, help='Directory with access tokens', default=Path.cwd() / '.bigsdb_tokens')
        parser.add_argument('--key-name', default='PubMLST', help='Key name')
        parser.add_argument('--site', default='PubMLST', help='Site')
        return parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the downloader.
        :return: None
        """
        downloader = DOWNLOADERS[self._args.downloader](
            dir_tokens=self._args.dir_tokens,
            key_name=self._args.key_name,
            site=self._args.site,
        )
        downloader.download_scheme(self._args.url, self._args.output, self._args.include_profiles)


def main() -> None:
    """
    Main script.
    :return: None
    """
    args = DownloadScheme.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    initialize_logging(args.output)
    download_scheme = DownloadScheme(args)
    download_scheme.run()


if __name__ == '__main__':
    main()
