import click
from furl import furl

from mist.app.dbs import BASE_URL_BY_DB, DOWNLOADERS


class MistList:
    """
    Commands to list the available schemes that can be downloaded.
    """

    def __init__(self, downloader: str, host: str, db: str | None = None) -> None:
        """
        Initializes the class.
        :param downloader: Key of the downloader class
        :param host: Host name
        :param db: Database name
        :return: None
        """
        download_key = downloader if downloader != 'bigsdb_auth' else 'bigsdb'
        self._downloader = DOWNLOADERS[download_key]()
        self._downloader_key = downloader
        self._host = host if host is not None else 'base'
        self._db = db

    def _build_bigsdb_url(self) -> str:
        """
        Returns the BIGSdb URL.
        :return: BIGSdb URL
        """
        if self._db is None:
            return BASE_URL_BY_DB['bigsdb'][self._host]
        else:
            api_root = BASE_URL_BY_DB['bigsdb'][self._host]
            return str(furl(api_root).add(path='db').add(path=self._db).add(path='schemes'))

    def print_available_schemes(self) -> None:
        """
        Prints the available schemes and the corresponding URL.
        :return: None
        """
        # Determine the base URL
        if self._downloader_key in ('bigsdb', 'bigsdb_auth'):
            base_url = furl(self._build_bigsdb_url())
        else:
            base_url = furl(BASE_URL_BY_DB[self._downloader_key][self._host])

        # Retrieve the available schemes
        schemes = self._downloader.get_available_schemes(furl(base_url), species=self._db is None)
        if len(schemes) ==0:
            click.echo('No schemes found')
            raise RuntimeError('No schemes found')

        # Log the results
        for scheme_name, url in schemes:
            click.echo(f"{scheme_name:<30} {url}")
        click.echo('-' * 10)

        if (self._db is None) and ('bigsdb' in self._downloader_key):
            click.echo(
                "IMPORTANT: Due to the large number of schemes hosted on PubMLST, they are grouped by species.\n"
                "To list the available databases for a species, use:\n"
                f"mist list --downloader {self._downloader_key} --host {self._host} --schemes --db {{DB}}\n"
                "where {DB} is the name of the database from the list above (with the '_seqdef' suffix)"
            )
        else:
            click.echo(
                'Download the scheme using:\n'
                f'mist download --downloader {self._downloader_key} --url {{URL}} --output {{DIR}}')
