import unittest

from click.testing import CliRunner

from mist.app.loggers.logger import initialize_logging
from mist.scripts.cli import cli


class TestList(unittest.TestCase):
    """
    Tests for the listing functionality.
    """

    @staticmethod
    def _parse_stdout(stdout: str) -> list[tuple[str, str]]:
        """
        Parses the schemes from the stdout.
        :param stdout: Tool stdout
        :return: List of schemes and URLs
        """
        schemes = []
        for row in stdout.splitlines():
            if row.startswith('----'):
                break
            scheme, second = row.rsplit(maxsplit=1)
            schemes.append((scheme, second))
        return schemes

    def test_list_cgmlst_org(self) -> None:
        """
        Tests listing the schemes from cgMLST.org.
        :return: None
        """
        runner = CliRunner()
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'cgmlstorg'
            ], catch_exceptions=False
        )
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 20, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)

    def test_list_enterobase(self) -> None:
        """
        Tests listing the schemes from EnteroBase.
        :return: None
        """
        runner = CliRunner()
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'enterobase'
            ], catch_exceptions=False
        )
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 10, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)

    def test_list_bigsdb_pubmlst_species(self) -> None:
        """
        Tests listing the species from PubMLST.
        :return: None
        """
        runner = CliRunner()
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pubmlst'
            ], catch_exceptions=False
        )
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 50, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)

    def test_list_bigsdb_pasteur_species(self) -> None:
        """
        Tests listing the species from BIGSdb Pasteur.
        :return: None
        """
        runner = CliRunner()
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pubmlst'
            ], catch_exceptions=False
        )
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 10, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)

    def test_list_bigsdb_pubmlst_schemes(self) -> None:
        """
        Tests listing the schemes from BIGSdb from a particular species.
        :return: None
        """
        runner = CliRunner()

        # Retrieve the species
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pubmlst',
            ], catch_exceptions=False)
        species = TestList._parse_stdout(result.output)

        # Retrieve the DB name for Neisseria
        neisseria_db_name = next(url for name, url in species if 'neisseria' in name.lower())

        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,
            [
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pubmlst',
                '--db', neisseria_db_name
            ], catch_exceptions=False)
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 20, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)

    def test_list_bigsdb_pasteur_schemes(self) -> None:
        """
        Tests listing the schemes from BIGSdb from a particular species.
        :return: None
        """
        runner = CliRunner()

        # Retrieve the species
        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,[
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pasteur',
            ], catch_exceptions=False)
        species = TestList._parse_stdout(result.output)

        # Retrieve the DB name for Neisseria
        db_name = next(url for name, url in species if 'klebsiella' in name.lower())

        # noinspection PyTypeChecker
        result = runner.invoke(
            cli,
            [
                'list',
                '--downloader', 'bigsdb',
                '--host', 'pasteur',
                '--db', db_name
            ], catch_exceptions=False)
        schemes = TestList._parse_stdout(result.output)
        self.assertGreater(len(schemes), 20, "Fewer schemes listed than expected")
        self.assertTrue(result.exit_code == 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
