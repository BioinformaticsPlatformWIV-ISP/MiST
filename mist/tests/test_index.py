import unittest
from importlib.resources import files
from pathlib import Path

from click.testing import CliRunner

from mist.app.loggers.logger import initialize_logging, logger
from mist.app.utils import dbutils, testingutils
from mist.scripts.cli import cli


class TestIndex(unittest.TestCase):
    """
    Tests for the indexing functionality.
    """

    def test_index_single_fasta(self) -> None:
        """
        Tests indexing on a single fasta file.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,[
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    '--output', str(dir_temp),
                    '--threads', '4'
                ], catch_exceptions=False
            )
            logger.info(result.output)
            self.assertTrue(result.exit_code == 0)
            self.assertTrue(dbutils.is_valid_db(Path(dir_temp)))

    def test_index_single_fasta_debug(self) -> None:
        """
        Tests indexing on a single fasta file in debug mode.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli, [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    '--output', str(dir_temp),
                    '--threads', '4',
                    '--debug',
                ], catch_exceptions=False
            )
            logger.info(result.output)
            self.assertTrue(dbutils.is_valid_db(Path(dir_temp)))
            self.assertEqual(result.exit_code, 0)

    def test_index_single_fasta_diff_fmt(self) -> None:
        """
        Tests indexing on a single fasta file in a different format.
        Each allele has just an identifier instead of the full allele name:
        >1
        [SEQ]
        >2
        [SEQ]
        ...
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli, [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-fmt.fasta')),
                    '--output', str(dir_temp),
                    '--threads', '4'
                ], catch_exceptions=False
            )
            logger.info(result.output)
            self.assertTrue(dbutils.is_valid_db(Path(dir_temp)))
            self.assertEqual(result.exit_code, 0)

    def test_index_from_list(self) -> None:
        """
        Tests indexing FASTA files from a list.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # Create a TXT file with the FASTA files
            path_txt = Path(dir_temp, 'fasta_in.txt')
            with path_txt.open('w') as handle:
                handle.write(str(files('mist').joinpath('resources/testdata/NEIS0140-fmt.fasta')) + '\n')
                handle.write(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')) + '\n')

            # Index the database
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli, [
                    'index',
                    '--fasta-list', str(path_txt),
                    '--output', str(dir_temp),
                    '--threads', '4',
                ], catch_exceptions=False
            )
            logger.info(result.output)
            self.assertTrue(dbutils.is_valid_db(Path(dir_temp)))
            self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
