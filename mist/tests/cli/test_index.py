import json
import shutil
import unittest
from importlib.resources import files
from pathlib import Path
from unittest.mock import Mock

from click.testing import CliRunner

from mist.app import model
from mist.app.loggers.logger import initialize_logging, logger
from mist.app.query.profileindex import ProfileIndex
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
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                ],
                catch_exceptions=False,
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
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                    '--debug',
                ],
                catch_exceptions=False,
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
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-fmt.fasta')),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                ],
                catch_exceptions=False,
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
                cli,
                [
                    'index',
                    '--fasta-list',
                    str(path_txt),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                ],
                catch_exceptions=False,
            )
            logger.info(result.output)
            self.assertTrue(dbutils.is_valid_db(Path(dir_temp)))
            self.assertEqual(result.exit_code, 0)

    def test_index_copies_db_info_from_profiles_dir(self) -> None:
        """
        Tests that db_info.json next to --profiles (but not next to the FASTA files) still gets copied to the
        output directory.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp, testingutils.get_temp_dir() as dir_profiles:
            dir_temp, dir_profiles = Path(dir_temp), Path(dir_profiles)

            # Place profiles.tsv and db_info.json together, separate from the FASTA files
            path_profiles_test = str(files('mist').joinpath('resources/testdata/profiles.tsv'))
            shutil.copyfile(path_profiles_test, dir_profiles / 'profiles.tsv')
            db_info = {'url': 'https://example.org/scheme', 'downloader': 'bigsdb', 'download_date': '2026-08-10'}
            with (dir_profiles / 'db_info.json').open('w') as handle:
                json.dump(db_info, handle)

            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    str(files('mist').joinpath('resources/testdata/NEIS0159-subset.fasta')),
                    '--profiles',
                    str(dir_profiles / 'profiles.tsv'),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                ],
                catch_exceptions=False,
            )
            logger.info(result.output)
            self.assertEqual(result.exit_code, 0)

            path_db_info_out = dir_temp / 'db_info.json'
            self.assertTrue(path_db_info_out.exists())
            with path_db_info_out.open() as handle:
                self.assertEqual(db_info, json.load(handle))

    def test_index_without_build_profile_index_skips_it(self) -> None:
        """
        Tests that no profile index is built unless --build-profile-index is passed, even with --profiles set.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    str(files('mist').joinpath('resources/testdata/NEIS0159-subset.fasta')),
                    '--profiles',
                    str(files('mist').joinpath('resources/testdata/profiles.tsv')),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                ],
                catch_exceptions=False,
            )
            logger.info(result.output)
            self.assertEqual(result.exit_code, 0)
            self.assertFalse((Path(dir_temp) / 'profile_index').exists())

    def test_index_with_build_profile_index_builds_queryable_index(self) -> None:
        """
        Tests that --build-profile-index builds a profile index that returns the correct matching ST.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,
                [
                    'index',
                    str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')),
                    str(files('mist').joinpath('resources/testdata/NEIS0159-subset.fasta')),
                    '--profiles',
                    str(files('mist').joinpath('resources/testdata/profiles.tsv')),
                    '--output',
                    str(dir_temp),
                    '--threads',
                    '4',
                    '--build-profile-index',
                ],
                catch_exceptions=False,
            )
            logger.info(result.output)
            self.assertEqual(result.exit_code, 0)

            dir_profile_index = Path(dir_temp) / 'profile_index'
            self.assertTrue(dir_profile_index.exists())

            # ST1 in profiles.tsv is NEIS0140-subset=10, NEIS0159-subset=67
            result_locus1 = Mock(spec=model.QueryResult)
            result_locus1.allele_str = '10'
            result_locus1.allele_results = [Mock()]
            result_locus2 = Mock(spec=model.QueryResult)
            result_locus2.allele_str = '67'
            result_locus2.allele_results = [Mock()]

            profile_idx = ProfileIndex(dir_index=dir_profile_index)
            profiles, nb_matches = profile_idx.query(
                {'NEIS0140-subset': result_locus1, 'NEIS0159-subset': result_locus2}
            )
            self.assertEqual(len(profiles), 1)
            self.assertEqual(profiles[0].name, 'ST1')
            self.assertEqual(nb_matches, 2)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
