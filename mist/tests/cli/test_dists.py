import unittest
from importlib.resources import files
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from mist.app.utils import testingutils
from mist.scripts.cli import cli


class TestDists(unittest.TestCase):
    """
    Tests for the mist dists script.
    """

    @staticmethod
    def get_output_files(ext: str) -> list[Path]:
        """
        Returns the output files for testing with the given extension.
        :param ext: Extension of the output files.
        :return: List of files
        """
        dir_test = Path(str(files('mist').joinpath('resources/testdata/output')))
        return sorted(list(dir_test.glob(f'*{ext}')))

    def _is_valid_output_file(self, path: Path) -> bool:
        """
        Checks if the (tabular) output file is valid.
        :param path: Output file
        :return: True if valid
        """
        df_matrix = pd.read_table(path)
        self.assertFalse(df_matrix.empty)
        return True

    def test_tsv_input(self) -> None:
        """
        Tests the mist dists script with TSV input.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            # noinspection PyTypeChecker
            result = runner.invoke(cli, [
                'dists',
                *(str(x) for x in TestDists.get_output_files('.tsv')),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ], catch_exceptions=False)

            # Verify output
            self.assertTrue(self._is_valid_output_file(path_out_dists))
            self.assertTrue(self._is_valid_output_file(path_out_matrix))
            self.assertEqual(result.exit_code, 0)

    def test_json_input(self) -> None:
        """
        Tests the mist dists script with JSON input.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            # noinspection PyTypeChecker
            result = runner.invoke(cli, [
                'dists',
                *(str(x) for x in TestDists.get_output_files('.json')),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])

            # Verify output
            self.assertTrue(self._is_valid_output_file(path_out_dists))
            self.assertTrue(self._is_valid_output_file(path_out_matrix))
            self.assertEqual(result.exit_code, 0)

    def test_mixed_input(self) -> None:
        """
        Tests the mist dists script with JSON input.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            # noinspection PyTypeChecker
            result = runner.invoke(cli, [
                'dists',
                *(str(x) for x in TestDists.get_output_files('.json')[:2]),
                *(str(x) for x in TestDists.get_output_files('.tsv')[2:]),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])

            # Verify output
            self.assertTrue(self._is_valid_output_file(path_out_dists))
            self.assertTrue(self._is_valid_output_file(path_out_matrix))
            self.assertEqual(result.exit_code, 0)

    def test_invalid_input(self) -> None:
        """
        Tests the mist dists script with invalid input (too few input files).
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            # noinspection PyTypeChecker
            result = runner.invoke(cli, [
                'dists',
                next(str(x) for x in TestDists.get_output_files('.json')),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])
            self.assertNotEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()
