import unittest
from importlib.resources import files
from pathlib import Path

from mist.app.utils import testingutils
from mist.scripts.mlst_phylogeny import main


class TestPhylogeny(unittest.TestCase):
    """
    Tests for the phylogeny script.
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

    def test_tsv_input(self) -> None:
        """
        Tests the MLST phylogeny script with TSV input.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            main([
                '--tsv', *(str(x) for x in TestPhylogeny.get_output_files('.tsv')),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])
            self.assertTrue(path_out_dists.exists())
            self.assertTrue(path_out_matrix.exists())

    def test_json_input(self) -> None:
        """
        Tests the MLST phylogeny script with JSON input.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            main([
                '--json', *(str(x) for x in TestPhylogeny.get_output_files('.json')),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])
            self.assertTrue(path_out_dists.exists())
            self.assertTrue(path_out_matrix.exists())

    def test_mixed_input(self) -> None:
        """
        Tests the MLST phylogeny script with JSON input.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            path_out_dists = Path(dir_temp, 'output_dists.tsv')
            path_out_matrix = Path(dir_temp, 'output_matrix.tsv')
            main([
                '--json', *(str(x) for x in TestPhylogeny.get_output_files('.json')[:2]),
                '--tsv', *(str(x) for x in TestPhylogeny.get_output_files('.tsv')[2:]),
                '--out-dists', str(path_out_dists),
                '--out-matrix', str(path_out_matrix)
            ])
            self.assertTrue(path_out_dists.exists())
            self.assertTrue(path_out_matrix.exists())
