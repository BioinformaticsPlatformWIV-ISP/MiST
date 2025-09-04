import unittest
from importlib.resources import files
from pathlib import Path

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import clusterutils, testingutils


class TestClusteringUtils(unittest.TestCase):
    """
    Tests for the clustering utilities.
    """

    def test_clustering_with_cdhit(self) -> None:
        """
        Tests the CD-HIT clustering function.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            path_in = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
            clusterutils.cluster_fasta(path_in, Path(dir_temp), 90)

    def test_clustering_with_cdhit_parsing(self) -> None:
        """
        Tests the CD-HIT clustering function with parsing.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            path_in = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
            path_fasta = clusterutils.cluster_fasta(path_in, Path(dir_temp), 90)
            clusters = clusterutils.parse_cluster_from_file(path_fasta)
            self.assertGreater(len(clusters), 0)
            self.assertEqual(sum(len(c['members']) for c in clusters), 5)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
