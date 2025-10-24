import unittest
from importlib.resources import files
from pathlib import Path

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import minimap2utils, testingutils


class TestMinimap2Utils(unittest.TestCase):
    """
    Tests for the Minimap2 utilities.
    """

    def test_minimap2_index(self) -> None:
        """
        Tests the create_index function.
        :return: None
        """
        path_db = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        with testingutils.get_temp_dir() as dir_temp:
            path_in = Path(dir_temp, path_db.name)
            path_in.symlink_to(path_db)
            path_index = minimap2utils.create_index(path_in)
            self.assertGreater(path_index.stat().st_size, 0)

    def test_minimap2_index_and_query(self) -> None:
        """
        Tests the create_index function followed by alignment.
        :return: None
        """
        path_db = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        path_query = Path(str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')))
        with testingutils.get_temp_dir() as dir_temp:
            # Create index
            path_in = Path(dir_temp, path_db.name)
            path_in.symlink_to(path_db)
            path_index = minimap2utils.create_index(path_in)
            self.assertGreater(path_index.stat().st_size, 0)

            # Align to index
            data_out = minimap2utils.align(path_query, path_in)
            self.assertGreater(len(data_out), 0)

    def test_minimap2_query(self) -> None:
        """
        Tests minimap2 query function.
        :return: None
        """
        path_db = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        path_query = Path(str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')))
        data_out = minimap2utils.align(path_query, path_db)
        self.assertGreater(len(data_out), 0)

    def test_minimap2_query_with_cigar(self) -> None:
        """
        Tests minimap2 query function with Cigar string.
        :return: None
        """
        path_db = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        path_query = Path(str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')))
        data_out = minimap2utils.align(path_query, path_db, include_cigar=True)
        self.assertGreater(len(data_out), 0)

if __name__ == '__main__':
    initialize_logging()
    unittest.main()
