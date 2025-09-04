import unittest
from importlib.resources import files
from pathlib import Path

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import testingutils, sequenceutils
from mist.scripts.mist_index import main as main_index
from mist.scripts.mist_query import main as main_query


class TestQuery(unittest.TestCase):
    """
    Tests the query functionality.
    """

    def setUp(self) -> None:
        """
        Sets up a temporary directory and builds a database there before each test.
        :return: None
        """
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        self.dir_temp = testingutils.get_temp_dir()
        self.db_path = Path(self.dir_temp.name)

        # Build the database once for the test
        main_index([
            '--fasta', str(path_fasta),
            '--output', str(self.dir_temp.name), '--threads', '4'
        ])

    def tearDown(self) -> None:
        """
        Clean up the temporary directory after the test.
        :return: None
        """
        self.dir_temp.cleanup()

    def test_query_with_hit(self) -> None:
        """
        Tests querying the database with a hit.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the script
            main_query([
                '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query.fasta')),
                '--db', str(self.db_path),
                '--out-json', str(path_json),
                '--threads', '4'
            ])

    def test_query_with_novel_allele(self) -> None:
        """
        Tests querying the database with a novel hit.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the script
            main_query([
                '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query_novel_allele.fasta')),
                '--db', str(self.db_path),
                '--out-json', str(path_json),
                '--out-dir', str(dir_out),
                '--threads', '4',
                '--export-novel'
            ])

            # Check for the FASTA output
            self.assertTrue((dir_out / 'novel_alleles').exists(), "novel_alleles directory not found")
            path_fasta = next((dir_out / 'novel_alleles').glob('*.fasta'))
            self.assertEqual(sequenceutils.count_sequences(path_fasta), 1)

    def test_query_with_novel_allele_rc(self) -> None:
        """
        Tests querying the database with a novel hit.
        Ensures that the same allele id is obtained regardless of strand
        :return: None
        """
        fasta_in = [
            str(files('mist').joinpath('resources/testdata/neiss_query_novel_allele.fasta')),
            str(files('mist').joinpath('resources/testdata/neiss_query_novel_allele_rc.fasta'))
        ]
        calls_out = []
        with testingutils.get_temp_dir() as dir_temp:
            # Run the calling
            for i, path_fasta in enumerate(fasta_in):
                # Output file(s)
                dir_out = Path(dir_temp, f'out_{i}')
                dir_out.mkdir(parents=True, exist_ok=True)
                path_json = dir_out / 'alleles.json'

                # Run the script
                main_query([
                    '--fasta', str(path_fasta),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--threads', '4'
                ])
                calls_out.append(path_json)

    def test_query_no_hit(self) -> None:
        """
        Tests querying the database without a hit.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the script
            main_query([
                '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query_no_hit.fasta')),
                '--db', str(self.db_path),
                '--out-json', str(path_json),
                '--threads', '4'
            ])


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
