import json
import unittest
from importlib.resources import files
from pathlib import Path

from click.testing import CliRunner

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import sequenceutils, testingutils
from mist.scripts.cli import cli
from mist.scripts.mistindex import MistIndex


class TestCall(unittest.TestCase):
    """
    Tests the allele calling script.
    """

    def setUp(self) -> None:
        """
        Sets up a temporary directory and builds a database there before each test.
        :return: None
        """
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta')))
        self.dir_temp = testingutils.get_temp_dir()
        self.db_path = Path(self.dir_temp.name)

        # Build the index once for each test
        mist_idx = MistIndex(paths_fasta=[path_fasta], path_profiles=None)
        mist_idx.create_index(dir_out=self.db_path, threads=4)

    def tearDown(self) -> None:
        """
        Clean up the temporary directory after the test.
        :return: None
        """
        self.dir_temp.cleanup()

    def test_call_with_hit(self) -> None:
        """
        Tests querying the database with a hit.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the allele calling
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,
                [
                    'call',
                    '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--threads', '4'
                ], catch_exceptions=False
            )
            self.assertEqual(0, result.exit_code)
            self.assertTrue(path_json.exists())

    def test_call_with_novel_allele(self) -> None:
        """
        Tests querying the database with a novel hit.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the script
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli, [
                    'call',
                    '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query_novel_allele.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--out-dir', str(dir_out),
                    '--threads', '4',
                    '--export-novel'
            ], catch_exceptions=False)

            # Check for the FASTA output
            self.assertTrue((dir_out / 'novel_alleles').exists(), "novel_alleles directory not found")
            path_fasta = next((dir_out / 'novel_alleles').glob('*.fasta'))
            self.assertEqual(sequenceutils.count_sequences(path_fasta), 1)
            self.assertEqual(result.exit_code, 0)

    def test_call_with_novel_allele_rc(self) -> None:
        """
        Tests querying the database with a novel hit.
        Ensures that the same allele id is obtained regardless of strand
        :return: None
        """
        runner = CliRunner()
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
                # noinspection PyTypeChecker
                result = runner.invoke(
                    cli, [
                        'call',
                        '--fasta', str(path_fasta),
                        '--db', str(self.db_path),
                        '--out-json', str(path_json),
                        '--threads', '4',
                    ],
                )
                self.assertEqual(result.exit_code, 0)
                with open(path_json) as handle:
                    calls_out.append(json.load(handle))

        import pprint
        pprint.pprint(calls_out)
        self.assertEqual(
            calls_out[0]['alleles']['NEIS0140-subset']['allele_str'],
            calls_out[1]['alleles']['NEIS0140-subset']['allele_str'],
            "Allele hashes do not match"
        )

    def test_call_no_hit(self) -> None:
        """
        Tests querying the database without a hit.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'

            # Run the script
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,[
                    'call',
                    '--fasta', str(files('mist').joinpath('resources/testdata/neiss_query_no_hit.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--threads', '4',
                ], catch_exceptions=False
            )
            self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
