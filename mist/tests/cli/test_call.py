import json
import unittest
from importlib.resources import files
from pathlib import Path

from click.testing import CliRunner

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import sequenceutils, testingutils
from mist.scripts.cli import cli
from mist.scripts.mistdists import MistDists
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
        paths_fasta = [
            Path(str(files('mist').joinpath('resources/testdata/NEIS0140-subset.fasta'))),
            Path(str(files('mist').joinpath('resources/testdata/NEIS0159-subset.fasta'))),
        ]
        self.dir_temp = testingutils.get_temp_dir()
        self.db_path = Path(self.dir_temp.name)

        # Build the index once for each test
        mist_idx = MistIndex(paths_fasta=paths_fasta, path_profiles=None)
        mist_idx.create_index(dir_out=self.db_path, threads=4)

    def tearDown(self) -> None:
        """
        Clean up the temporary directory after the test.
        :return: None
        """
        self.dir_temp.cleanup()

    def test_call_with_hits(self) -> None:
        """
        Tests calling the alleles with all perfect hits.
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
                    '--fasta', str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--threads', '4'
                ], catch_exceptions=False
            )

            # Check if the command was executed successfully
            self.assertEqual(0, result.exit_code)
            self.assertTrue(path_json.exists(), "Output JSON file not generated")

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
                    '--fasta', str(files('mist').joinpath('resources/testdata/query-novel.fasta')),
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
            str(files('mist').joinpath('resources/testdata/query-novel.fasta')),
            str(files('mist').joinpath('resources/testdata/query-novel_rc.fasta'))
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
                    '--fasta', str(files('mist').joinpath('resources/testdata/query-both_absent.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--threads', '4',
                ], catch_exceptions=False
            )
            self.assertEqual(result.exit_code, 0)

    def test_call_sample_id(self) -> None:
        """
        Tests querying the database without a sample id.
        :return: None
        """
        runner = CliRunner()
        with testingutils.get_temp_dir() as dir_temp:
            # Output file(s)
            dir_out = Path(dir_temp, 'out')
            dir_out.mkdir(parents=True, exist_ok=True)
            path_json = dir_out / 'alleles.json'
            path_tsv = dir_out / 'alleles.tsv'

            # Run the script
            # noinspection PyTypeChecker
            result = runner.invoke(
                cli,[
                    'call',
                    '--fasta', str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')),
                    '--db', str(self.db_path),
                    '--out-json', str(path_json),
                    '--out-tsv', str(path_tsv),
                    '--threads', '4',
                    '--sample-id', 'my_sample_id',
                ], catch_exceptions=False
            )
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(MistDists.parse_tsv(path_tsv)[0], 'my_sample_id')
            self.assertEqual(MistDists.parse_json(path_json)[0], 'my_sample_id')


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
