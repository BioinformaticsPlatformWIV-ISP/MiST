import json
import random
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import nucmerutils, testingutils
from mist.app.utils.clustersplit import ClusterSplit


class TestClusterSplit(unittest.TestCase):
    """
    Tests the splitting of clusters with sequences that start or end at a different position than the representative.
    Based on NEIS0992, where allele 178 lacks the alternative start codon (TTG) of allele 6, while the representative
    of its cluster (allele 37) contains the start codon and a 3 bp deletion.
    """

    def setUp(self) -> None:
        """
        Creates the test sequences.
        :return: None
        """
        rng = random.Random(1)
        core = 'ATG' + ''.join(rng.choice('ACGT') for _ in range(567)) + 'TAA'
        self.seq_by_id = {
            'locus_37': 'TTG' + core[:300] + core[303:],
            'locus_42': 'TTG' + core[:300] + core[303:450] + ('A' if core[450] != 'A' else 'C') + core[451:],
            'locus_178': core,
        }
        self.dir_temp = testingutils.get_temp_dir()
        self.dir_in = Path(self.dir_temp.name)

    def tearDown(self) -> None:
        """
        Cleans up the temporary directory.
        :return: None
        """
        self.dir_temp.cleanup()

    def _write_fasta(self, seq_ids: list[str], path_out: Path) -> Path:
        """
        Writes the given sequences to a FASTA file.
        :param seq_ids: Sequence ids
        :param path_out: Output path
        :return: Output path
        """
        with path_out.open('w') as handle:
            for seq_id in seq_ids:
                handle.write(f'>{seq_id}\n{self.seq_by_id[seq_id]}\n')
        return path_out

    def test_show_coords(self) -> None:
        """
        Tests the parsing of the show-coords output.
        :return: None
        """
        path_ref = self._write_fasta(['locus_37'], self.dir_in / 'ref.fasta')
        path_query = self._write_fasta(['locus_178'], self.dir_in / 'query.fasta')
        data_coords = nucmerutils.show_coords(nucmerutils.nucmer(path_ref, path_query, self.dir_in))
        row = data_coords.iloc[0]
        self.assertEqual((row['[S1]'], row['[E1]'], row['[S2]'], row['[E2]']), (4, 573, 1, 573))
        self.assertEqual((row['[TAG R]'], row['[TAG Q]']), ('locus_37', 'locus_178'))

    def test_offsets_inverted_repeat(self) -> None:
        """
        Tests that the alignments of an inverted repeat do not determine the offsets. The reference (660 bp) ends with
        the reverse complement of its region 101-200, 'locus_2' has 3 additional bases at the start, and 'locus_3' is
        the reverse complement of 'locus_2'.
        :return: None
        """
        rows = [
            (1, 660, 1, 660, 660, 660, 'locus_1'),
            (101, 200, 660, 561, 100, 660, 'locus_1'),
            (561, 660, 200, 101, 100, 660, 'locus_1'),
            (1, 660, 4, 663, 660, 663, 'locus_2'),
            (101, 200, 663, 564, 100, 663, 'locus_2'),
            (561, 660, 203, 104, 100, 663, 'locus_2'),
            (1, 660, 660, 1, 660, 663, 'locus_3'),
            (101, 200, 1, 100, 100, 663, 'locus_3'),
        ]
        data_coords = pd.DataFrame(rows, columns=['[S1]', '[E1]', '[S2]', '[E2]', '[LEN 2]', '[LEN Q]', '[TAG Q]'])
        data_coords['[LEN R]'] = 660
        self.assertEqual(
            ClusterSplit.calculate_offsets(data_coords),
            {'locus_1': (0, 0), 'locus_2': (3, 0), 'locus_3': (3, 0)},
        )

    def test_split_start_offset(self) -> None:
        """
        Tests that a sequence with a different start position is added as representative.
        :return: None
        """
        members = ['locus_37', 'locus_42', 'locus_178']
        self._write_fasta(members, self.dir_in / 'locus.fasta')
        self._write_fasta(['locus_37'], self.dir_in / 'locus-clustered.fasta')
        with open(self.dir_in / 'mist_db.json', 'w') as handle:
            json.dump(
                {
                    'name': 'locus',
                    'fasta_full': 'locus.fasta',
                    'fasta_clustered': 'locus-clustered.fasta',
                    'clusters': {
                        'locus_37': {
                            'name': 'Cluster_0',
                            'id': 0,
                            'members': [{'seq_id': seq_id, 'ori': '+'} for seq_id in members],
                        }
                    },
                },
                handle,
            )
        ClusterSplit(self.dir_in, debug=False).run()
        with open(self.dir_in / 'locus-clustered.fasta') as handle:
            self.assertEqual([seq.id for seq in SeqIO.parse(handle, 'fasta')], ['locus_37', 'locus_178'])


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
