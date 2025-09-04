import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mist.app.loggers.logger import initialize_logging
from mist.app.utils import sequenceutils, testingutils


class TestSequenceUtils(unittest.TestCase):
    """
    Tests for the sequence utilities.
    """

    def test_hash_sequence(self) -> None:
        """
        Tests the hash sequence function.
        :return: None
        """
        hash_a = sequenceutils.hash_sequence('ACTGA', normalize=False)
        hash_b = sequenceutils.hash_sequence('CGATG', normalize=False)
        self.assertNotEqual(hash_a, hash_b)

    def test_rev_complement(self) -> None:
        """
        Tests the reverse complement function.
        :return: None
        """
        seq_in = 'ACTAGATAGAC'
        self.assertEqual(
            seq_in, sequenceutils.rev_complement(sequenceutils.rev_complement(seq_in))
        )

    def test_merge_fasta(self) -> None:
        """
        Tests the merge fasta function.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            # Create FASTA files
            path_fasta_a = Path(dir_temp, 'seqs_a.fasta')
            with path_fasta_a.open('w') as handle:
                SeqIO.write([SeqRecord(Seq('ACTG'), 'a1'), SeqRecord(Seq('CTGA'), 'a2')], handle, 'fasta')
            path_fasta_b = Path(dir_temp, 'seqs_b.fasta')
            with path_fasta_b.open('w') as handle:
                SeqIO.write([
                    SeqRecord(Seq('CCTG'), 'b1'),
                    SeqRecord(Seq('AAAA'), 'b2'),
                    SeqRecord(Seq('TTTT'), 'b3')],
                handle, 'fasta')

            # Merge FASTA files
            path_merged = Path(dir_temp, 'merged.fasta')
            nb_seqs = sequenceutils.merge_fasta_files([path_fasta_a, path_fasta_b], path_merged)
            self.assertEqual(nb_seqs, 5)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
