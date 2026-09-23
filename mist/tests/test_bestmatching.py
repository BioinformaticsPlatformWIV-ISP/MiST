import unittest
from pathlib import Path

from mist.app.query.bestmatching import ImperfectMatchDetector, InvalidLengthException
from mist.app.utils import testingutils


class TestImperfectMatchDetector(unittest.TestCase):
    """
    Tests the detection of the best matching imperfect hit.
    """

    ALLELES = {
        'LOCUS_1': 'ACGT' * 25,
        'LOCUS_2': 'ACGT' * 24 + 'ACGA',
        'LOCUS_3': 'ACGTACGTAC',
    }

    def setUp(self) -> None:
        """
        Creates a locus directory with the test alleles.
        :return: None
        """
        self.dir_temp = testingutils.get_temp_dir()
        self.dir_locus = Path(self.dir_temp.name, 'LOCUS')
        self.dir_locus.mkdir()
        with open(self.dir_locus / 'LOCUS.fasta', 'w') as handle:
            for seq_id, seq in TestImperfectMatchDetector.ALLELES.items():
                handle.write(f'>{seq_id}\n{seq}\n')

    def tearDown(self) -> None:
        """
        Clean up the temporary directory after the test.
        :return: None
        """
        self.dir_temp.cleanup()

    def test_best_matching(self) -> None:
        """
        Tests that the closest allele of the same length is returned.
        :return: None
        """
        seq = TestImperfectMatchDetector.ALLELES['LOCUS_1'][:-1] + 'C'
        detector = ImperfectMatchDetector(self.dir_locus)
        self.assertEqual(detector.retrieve_best_matching(seq, min_id=98), ['LOCUS_1', 'LOCUS_2'])

    def test_best_matching_case_insensitive(self) -> None:
        """
        Tests that lowercase (soft-masked) bases do not count as mismatches.
        :return: None
        """
        seq = TestImperfectMatchDetector.ALLELES['LOCUS_1'][:-1] + 'C'
        seq = seq[:50] + seq[50:].lower()
        detector = ImperfectMatchDetector(self.dir_locus)
        self.assertEqual(detector.retrieve_best_matching(seq, min_id=98), ['LOCUS_1', 'LOCUS_2'])

    def test_below_min_id(self) -> None:
        """
        Tests that no alleles are returned when the identity is below the threshold.
        :return: None
        """
        seq = TestImperfectMatchDetector.ALLELES['LOCUS_1'][:-5] + 'CCCCC'
        detector = ImperfectMatchDetector(self.dir_locus)
        self.assertEqual(detector.retrieve_best_matching(seq, min_id=99), [])

    def test_invalid_length(self) -> None:
        """
        Tests that an exception is raised when no alleles have the same length.
        :return: None
        """
        detector = ImperfectMatchDetector(self.dir_locus)
        with self.assertRaises(InvalidLengthException):
            detector.retrieve_best_matching('ACGT', min_id=99)


if __name__ == '__main__':
    unittest.main()
