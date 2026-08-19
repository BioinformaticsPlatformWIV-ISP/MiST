import unittest
from unittest.mock import Mock

from mist.app import model
from mist.app.utils import alleleutils


class TestAlleleUtils(unittest.TestCase):
    """
    Tests for the allele utilities.
    """

    def test_encode_decode_allele_round_trip(self) -> None:
        """
        Tests that encoding and decoding an allele returns the original value.
        :return: None
        """
        for allele in ['0', '1', '42', model.ALLELE_WILDCARD]:
            self.assertEqual(allele, alleleutils.decode_allele(alleleutils.encode_allele(allele)))

    def test_encode_allele_unreadable(self) -> None:
        """
        Tests that a non-numeric, non-wildcard allele is encoded as unreadable.
        :return: None
        """
        self.assertEqual(alleleutils.CODE_UNREADABLE, alleleutils.encode_allele('n/a'))

    def test_candidate_codes_no_result(self) -> None:
        """
        Tests that a locus without a detected result yields the absent code.
        :return: None
        """
        self.assertEqual({alleleutils.CODE_ABSENT}, alleleutils.candidate_codes(None))

    def test_candidate_codes_single_hit(self) -> None:
        """
        Tests that a single detected allele yields its own code.
        :return: None
        """
        res = Mock(spec=model.QueryResult)
        res.allele_str = '5'
        res.allele_results = [Mock()]
        self.assertEqual({5}, alleleutils.candidate_codes(res))

    def test_candidate_codes_multi_hit_drops_novel_marker(self) -> None:
        """
        Tests that multiple detected alleles yield each code, and a novel-allele marker is dropped.
        :return: None
        """
        res = Mock(spec=model.QueryResult)
        res.allele_str = '5__12*'
        res.allele_results = [Mock(), Mock()]
        self.assertEqual({5}, alleleutils.candidate_codes(res))
