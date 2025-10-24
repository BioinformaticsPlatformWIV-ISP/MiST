import unittest
from importlib.resources import files
from pathlib import Path

from mist.app import model
from mist.app.query.allelequeryminimap import AlleleQueryMinimap2, MultiStrategy
from mist.app.utils import testingutils
from mist.scripts.mistindex import MistIndex


class TestAlleleQuery(unittest.TestCase):
    """
    Tests the allele calling method.
    """

    @staticmethod
    def get_expected_tag(allele_result: str | model.Tag) -> model.Tag:
        """
        Returns the expected tag for the given allele result.
        :param allele_result: Allele result
        :return: Expected tag
        """
        if isinstance(allele_result, model.Tag):
            return allele_result
        if allele_result == model.ALLELE_MISSING:
            return model.Tag.ABSENT
        if allele_result.startswith('*'):
            return model.Tag.NOVEL
        if '__' in allele_result:
            return model.Tag.MULTI
        return model.Tag.EXACT

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

    def test_perfect_matches(self) -> None:
        """
        Tests allele calling with perfect matches for both loci.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-perfect_hits.fasta')))
        allele_by_locus = {'NEIS0140-subset': '10', 'NEIS0159-subset': '67'}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            self.assertIn(model.Tag.EXACT, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_novel_allele(self) -> None:
        """
        Tests allele calling with a novel allele.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-novel.fasta')))
        allele_by_locus = {'NEIS0140-subset': '*dd0d', 'NEIS0159-subset': model.ALLELE_MISSING}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_novel_allele_rc(self) -> None:
        """
        Tests allele calling with a novel allele on the reverse strand.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-novel_rc.fasta')))
        allele_by_locus = {'NEIS0140-subset': '*dd0d', 'NEIS0159-subset': model.ALLELE_MISSING}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_multi_hit(self) -> None:
        """
        Tests allele calling with multiple exact hits.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path, multi_strategy=MultiStrategy.ALL)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-multi_hit.fasta')))
        allele_by_locus = {'NEIS0140-subset': '1__108', 'NEIS0159-subset': model.ALLELE_MISSING}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_all_absent(self) -> None:
        """
        Tests allele calling with no hits.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path, multi_strategy=MultiStrategy.ALL)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-both_absent.fasta')))
        allele_by_locus = {'NEIS0140-subset': model.ALLELE_MISSING, 'NEIS0159-subset': model.ALLELE_MISSING}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_one_absent(self) -> None:
        """
        Tests allele calling with no hits.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path, multi_strategy=MultiStrategy.ALL)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-one_absent.fasta')))
        allele_by_locus = {'NEIS0140-subset': model.ALLELE_MISSING, 'NEIS0159-subset': '42'}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")
            self.assertEqual(res.allele_str, allele_by_locus[locus], f"Wrong allele for {locus}")

    def test_edge(self) -> None:
        """
        Tests allele calling with a partial hit on a contig edge.
        :return: None
        """
        caller = AlleleQueryMinimap2(dir_db=self.db_path, multi_strategy=MultiStrategy.ALL)
        path_fasta = Path(str(files('mist').joinpath('resources/testdata/query-edge.fasta')))
        allele_by_locus = {'NEIS0140-subset': model.Tag.EDGE, 'NEIS0159-subset': model.ALLELE_MISSING}
        result_by_locus = caller.query(path_fasta=path_fasta)
        for locus, res in result_by_locus.items():
            expected_tag = TestAlleleQuery.get_expected_tag(allele_by_locus[locus])
            self.assertIn(expected_tag, res.tags, f"No exact match for {locus}")

if __name__ == '__main__':
    unittest.main()
