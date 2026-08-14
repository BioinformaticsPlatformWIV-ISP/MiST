from pathlib import Path
from unittest.mock import Mock

import pandas as pd

from mist.app import model
from mist.app.query.profileindex import ProfileIndex
from mist.app.utils import testingutils
import unittest


class TestProfileIndex(unittest.TestCase):
    """
    Contains tests for performing profile queries using the index.
    """

    TEST_PROFILES = [
        {'ST': 'ST1', 'locus1': '1', 'locus2': '1', 'locus3': '1'},
        {'ST': 'ST2', 'locus1': '2', 'locus2': 'N', 'locus3': '2'},
        {'ST': 'ST2_alt', 'locus1': '2', 'locus2': '2', 'locus3': 'N'},
        {'ST': 'ST3', 'locus1': '3', 'locus2': 'N', 'locus3': 'N'},
        {'ST': 'ST4', 'locus1': '4', 'locus2': 'N', 'locus3': 'N'},
    ]

    def test_query_returns_best_matching_profiles(self) -> None:
        """
        Tests the query method returns the best matching profiles.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)

            # Create the profiles file
            path_profiles = dir_temp / "profiles.tsv"
            pd.DataFrame(TestProfileIndex.TEST_PROFILES).to_csv(path_profiles, sep='\t', index=False)

            # Create query results
            result_locus1 = Mock(spec=model.QueryResult)
            result_locus1.allele_str = "1"
            result_locus1.allele_results = [Mock()]

            result_locus2 = Mock(spec=model.QueryResult)
            result_locus2.allele_str = "1"
            result_locus2.allele_results = [Mock()]

            result_locus3 = Mock(spec=model.QueryResult)
            result_locus3.allele_str = "1"
            result_locus3.allele_results = [Mock()]

            result_by_locus = {
                "locus1": result_locus1,
                "locus2": result_locus2,
                "locus3": result_locus3,
            }

            # Query the profiles
            pi = ProfileIndex(
                path_profiles=path_profiles, loci=list(result_by_locus.keys()), dir_out=dir_temp / "index"
            )
            profiles, nb_matches = pi.query(result_by_locus)

            # Assertions
            self.assertEqual(len(profiles), 1)
            self.assertEqual(profiles[0].name, "ST1")
            self.assertEqual(nb_matches, 3)

    def test_query_returns_multiple_best_matching_profiles(self) -> None:
        """
        Tests the query method returns the best matching profiles (multiple).
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)

            # Create the profiles file
            path_profiles = dir_temp / "profiles.tsv"
            pd.DataFrame(TestProfileIndex.TEST_PROFILES).to_csv(path_profiles, sep='\t', index=False)

            # Create query results
            result_locus1 = Mock(spec=model.QueryResult)
            result_locus1.allele_str = "2"
            result_locus1.allele_results = [Mock()]

            result_locus2 = Mock(spec=model.QueryResult)
            result_locus2.allele_str = "2"
            result_locus2.allele_results = [Mock()]

            result_locus3 = Mock(spec=model.QueryResult)
            result_locus3.allele_str = "2"
            result_locus3.allele_results = [Mock()]

            result_by_locus = {
                "locus1": result_locus1,
                "locus2": result_locus2,
                "locus3": result_locus3,
            }

            # Query the profiles
            pi = ProfileIndex(
                path_profiles=path_profiles, loci=list(result_by_locus.keys()), dir_out=dir_temp / "index"
            )
            profiles, nb_matches = pi.query(result_by_locus)

            # Assertions
            self.assertEqual(len(profiles), 2)
            self.assertIn("ST2", [p.name for p in profiles])
            self.assertIn("ST2_alt", [p.name for p in profiles])
            self.assertEqual(nb_matches, 3)

    def test_build_and_reload_index_returns_same_result(self) -> None:
        """
        Tests that a built index, once reloaded, returns the same query result as the original.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)

            # Create the profiles file
            path_profiles = dir_temp / "profiles.tsv"
            pd.DataFrame(TestProfileIndex.TEST_PROFILES).to_csv(path_profiles, sep='\t', index=False)

            # Create query results
            result_locus1 = Mock(spec=model.QueryResult)
            result_locus1.allele_str = "1"
            result_locus1.allele_results = [Mock()]

            result_by_locus = {"locus1": result_locus1}

            # Build and reload the index
            dir_index = dir_temp / "index"
            ProfileIndex(path_profiles=path_profiles, loci=["locus1", "locus2", "locus3"], dir_out=dir_index)
            pi_reloaded = ProfileIndex(dir_index=dir_index)

            # Assertions
            profiles, nb_matches = pi_reloaded.query(result_by_locus)
            self.assertEqual(len(profiles), 1)
            self.assertEqual(profiles[0].name, "ST1")
            self.assertEqual(nb_matches, 1)
