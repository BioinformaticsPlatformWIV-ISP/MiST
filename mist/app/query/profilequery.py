from pathlib import Path

import pandas as pd

from mist.app import model
from mist.app.loggers.logger import logger


class ProfileQuery:
    """
    Class to query profiles files with a given allele combination.
    """

    ALLELE_ABSENT = '0'
    ALLELE_WILDCARD = 'N'

    def _parse_profiles(self, path: Path, locus_names: set[str]) -> list[model.Profile]:
        """
        Parses the profile file.
        :param path: Path to the TSV file
        :param locus_names: Locus names
        :return: List of profiles
        """
        # Parse input data
        data_in = pd.read_table(path, dtype=str)
        cols_metadata = [c for c in data_in.columns if c not in locus_names]
        logger.debug(f'Metadata columns: {cols_metadata}')
        cols_alleles = [c for c in data_in.columns if c in locus_names]
        logger.debug(f'Gene columns: {cols_alleles}')

        # Construct the profiles
        profiles = []
        for row in data_in.fillna('n/a').to_dict('records'):
            profiles.append(model.Profile(
                name=row[data_in.columns[0]],
                alleles={c: row[c] for c in cols_alleles},
                metadata=[(c, row[c] if not pd.isna(row[c]) else '-') for c in cols_metadata],
            ))
        logger.debug(f'Parsed {len(profiles):,} profiles')
        return profiles

    @staticmethod
    def _alleles_match(res: model.QueryResult, profile_allele: str) -> bool:
        """
        Checks whether two alleles match.
        :param res: Result
        :param profile_allele: Profile allele
        :return: True if the alleles match, False otherwise
        """
        if profile_allele == ProfileQuery.ALLELE_WILDCARD:
            return True
        if profile_allele == ProfileQuery.ALLELE_ABSENT:
            return res is None

        # No allele detected
        if res is None:
            return False

        # Single allele detected
        if len(res.allele_results) == 1:
            return res.allele_str == profile_allele

        # Multiple alleles detected (-> match if one of them matches)
        return any(allele == profile_allele for allele in res.allele_str.split('__'))

    def __init__(self, path_profiles: Path) -> None:
        """
        Initializes the profiles query class.
        :param path_profiles: Path to profiles file
        :return: None
        """
        locus_names = set(d.name for d in path_profiles.parent.iterdir() if not d.name.startswith('.'))
        self._profiles_by_name = {p.name: p for p in self._parse_profiles(path_profiles, locus_names)}

    def query(self, result_by_locus: dict[str, model.QueryResult | None]) -> tuple[model.Profile, int]:
        """
        Queries the profiles using the detected alleles.
        :param result_by_locus: Detected allele(s) by locus
        :return: Best matching profile, nb. of matching loci
        """
        data_profiles = pd.DataFrame({'profile': [name for name in self._profiles_by_name.keys()]})
        data_profiles['nb_matches'] = data_profiles['profile'].apply(
            lambda p: sum(
                ProfileQuery._alleles_match(res, self._profiles_by_name[p].alleles[locus])
                for locus, res in result_by_locus.items()))
        best_match = data_profiles.sort_values(by='nb_matches', ascending=False).iloc[0]
        return self._profiles_by_name[best_match['profile']], int(best_match['nb_matches'])
