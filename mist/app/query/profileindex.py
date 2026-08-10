import pickle
from collections import Counter
from pathlib import Path

import pandas as pd

from mist.app import model
from mist.app.loggers.logger import logger

NAME_INDEX = 'profile_index.pkl'


class ProfileIndex:
    """
    Class to query profiles files with a given allele combination, using an inverted index for speed.

    Groups profile indices by (locus, allele) at construction time, so a query only has to look up the profiles that
    share a detected allele at each locus, instead of comparing every profile against every queried locus. Since
    building the index means reading through the whole profiles file, it can be persisted with `save` and reloaded
    with `dir_index` instead of rebuilt on every query.
    """

    ALLELE_ABSENT = '0'
    ALLELE_WILDCARD = 'N'

    def __init__(
        self, path_profiles: Path | None = None, loci: list[str] | None = None, dir_index: Path | None = None
    ) -> None:
        """
        Initializes the profile index, either by parsing a profiles file and building the index from scratch,
        or by reloading a previously saved one.
        :param path_profiles: Path to profiles file (builds the index from scratch, together with `loci`)
        :param loci: List of loci (builds the index from scratch, together with `path_profiles`)
        :param dir_index: Directory containing a previously saved index (see `save`), reloaded instead of rebuilt
        :return: None
        """
        if dir_index is not None:
            with open(dir_index / NAME_INDEX, 'rb') as handle:
                state = pickle.load(handle)
            self._profiles: list[model.Profile] = state['profiles']
            self._buckets: dict[str, dict[str, list[int]]] = state['buckets']
            logger.debug(f'Loaded index: {len(self._profiles):,} profiles ({dir_index})')
        else:
            self._profiles = self._parse_profiles(path_profiles, set(loci))
            self._buckets = self._build_buckets(self._profiles, set(loci))

    def save(self, dir_out: Path) -> None:
        """
        Persists the index to disk, so it can be reloaded with `dir_index` instead of rebuilt on every query.
        :param dir_out: Output directory
        :return: None
        """
        dir_out.mkdir(parents=True, exist_ok=True)
        with open(dir_out / NAME_INDEX, 'wb') as handle:
            pickle.dump({'profiles': self._profiles, 'buckets': self._buckets}, handle)
        logger.info(f'Profiles index saved: {dir_out / NAME_INDEX}')

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
            profiles.append(
                model.Profile(
                    name=row[data_in.columns[0]],
                    alleles={c: row[c] for c in cols_alleles},
                    metadata=[(c, row[c] if not pd.isna(row[c]) else '-') for c in cols_metadata],
                )
            )
        logger.debug(f'Parsed {len(profiles):,} profiles')
        return profiles

    @staticmethod
    def _build_buckets(profiles: list[model.Profile], loci: set[str]) -> dict[str, dict[str, list[int]]]:
        """
        Groups profile indices by (locus, allele), so a query only needs to look up the profiles that share a detected
        allele at a given locus, instead of scanning every profile.
        :param profiles: Profiles
        :param loci: Locus names
        :return: Mapping of locus -> allele -> profile indices with that allele
        """
        buckets: dict[str, dict[str, list[int]]] = {locus: {} for locus in loci}
        for idx, profile in enumerate(profiles):
            for locus in loci:
                bucket = buckets[locus].setdefault(profile.alleles[locus], [])
                bucket.append(idx)
        logger.debug(f'Built index: {len(profiles):,} profiles, {len(loci):,} loci')
        return buckets

    @staticmethod
    def _candidate_alleles(res: model.QueryResult | None) -> set[str]:
        """
        Determines which stored profile alleles a detected result could match at a single locus.
        :param res: Detected result for the locus (None if nothing was detected)
        :return: Candidate allele strings
        """
        if res is None:
            return {ProfileIndex.ALLELE_ABSENT}
        if len(res.allele_results) == 1:
            return {res.allele_str}
        return set(res.allele_str.split('__'))

    def query(self, result_by_locus: dict[str, model.QueryResult | None]) -> tuple[list[model.Profile], int]:
        """
        Queries the index using the detected alleles.
        :param result_by_locus: Detected allele(s) by locus
        :return: Best matching profiles, nb. of matching loci
        """
        counts: Counter[int] = Counter()
        for locus, res in result_by_locus.items():
            bucket = self._buckets[locus]
            # Profiles with a wildcard at this locus always match, regardless of the query
            for allele in self._candidate_alleles(res) | {ProfileIndex.ALLELE_WILDCARD}:
                counts.update(bucket.get(allele, []))

        if not counts:
            # No profile matched anything (incl. wildcards) at any queried locus - nothing useful to report.
            return [], -1 if len(self._profiles) == 0 else 0

        best_matches = max(counts.values())
        best_indices = sorted(idx for idx, nb in counts.items() if nb == best_matches)
        return [self._profiles[idx] for idx in best_indices], best_matches
