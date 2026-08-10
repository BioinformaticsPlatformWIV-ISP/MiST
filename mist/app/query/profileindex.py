import pickle
from collections import Counter
from pathlib import Path

import pandas as pd

from mist.app import model
from mist.app.loggers.logger import logger
from mist.app.utils import alleleutils

NAME_INDEX = 'profile_index.pkl'


class ProfileIndex:
    """
    Class to query profiles files with a given allele combination, using an inverted index for speed.

    Groups profile indices by (locus, allele) at construction time, so a query only has to look up the profiles that
    share a detected allele at each locus, instead of comparing every profile against every queried locus. Since
    building the index means reading through the whole profiles file, it can be persisted with `save` and reloaded
    with `dir_index` instead of rebuilt on every query.

    Alleles are stored internally as compact int codes rather than strings.
    """

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
            self._loci: list[str] = state['loci']
            self._names: list[str] = state['names']
            self._metadata: list[list[tuple[str, str]]] = state['metadata']
            self._allele_codes: list[list[int]] = state['allele_codes']
            self._buckets: dict[str, dict[int, list[int]]] = state['buckets']
            logger.debug(f'Loaded index: {len(self._names):,} profiles ({dir_index})')
        else:
            self._loci, self._names, self._metadata, self._allele_codes = self._parse_profiles(
                path_profiles, set(loci)
            )
            self._buckets = self._build_buckets(self._loci, self._allele_codes)

    def save(self, dir_out: Path) -> None:
        """
        Persists the index to disk, so it can be reloaded with `dir_index` instead of rebuilt on every query.
        :param dir_out: Output directory
        :return: None
        """
        dir_out.mkdir(parents=True, exist_ok=True)
        with open(dir_out / NAME_INDEX, 'wb') as handle:
            pickle.dump(
                {
                    'loci': self._loci,
                    'names': self._names,
                    'metadata': self._metadata,
                    'allele_codes': self._allele_codes,
                    'buckets': self._buckets,
                },
                handle,
            )
        logger.info(f'Profiles index saved: {dir_out / NAME_INDEX}')

    def _parse_profiles(
        self, path: Path, locus_names: set[str]
    ) -> tuple[list[str], list[str], list[list[tuple[str, str]]], list[list[int]]]:
        """
        Parses the profile file.
        :param path: Path to the TSV file
        :param locus_names: Locus names
        :return: Locus names (in file column order), profile names, metadata, and encoded allele codes (one list
            per profile, aligned to the returned locus names)
        """
        data_in = pd.read_table(path, dtype=str).fillna('n/a')
        cols_metadata = [c for c in data_in.columns if c not in locus_names]
        logger.debug(f'Metadata columns: {cols_metadata}')
        cols_alleles = [c for c in data_in.columns if c in locus_names]
        logger.debug(f'Gene columns: {cols_alleles}')

        names = data_in[data_in.columns[0]].tolist()
        metadata = [list(zip(cols_metadata, row)) for row in data_in[cols_metadata].to_numpy()]
        allele_codes = [[alleleutils.encode_allele(val) for val in row] for row in data_in[cols_alleles].to_numpy()]
        logger.debug(f'Parsed {len(names):,} profiles')
        return cols_alleles, names, metadata, allele_codes

    @staticmethod
    def _build_buckets(loci: list[str], allele_codes: list[list[int]]) -> dict[str, dict[int, list[int]]]:
        """
        Groups profile indices by (locus, allele code), so a query only needs to look up the profiles that share a
        detected allele at a given locus, instead of scanning every profile.
        :param loci: Locus names, aligned to the columns of `allele_codes`
        :param allele_codes: Encoded allele codes per profile (one list per profile, aligned to `loci`)
        :return: Mapping of locus -> allele code -> profile indices with that code
        """
        buckets: dict[str, dict[int, list[int]]] = {locus: {} for locus in loci}
        for idx, codes in enumerate(allele_codes):
            for locus, code in zip(loci, codes):
                buckets[locus].setdefault(code, []).append(idx)
        logger.debug(f'Built index: {len(allele_codes):,} profiles, {len(loci):,} loci')
        return buckets

    def _build_profile(self, idx: int) -> model.Profile:
        """
        Builds a `Profile` object for the given row index, decoding its allele codes back into strings.
        :param idx: Row index
        :return: Profile
        """
        return model.Profile(
            name=self._names[idx],
            alleles={
                locus: alleleutils.decode_allele(code) for locus, code in zip(self._loci, self._allele_codes[idx])
            },
            metadata=self._metadata[idx],
        )

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
            for code in alleleutils.candidate_codes(res) | {alleleutils.CODE_WILDCARD}:
                counts.update(bucket.get(code, []))

        if not counts:
            # No profile matched anything (incl. wildcards) at any queried locus - nothing useful to report.
            return [], -1 if len(self._names) == 0 else 0

        best_matches = max(counts.values())
        best_indices = sorted(idx for idx, nb in counts.items() if nb == best_matches)
        return [self._build_profile(idx) for idx in best_indices], best_matches
