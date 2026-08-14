import json
import re
from pathlib import Path
from typing import Any

from furl import furl

from mist.app import NAME_DB_INFO, errors
from mist.app.loggers.logger import logger
from mist.app.utils import restutils

URL_ENTEROBASE_API = 'https://enterobase.warwick.ac.uk/api/v2.0'


class MistLinCode:
    """
    Extracts LIN codes for the best-matching profile in a `mist call` JSON output.

    LIN codes come from profile metadata for BIGSdb databases or the EnteroBase API for EnteroBase databases.
    Cached bin thresholds determine how many LIN-code positions can be reliably assigned, returning both the full and
    partially masked codes.
    """

    def __init__(
        self,
        dir_db: Path,
        entero_token: str | None = None,
        entero_species: str | None = None,
        entero_scheme: str | None = None,
    ) -> None:
        """
        Initializes the extractor.
        :param dir_db: Database directory (indexed by `mist index`)
        :param entero_token: EnteroBase API token
        :param entero_species: EnteroBase API species name
        :param entero_scheme: EnteroBase API scheme name
        :return: None
        """
        self._entero_token = entero_token
        self._entero_species = entero_species
        self._entero_scheme = entero_scheme
        path_db_info = dir_db / NAME_DB_INFO
        if not path_db_info.exists():
            raise errors.LinCodeError(
                f"'{NAME_DB_INFO}' not found in {dir_db} - LIN-code extraction needs to know which source "
                f"the database came from, so it must have been downloaded with 'mist download'"
            )
        with path_db_info.open() as handle:
            self._db_info: dict[str, Any] = json.load(handle)

    def run(self, path_mist_json: Path, path_out: Path) -> None:
        """
        Extracts the LIN code for the best-matching profile in a `mist call` JSON output file, and writes the
        result to a JSON file.
        :param path_mist_json: Path to the `mist call` JSON output
        :param path_out: LIN-code output JSON file
        :return: None
        """
        with path_mist_json.open() as handle:
            data_mist = json.load(handle)
        result = self.extract(data_mist)

        with path_out.open('w') as handle:
            json.dump(result, handle, indent=2)

        if result['lincode_partial'] is not None:
            rendered = '-'.join(v if v is not None else '*' for v in result['lincode_partial'])
            logger.info(f"Extracted LIN code for ST {result['st']}: {rendered}")
        else:
            logger.info(f"No LIN code for ST {result['st']}")

    def extract(self, data_mist: dict[str, Any]) -> dict[str, Any]:
        """
        Extracts the LIN code for the best-matching profile in a `mist call` result.
        :param data_mist: Parsed `mist call` JSON output
        :return: LIN-code result
        """
        profiles = data_mist.get('profiles')
        if not profiles:
            raise errors.LinCodeError('No matching profile found in the input - LIN code cannot be determined')
        if len(profiles) > 1:
            logger.warning('Multiple equivalent matching profiles found - using the first one')
        profile = profiles[0]

        downloader = self._db_info.get('downloader')
        if downloader in ('bigsdb', 'bigsdb_auth'):
            return self._extract_bigsdb(profile)
        if downloader == 'enterobase':
            return self._extract_enterobase(profile)
        raise errors.LinCodeError(f"LIN-code extraction is not supported for downloader '{downloader}'")

    def _extract_bigsdb(self, profile: dict[str, Any]) -> dict[str, Any]:
        """
        Extracts the LIN code for a BIGSdb-sourced profile.
        :param profile: Best-matching profile (from the `mist call` JSON output)
        :return: LIN-code result
        """
        metadata = dict(profile['metadata'])
        lincode_full = _split_lincode(_find_metadata_value(metadata, ['lincode', 'lin_code', 'lin code']))

        lincodes = self._db_info.get('lincodes')
        if not lincodes:
            raise errors.LinCodeError(
                f"'{NAME_DB_INFO}' has no cached LIN-code thresholds. The scheme at {self._db_info['url']} "
                "may not define them, or this database was downloaded with an older MiST version. "
                f"Re-run 'mist download' to refresh {NAME_DB_INFO}."
            )
        thresholds = [int(t) for t in lincodes['thresholds'].split(';')]

        nb_matches, nb_loci = profile['nb_matches'], len(profile['alleles'])
        nb_diff = nb_loci - nb_matches
        bin_assigned = _determine_bin(nb_diff, thresholds)
        lincode_partial = _mask_lincode(lincode_full, bin_assigned)

        fields = {}
        for field in lincodes.get('fields', []):
            name = field['field']
            fields[name] = _find_metadata_value(metadata, [name], required=False)

        return {
            'st': profile['name'],
            'nb_matches': nb_matches,
            'nb_loci': nb_loci,
            'lincode_full': lincode_full,
            'lincode_partial': lincode_partial,
            'fields': fields,
        }

    def _extract_enterobase(self, profile: dict[str, Any]) -> dict[str, Any]:
        """
        Extracts the LIN code for an EnteroBase-sourced profile via the EnteroBase API. The bin thresholds needed to
        mask the LIN code come from the ST's hierCC distances.
        :param profile: Best-matching profile
        :return: LIN-code result
        """
        if not self._entero_token:
            raise errors.LinCodeError('An EnteroBase API token is required to extract LIN codes for this database')

        species_guess, scheme_guess = _entero_species_scheme(self._db_info['url'])
        species = self._entero_species or species_guess
        scheme = self._entero_scheme or scheme_guess
        st_id = profile['name']
        url = f'{URL_ENTEROBASE_API}/{species}/{scheme}/sts?limit=1&offset=0&st_id={st_id}&scheme={scheme}'
        sts = restutils.retrieve_page_data(url, auth=(self._entero_token, '')).json().get('STs', [])
        if not sts:
            raise errors.LinCodeError(
                f"No LIN code found for ST {st_id} at {url} - the species/scheme slug ('{species}'/'{scheme}') "
                f"may be wrong for this database; override with --entero-species/--entero-scheme if so "
                f"(EnteroBase's API naming doesn't reliably match its download URLs)"
            )
        info = sts[0]['info']
        lincode_full = _split_lincode(info['lin_code'])
        nb_matches, nb_loci = profile['nb_matches'], len(profile['alleles'])

        # Prefer 'hierCCv0' over 'hierCC'
        hiercc_raw = info.get('hierCCv0') or info.get('hierCC', {})
        thresholds = _entero_hiercc_thresholds(hiercc_raw)
        if len(thresholds) == len(lincode_full):
            nb_assigned = _determine_bin(nb_loci - nb_matches, thresholds)
            lincode_partial = _mask_lincode(lincode_full, nb_assigned)
        else:
            logger.warning(
                f'hierCC thresholds ({len(thresholds)}) do not match the LIN code length '
                f'({len(lincode_full)}) for ST {st_id} - returning the LIN code unmasked'
            )
            lincode_partial = None

        return {
            'st': st_id,
            'nb_matches': nb_matches,
            'nb_loci': nb_loci,
            'lincode_full': lincode_full,
            'lincode_partial': lincode_partial,
        }


def _find_metadata_value(metadata: dict[str, str], candidate_keys: list[str], required: bool = True) -> str | None:
    """
    Looks up a metadata value by a case-insensitive match against a list of candidate keys.
    :param metadata: Profile metadata
    :param candidate_keys: Candidate keys to try, case-insensitively
    :param required: If True, raises when no candidate key is found
    :return: Metadata value, or None if not found and not required
    """
    by_lower = {k.lower(): v for k, v in metadata.items()}
    for candidate in candidate_keys:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    if required:
        raise errors.LinCodeError(f'None of the expected metadata field(s) {candidate_keys} found in profile metadata')
    return None


def _split_lincode(lincode: str) -> list[str]:
    """
    Splits a LIN code string into its positions, regardless of whether the source uses '_' or '-' separators.
    :param lincode: LIN code string
    :return: LIN code positions
    """
    return re.split(r'[_-]', lincode)


def _determine_bin(nb_diff: int, thresholds: list[int]) -> int:
    """
    Determines the number of leading LIN-code positions that can be reliably assigned, given a descending
    list of per-position mismatch thresholds (thresholds[i] is the max. nb. of mismatches for which position
    i remains reliably assignable) and the number of mismatches to the best-matching profile.
    :param nb_diff: Number of mismatching loci
    :param thresholds: Per-position mismatch-count thresholds, descending
    :return: Nb. of reliably assignable leading positions
    """
    for i, threshold in enumerate(thresholds):
        if threshold < nb_diff:
            return i
    return len(thresholds)


def _mask_lincode(lincode_full: list[str], nb_assigned: int) -> list[str | None]:
    """
    Masks LIN-code positions beyond the reliably assignable ones with None.
    :param lincode_full: Full LIN code, as a list of position values
    :param nb_assigned: Nb. of reliably assignable leading positions (from `_determine_bin`)
    :return: Partial LIN code, same length as `lincode_full`
    """
    return [(value if i < nb_assigned else None) for i, value in enumerate(lincode_full)]


def _entero_hiercc_thresholds(hiercc: dict[str, Any]) -> list[int]:
    """
    Parses EnteroBase's hierCC distances (the 'dNNN' keys of a ST's 'hierCC'/'hierCCv0' field, e.g. 'd0',
    'd2', 'd2850') into a descending list of mismatch-count thresholds, equivalent to BIGSdb's LIN-code
    thresholds. 'd0' (identical profile, i.e. 100% LIN-code similarity) is kept - EnteroBase's documented
    per-species threshold tables (see the 'LIN codes' wiki page) list it as the deepest real level, not a
    sentinel to discard.
    :param hiercc: The 'hierCC' (or 'hierCCv0') object from an EnteroBase ST API response
    :return: Mismatch-count thresholds, descending
    """
    distances = []
    for key in hiercc:
        match = re.fullmatch(r'd(\d+)', key)
        if match:
            distances.append(int(match.group(1)))
    return sorted(distances, reverse=True)


def _entero_species_scheme(url: str) -> tuple[str, str]:
    """
    Attempts to guess the EnteroBase species and scheme name from the URL (unreliable!).
    :param url: Scheme URL
    :return: (species, scheme) API slugs
    """
    segment = next(s for s in reversed(furl(url).path.segments) if s)
    species_part, scheme_part = segment.split('.', 1)
    return species_part.lower(), re.sub(r'v(\d+)$', r'_v\1', scheme_part)
