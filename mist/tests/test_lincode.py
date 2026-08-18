import json
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from mist.app.errors import LinCodeError
from mist.app.utils import testingutils
from mist.scripts.mistlincode import (
    MistLinCode,
    _determine_bin,
    _entero_hiercc_thresholds,
    _entero_species_scheme,
    _mask_lincode,
)


class TestLinCodeHelpers(unittest.TestCase):
    """
    Tests for the pure LIN-code helper functions.
    """

    # Real thresholds/LIN code for the Klebsiella scgMLST629_S scheme (10 positions).
    THRESHOLDS = [629, 610, 585, 190, 43, 10, 7, 4, 2, 1]
    LINCODE_FULL = ['0', '0', '369', '0', '0', '0', '0', '35', '0', '0']

    def test_determine_bin(self) -> None:
        """
        Tests that the number of assignable positions matches the known Klebsiella example (5 mismatches):
        position 6 (threshold 7) still tolerates 5 mismatches, position 7 (threshold 4) does not, so 7
        leading positions (indices 0-6) are assignable.
        :return: None
        """
        self.assertEqual(7, _determine_bin(5, self.THRESHOLDS))

    def test_determine_bin_exact_match(self) -> None:
        """
        Tests that zero mismatches assigns every position, including the deepest one.
        :return: None
        """
        self.assertEqual(len(self.THRESHOLDS), _determine_bin(0, self.THRESHOLDS))

    def test_mask_lincode_matches_known_example(self) -> None:
        """
        Tests that masking with the known example's assignable-position count reproduces the documented
        partial LIN code.
        :return: None
        """
        nb_assigned = _determine_bin(5, self.THRESHOLDS)
        masked = _mask_lincode(self.LINCODE_FULL, nb_assigned)
        self.assertEqual(['0', '0', '369', '0', '0', '0', '0', None, None, None], masked)

    def test_entero_species_scheme(self) -> None:
        """
        Tests deriving the EnteroBase API species/scheme slugs from a scheme download URL.
        :return: None
        """
        species, scheme = _entero_species_scheme('https://enterobase.warwick.ac.uk/schemes/Senterica.cgMLSTv2/')
        self.assertEqual('senterica', species)
        self.assertEqual('cgMLST_v2', scheme)

    def test_entero_hiercc_thresholds(self) -> None:
        """
        Tests parsing hierCC distances (real EnteroBase example) into descending thresholds, keeping d0 as
        the deepest (100%-similarity) level.
        :return: None
        """
        hiercc = {
            'd0': '102245', 'd2': '102245', 'd5': '102245', 'd10': '102245', 'd20': '78391', 'd50': '4323',
            'd100': '36', 'd150': '36', 'd200': '36', 'd400': '36', 'd900': '36', 'd2000': '36', 'd2600': '2',
            'd2850': '2',
        }
        self.assertEqual(
            [2850, 2600, 2000, 900, 400, 200, 150, 100, 50, 20, 10, 5, 2, 0], _entero_hiercc_thresholds(hiercc)
        )


class TestMistLinCode(unittest.TestCase):
    """
    Tests for the `MistLinCode` class.
    """

    def _write_db_info(
        self, dir_db: Path, downloader: str, url: str, lincodes: dict[str, object] | None = None
    ) -> None:
        """
        Writes a minimal `db_info.json` file, as produced by `mist download`.
        :param dir_db: Database directory
        :param downloader: Downloader key
        :param url: Scheme URL
        :param lincodes: Cached LIN-code thresholds/fields, as fetched at download time for BIGSdb sources
        :return: None
        """
        data = {'url': url, 'downloader': downloader, 'download_date': '2026-08-10T00:00:00'}
        if lincodes is not None:
            data['lincodes'] = lincodes
        with (dir_db / 'db_info.json').open('w') as handle:
            json.dump(data, handle)

    def test_missing_db_info_raises(self) -> None:
        """
        Tests that a database without db_info.json raises a clear error.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            with self.assertRaises(LinCodeError):
                MistLinCode(dir_db=Path(dir_temp))

    def test_no_matching_profile_raises(self) -> None:
        """
        Tests that an input with no matching profile raises a clear error.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(dir_temp, downloader='bigsdb', url='https://example.org/scheme')
            extractor = MistLinCode(dir_db=dir_temp)
            with self.assertRaises(LinCodeError):
                extractor.extract({'profiles': None})

    def test_unsupported_downloader_raises(self) -> None:
        """
        Tests that a downloader without LIN-code support raises a clear error.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(dir_temp, downloader='cgmlstorg', url='https://example.org/scheme')
            extractor = MistLinCode(dir_db=dir_temp)
            with self.assertRaises(LinCodeError):
                extractor.extract({'profiles': [{'name': 'ST1', 'metadata': [], 'alleles': {}}]})

    def test_extract_bigsdb(self) -> None:
        """
        Tests extracting a partial LIN code for a BIGSdb-sourced profile end-to-end, using LIN-code thresholds
        cached in `db_info.json` at download time (no live API call needed for BIGSdb extraction).
        :return: None
        """
        lincodes = {
            'thresholds': '629;610;585;190;43;10;7;4;2;1',
            'fields': [{'field': 'Phylogroup'}],
        }
        profile = {
            'name': '4362',
            'nb_matches': 624,
            'alleles': {f'locus{i}': '1' for i in range(629)},
            'metadata': [['LINcode', '0_0_369_0_0_0_0_35_0_0'], ['Phylogroup', 'KpI']],
        }

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            url = 'https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18'
            self._write_db_info(dir_temp, downloader='bigsdb', url=url, lincodes=lincodes)
            result = MistLinCode(dir_db=dir_temp).extract({'profiles': [profile]})

        self.assertEqual('4362', result['st'])
        self.assertEqual(['0', '0', '369', '0', '0', '0', '0', '35', '0', '0'], result['lincode_full'])
        self.assertEqual(['0', '0', '369', '0', '0', '0', '0', None, None, None], result['lincode_partial'])
        self.assertEqual([629, 610, 585, 190, 43, 10, 7, 4, 2, 1], result['thresholds'])
        self.assertEqual({'Phylogroup': 'KpI'}, result['fields'])

    def test_run_writes_output_file(self) -> None:
        """
        Tests that `run` reads the `mist call` JSON input, extracts the LIN code, and writes the result to
        the output JSON file (the same result `extract` would return directly).
        :return: None
        """
        lincodes = {
            'thresholds': '629;610;585;190;43;10;7;4;2;1',
            'fields': [{'field': 'Phylogroup'}],
        }
        profile = {
            'name': '4362',
            'nb_matches': 624,
            'alleles': {f'locus{i}': '1' for i in range(629)},
            'metadata': [['LINcode', '0_0_369_0_0_0_0_35_0_0'], ['Phylogroup', 'KpI']],
        }

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            url = 'https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18'
            self._write_db_info(dir_temp, downloader='bigsdb', url=url, lincodes=lincodes)

            path_mist_json = dir_temp / 'mist.json'
            with path_mist_json.open('w') as handle:
                json.dump({'profiles': [profile]}, handle)
            path_out = dir_temp / 'lincode.json'

            MistLinCode(dir_db=dir_temp).run(path_mist_json, path_out)

            with path_out.open() as handle:
                result = json.load(handle)

        self.assertEqual('4362', result['st'])
        self.assertEqual(['0', '0', '369', '0', '0', '0', '0', '35', '0', '0'], result['lincode_full'])
        self.assertEqual(['0', '0', '369', '0', '0', '0', '0', None, None, None], result['lincode_partial'])
        self.assertEqual([629, 610, 585, 190, 43, 10, 7, 4, 2, 1], result['thresholds'])

    def test_extract_bigsdb_missing_cached_thresholds_raises(self) -> None:
        """
        Tests that a BIGSdb-sourced database without cached LIN-code thresholds (e.g. downloaded before this
        was cached, or a scheme that doesn't define them) raises a clear, actionable error instead of trying
        a live API call.
        :return: None
        """
        profile = {
            'name': '4362',
            'nb_matches': 624,
            'alleles': {f'locus{i}': '1' for i in range(629)},
            'metadata': [['LINcode', '0_0_369_0_0_0_0_35_0_0']],
        }
        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(dir_temp, downloader='bigsdb', url='https://example.org/scheme')
            with self.assertRaises(LinCodeError):
                MistLinCode(dir_db=dir_temp).extract({'profiles': [profile]})

    # Real EnteroBase example response (ST 102245, senterica/cgMLST_v2). 'hierCC' is the current pass, which
    # has since picked up an extra 'd150' level (14 raw keys) not part of EnteroBase's documented 13-level
    # LIN-code mapping; 'hierCCv0' is the original pass LIN codes were mapped from and matches that mapping
    # exactly (13 raw keys, no 'd150').
    ENTERO_INFO = {
        'lin_code': '0-0-3-0-0-0-0-3-3-7-0-0-0',
        'hierCC': {
            'd0': '102245', 'd2': '102245', 'd5': '102245', 'd10': '102245', 'd20': '78391', 'd50': '4323',
            'd100': '36', 'd150': '36', 'd200': '36', 'd400': '36', 'd900': '36', 'd2000': '36', 'd2600': '2',
            'd2850': '2',
        },
        'hierCCv0': {
            'd0': '102245', 'd2': '102245', 'd5': '102245', 'd10': '102245', 'd20': '78391', 'd50': '4323',
            'd100': '36', 'd200': '36', 'd400': '36', 'd900': '36', 'd2000': '36', 'd2600': '2', 'd2850': '2',
        },
    }

    @patch('mist.scripts.mistlincode.restutils.retrieve_page_data')
    def test_extract_enterobase_with_hiercc_thresholds(self, mock_retrieve: Mock) -> None:
        """
        Tests extracting a partial LIN code for an EnteroBase-sourced profile, masked using hierCCv0-derived
        thresholds - 'hierCCv0' is preferred over the longer 'hierCC' (see ENTERO_INFO), and its 13 thresholds
        (including d0 as the deepest, 100%-similarity level) match the 13-position LIN code, so masking
        applies. Real-world example (ST 102245, senterica/cgMLST_v2, 3 mismatches): the threshold for the
        second-deepest position (2) is exceeded, so the two deepest positions are masked.
        :return: None
        """
        mock_retrieve.return_value.json.return_value = {'STs': [{'ST_id': '102245', 'info': self.ENTERO_INFO}]}
        profile = {'name': '102245', 'nb_matches': 2999, 'alleles': {f'locus{i}': '1' for i in range(3002)}}

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(
                dir_temp,
                downloader='enterobase',
                url='https://enterobase.warwick.ac.uk/schemes/Senterica.cgMLSTv2/',
            )
            result = MistLinCode(dir_db=dir_temp, entero_token='dummy-token').extract({'profiles': [profile]})

        self.assertEqual('102245', result['st'])
        self.assertEqual(['0', '0', '3', '0', '0', '0', '0', '3', '3', '7', '0', '0', '0'], result['lincode_full'])
        expected_partial = ['0', '0', '3', '0', '0', '0', '0', '3', '3', '7', '0', None, None]
        self.assertEqual(expected_partial, result['lincode_partial'])
        expected_thresholds = [2850, 2600, 2000, 900, 400, 200, 100, 50, 20, 10, 5, 2, 0]
        self.assertEqual(expected_thresholds, result['thresholds'])

        # Confirm the mocked call targeted the derived species/scheme slugs, authenticated with the token
        called_args, called_kwargs = mock_retrieve.call_args
        self.assertIn('/senterica/cgMLST_v2/sts', called_args[0])
        self.assertEqual(('dummy-token', ''), called_kwargs['auth'])

    @patch('mist.scripts.mistlincode.restutils.retrieve_page_data')
    def test_extract_enterobase_species_scheme_override(self, mock_retrieve: Mock) -> None:
        """
        Tests that explicit entero_species/entero_scheme override the (unreliable) auto-derived guess - e.g.
        E. coli needs 'ecoli'/'cgMLST', not the 'escherichia'/'cgMLST_v1' guessed from its download URL.
        :return: None
        """
        mock_retrieve.return_value.json.return_value = {'STs': [{'ST_id': '146609', 'info': self.ENTERO_INFO}]}
        profile = {'name': '146609', 'nb_matches': 3000, 'alleles': {f'locus{i}': '1' for i in range(3002)}}

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(
                dir_temp,
                downloader='enterobase',
                url='https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/',
            )
            MistLinCode(
                dir_db=dir_temp, entero_token='dummy-token', entero_species='ecoli', entero_scheme='cgMLST'
            ).extract({'profiles': [profile]})

        called_url = mock_retrieve.call_args[0][0]
        self.assertIn('/ecoli/cgMLST/sts', called_url)

    @patch('mist.scripts.mistlincode.restutils.retrieve_page_data')
    def test_extract_enterobase_falls_back_to_hiercc_without_v0(self, mock_retrieve: Mock) -> None:
        """
        Tests the real E. coli case: no 'hierCCv0' field at all (unlike Salmonella), so 'hierCC' is used
        directly. E. coli's 'hierCC' has 13 raw keys already matching EnteroBase's documented 13-level
        mapping (no extra level like Salmonella's 'd150'), so once d0 is kept as the deepest threshold, it
        lines up with the 13-position LIN code and masking applies - this used to always fall back to
        unmasked before d0 was included. Real-world example (ST 268260, ecoli/cgMLST, 11 mismatches).
        :return: None
        """
        info = {
            'lin_code': '0-0-0-0-2-89-0-0-4-0-2-1-36',
            'hierCC': {
                'd0': '268260', 'd2': '151557', 'd5': '142807', 'd10': '109409', 'd20': '109409',
                'd50': '71026', 'd100': '71026', 'd200': '71026', 'd400': '37', 'd1100': '13', 'd1500': '5',
                'd2000': '2', 'd2350': '1',
            },
        }
        mock_retrieve.return_value.json.return_value = {'STs': [{'ST_id': '268260', 'info': info}]}
        profile = {'name': '268260', 'nb_matches': 2502, 'alleles': {f'locus{i}': '1' for i in range(2513)}}

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(
                dir_temp,
                downloader='enterobase',
                url='https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/',
            )
            result = MistLinCode(
                dir_db=dir_temp, entero_token='dummy-token', entero_species='ecoli', entero_scheme='cgMLST'
            ).extract({'profiles': [profile]})

        expected_full = ['0', '0', '0', '0', '2', '89', '0', '0', '4', '0', '2', '1', '36']
        self.assertEqual(expected_full, result['lincode_full'])
        expected_partial = ['0', '0', '0', '0', '2', '89', '0', '0', '4', None, None, None, None]
        self.assertEqual(expected_partial, result['lincode_partial'])
        expected_thresholds = [2350, 2000, 1500, 1100, 400, 200, 100, 50, 20, 10, 5, 2, 0]
        self.assertEqual(expected_thresholds, result['thresholds'])

    @patch('mist.scripts.mistlincode.restutils.retrieve_page_data')
    def test_extract_enterobase_without_matching_hiercc_falls_back_unmasked(self, mock_retrieve: Mock) -> None:
        """
        Tests that a hierCC/LIN-code length mismatch falls back to an unmasked LIN code, rather than masking
        with a threshold count that doesn't correspond to the actual LIN-code positions.
        :return: None
        """
        info = {'lin_code': '0-0-3-0-0-0-0-3-3-7-0-0-0', 'hierCC': {'d0': '102245', 'd2': '102245'}}
        mock_retrieve.return_value.json.return_value = {'STs': [{'ST_id': '102245', 'info': info}]}
        profile = {'name': '102245', 'nb_matches': 3000, 'alleles': {f'locus{i}': '1' for i in range(3002)}}

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(
                dir_temp,
                downloader='enterobase',
                url='https://enterobase.warwick.ac.uk/schemes/Senterica.cgMLSTv2/',
            )
            result = MistLinCode(dir_db=dir_temp, entero_token='dummy-token').extract({'profiles': [profile]})

        self.assertIsNone(result['lincode_partial'])
        self.assertEqual([2, 0], result['thresholds'])

    def test_extract_enterobase_without_token_raises(self) -> None:
        """
        Tests that extracting a LIN code for an EnteroBase-sourced profile without a token raises clearly.
        :return: None
        """
        profile = {'name': '102245', 'nb_matches': 3000, 'alleles': {f'locus{i}': '1' for i in range(3002)}}

        with testingutils.get_temp_dir() as dir_temp:
            dir_temp = Path(dir_temp)
            self._write_db_info(
                dir_temp,
                downloader='enterobase',
                url='https://enterobase.warwick.ac.uk/schemes/Senterica.cgMLSTv2/',
            )
            with self.assertRaises(LinCodeError):
                MistLinCode(dir_db=dir_temp).extract({'profiles': [profile]})