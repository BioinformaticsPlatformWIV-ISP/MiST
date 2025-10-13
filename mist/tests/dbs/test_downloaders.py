import os
import unittest
from pathlib import Path

from mist.app.dbs.bigsdbauthdownloader import BIGSDbAuthDownloader
from mist.app.dbs.bigsdbdownloader import BIGSDbDownloader
from mist.app.dbs.cgmlstorgdownloader import CgMLSTOrgDownloader
from mist.app.dbs.enterobasedownloader import EnteroBaseDownloader
from mist.app.loggers.logger import initialize_logging
from mist.app.utils import testingutils


class TestDownloaders(unittest.TestCase):
    """
    Tests for the downloading MLST schemes.
    """

    @staticmethod
    def download_succeeded(dir_scheme: Path, with_profiles: bool) -> bool:
        """
        Checks if the scheme in the target directory is valid.
        :param dir_scheme: Directory to check
        :param with_profiles: If true, profiles are checked
        :return: True if valid
        """
        path_fasta_list = dir_scheme / 'fasta_list.txt'
        if not path_fasta_list.exists():
            raise FileNotFoundError(f'{path_fasta_list.name} is missing')
        with path_fasta_list.open() as handle:
            for p in [Path(x.strip()) for x in handle.readlines()]:
                if not p.exists():
                    raise FileNotFoundError(f'{p.name} is missing')
        if with_profiles:
            if not (dir_scheme / 'profiles.tsv').exists():
                raise FileNotFoundError('profiles.tsv file is missing')
        return True

    def setUp(self) -> None:
        """
        Sets up a temporary directory and builds a database there before each test.
        :return: None
        """
        self.dir_temp = testingutils.get_temp_dir()

    def tearDown(self) -> None:
        """
        Clean up the temporary directory after the test.
        :return: None
        """
        self.dir_temp.cleanup()

    def test_download_cgmlst_org(self) -> None:
        """
        Tests downloading a scheme from cgMLST.org.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            downloader = CgMLSTOrgDownloader()
            downloader.download_scheme(
                url='https://www.cgmlst.org/ncs/schema/Cjejuni22/',
                dir_out=Path(dir_temp),
                include_profiles=True)
            self.assertTrue(TestDownloaders.download_succeeded(Path(dir_temp), with_profiles=True))

    def test_download_enterobase(self) -> None:
        """
        Tests downloading a scheme from EnteroBase.
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            downloader = EnteroBaseDownloader()
            downloader.download_scheme(
                url='https://enterobase.warwick.ac.uk/schemes/Yersinia.Achtman7GeneMLST/',
                dir_out=Path(dir_temp),
                include_profiles=True)
            self.assertTrue(TestDownloaders.download_succeeded(Path(dir_temp), with_profiles=True))

    def test_download_pubmlst(self) -> None:
        """
        Tests downloading a scheme from PubMLST (without authentication).
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            downloader = BIGSDbDownloader()
            downloader.download_scheme(
                url='https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/1',
                dir_out=Path(dir_temp),
                include_profiles=True)
            self.assertTrue(TestDownloaders.download_succeeded(Path(dir_temp), with_profiles=True))

    @unittest.skipUnless(os.environ.get('MIST_AUTH'), "Skipping test (MIST_AUTH not set).")
    def test_download_pubmlst_with_auth(self) -> None:
        """
        Tests downloading a scheme from PubMLST (with authentication).
        :return: None
        """
        with testingutils.get_temp_dir() as dir_temp:
            downloader = BIGSDbAuthDownloader(
                site='PubMLST',
                key_name='PubMLST',
                dir_tokens=Path(os.environ.get('MIST_AUTH'))
            )
            downloader.download_scheme(
                url='https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/1',
                dir_out=Path(dir_temp),
                include_profiles=True
            )
            self.assertTrue(TestDownloaders.download_succeeded(Path(dir_temp), with_profiles=True))


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
