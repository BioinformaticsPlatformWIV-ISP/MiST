import unittest

from click.testing import CliRunner

from mist.app.loggers.logger import logger
from mist.scripts.cli import cli
from mist.version import __version__


class TestVersion(unittest.TestCase):
    """
    Tests the version command.
    """

    def test_version(self) -> None:
        """
        Tests the version command.
        :return: None
        """
        runner = CliRunner()
        # noinspection PyTypeChecker
        result = runner.invoke(cli,['--version'])
        logger.info(result.output)
        self.assertIn(__version__, result.output)
        self.assertTrue(result.exit_code == 0)
