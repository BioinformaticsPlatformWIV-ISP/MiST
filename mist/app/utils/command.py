import subprocess
from pathlib import Path
from typing import Optional

from mist.app.loggers.logger import logger


class Command:
    """
    Class to handle the execution of commands
    """

    def __init__(self, command: str = None) -> None:
        """
        Initializes the command object.
        :param command: (optional) Command line call
        """
        self._stdout = None
        self._stderr = None
        self._procedure = None
        self._exit_code = None
        self._command = command

    @property
    def stderr(self) -> str:
        """
        Returns the stderr from the command execution.
        :return: Standard error
        """
        return self._stderr

    @property
    def stdout(self) -> str:
        """
        Returns the stdout from the command execution.
        :return: Standard error
        """
        return self._stdout

    @property
    def exit_code(self) -> int:
        """
        Returns the exit code from the command execution.
        :return: Exit code
        """
        return self._exit_code

    @property
    def command(self) -> str:
        """
        Returns the command line call.
        :return: Command line call
        """
        return self._command

    @command.setter
    def command(self, cmd: str) -> None:
        """
        Sets the command line call.
        :param cmd: Command
        :return: None
        """
        self._command = cmd

    def run(self, folder: Path, stderr_handle=subprocess.PIPE, disable_logging: bool = False,
            timeout: Optional[int] = None) -> None:
        """
        Runs the command given at command initialization
        :param folder: Folder where the command is executed
        :param stderr_handle: Handle for the standard error (e.g. PIPE or STDOUT)
        :param disable_logging: If True, logging is disabled
        :param timeout: Timeout for the command execution
        :return: None
        """
        if disable_logging is False:
            logger.info(f'Executing command: {self.command}')
        if self.command is None:
            raise ValueError("Invalid command 'None'")
        self._procedure = subprocess.run(
            self._command,
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
            shell=True,
            executable='/bin/bash',
            cwd=folder,
            timeout=timeout,
            text=True)
        self._stdout = self._procedure.stdout or ''
        self._stderr = self._procedure.stderr or ''
        self._exit_code = self._procedure.returncode
        if disable_logging is False:
            logger.debug(f'stdout: {self._stdout}')
            logger.debug(f'stderr: {self._stderr}')
