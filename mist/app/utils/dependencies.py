from pathlib import Path

from mist.app.errors import DependencyError
from mist.app.loggers.logger import logger
from mist.app.utils.command import Command

VERSION_COMMANDS = {
    'CD-HIT': {'command': 'cd-hit -h', 'exit_code': 1},
    'minimap2': {'command': 'minimap2 --version', 'exit_code': 0},
    'nucmer': {'command': 'nucmer --version', 'exit_code': 0},
}

def _check_dependency(key: str) -> None:
    """
    Checks if the given dependency is available.
    :param key: Dependency key
    :return: None
    """
    if key not in VERSION_COMMANDS:
        raise ValueError(f'Unknown dependency key: {key}')
    command = Command(VERSION_COMMANDS[key]['command'])
    command.run(Path.cwd(), disable_logging=True)
    if command.exit_code != VERSION_COMMANDS[key]['exit_code']:
        raise DependencyError(f"Dependency '{key}' is not available")

def check_dependencies(keys: list[str]) -> None:
    """
    Checks if the given dependencies are available.
    :param keys: Keys to check
    :return: None
    """
    logger.info('Checking if dependencies are available')

    for key in keys:
        _check_dependency(key)
        logger.info(f'  {key}: OK')
