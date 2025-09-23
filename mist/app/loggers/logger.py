import logging
from pathlib import Path

logger = logging.getLogger('MiST')


def initialize_logging(log_path: Path | None = None, debug: bool = False) -> None:
    """
    Initializes the logging.
    :param log_path: Path to store the log
    :param debug: If true, enable debug mode
    :return: None
    """
    formatter = logging.Formatter('%(asctime)s - %(module)15s - %(levelname)7s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO if not debug else logging.DEBUG)
    console_handler.name = 'console'
    logger.addHandler(console_handler)

    # File handler (pacu.log file)
    if log_path is not None:
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    # General logging level
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
