import time
from pathlib import Path

import requests

from mist.app.loggers.logger import logger
from mist.app.utils import sequenceutils

REQUEST_SLEEP_SECONDS = 1


def retrieve_page_data(url: str, retries: int = 3, timeout: int = 20) -> requests.Response:
    """
    Retrieves data from the given URL.
    :param url: URL
    :param retries: Number of retries
    :param timeout: Timeout (in seconds)
    :return: URL data
    """
    error = None
    for retry in range(retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as err:
            logger.warning(f'Error retrieving {url} ({err!r}), retrying ({retry + 1}/{retries})')
            error = err
            time.sleep(min(10, 1 * (2**retry)))
    raise RuntimeError(f'Failed to retrieve {url} after {retries} attempts: {type(error).__name__}: {error}')


def download_fasta(locus_name: str, url: str, dir_out: Path) -> Path:
    """
    Downloads the FASTA for the given locus.
    :param locus_name: Locus name
    :param url: URL
    :param dir_out: Output directory
    :return: FASTA filename
    """
    response = retrieve_page_data(url)
    path_out = dir_out / f'{locus_name}.fasta'
    with path_out.open('w', encoding='utf-8') as handle:
        handle.write(response.text)
    logger.debug(f'FASTA file downloaded: {path_out.name} ({sequenceutils.count_sequences(path_out)})')
    time.sleep(REQUEST_SLEEP_SECONDS)
    return path_out
