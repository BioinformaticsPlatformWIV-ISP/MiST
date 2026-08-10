from pathlib import Path


def count_rows(path: Path) -> int:
    """
    Counts the number of data rows (excluding the header) in a TSV file.
    :param path: Path to the TSV file
    :return: Number of data rows
    """
    nb_lines = 0
    with open(path, 'rb') as handle:
        while True:
            buf = handle.read(1 << 20)
            if not buf:
                break
            nb_lines += buf.count(b'\n')
    return max(nb_lines - 1, 0)


def parse_header(path: Path) -> list[str]:
    """
    Returns the columns for the TSV header.
    :param path: Path
    :return: Header
    """
    with path.open() as handle:
        return handle.readline().strip().split('\t')
