from io import StringIO
from pathlib import Path

import pandas as pd

from mist.app.utils.command import Command

# Columns of the show-coords output ('-rcl -T'), the [TAGS] column in the header spans two fields
COLS_COORDS = [
    '[S1]',
    '[E1]',
    '[S2]',
    '[E2]',
    '[LEN 1]',
    '[LEN 2]',
    '[% IDY]',
    '[LEN R]',
    '[LEN Q]',
    '[COV R]',
    '[COV Q]',
    '[TAG R]',
    '[TAG Q]',
]


def nucmer(path_fasta_ref: Path, path_fasta_in: Path, dir_out: Path, threads: int = 1, debug: bool = False) -> Path:
    """
    Runs Minimap2 on the input sequence and DB.
    :param path_fasta_ref: Reference sequence FASTA file
    :param path_fasta_in: Input FASTA path
    :param dir_out: Output directory
    :param threads: Number of threads to use
    :param debug: If True, enable debug mode
    :return: None
    """
    basename = path_fasta_in.name.replace('.fasta', '')
    path_out = dir_out / basename
    command = Command(
        ' '.join(
            [
                'nucmer',
                '--maxmatch',
                f'-p {path_out}',
                str(path_fasta_ref),
                str(path_fasta_in),
                '--threads',
                str(threads),  # nb. of threads
            ]
        )
    )
    command.run(Path().cwd(), disable_logging=not debug)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running nucmer: {command.stderr}')
    return path_out.parent / f'{basename}.delta'


def show_coords(path_tsv: Path, debug: bool = False) -> pd.DataFrame:
    """
    Extracts the coordinates from the delta file.
    :param path_tsv: Path to delta file
    :param debug: If True, enable debug mode
    :return: Parsed coordinates
    """
    command = Command(' '.join(['show-coords', '-rcl', '-T', str(path_tsv)]))
    command.run(Path().cwd(), disable_logging=not debug)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running show-coords: {command.stderr}')
    return pd.read_table(StringIO(command.stdout), skiprows=2, header=0, names=COLS_COORDS)
