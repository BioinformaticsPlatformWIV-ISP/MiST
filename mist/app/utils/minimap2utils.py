from io import StringIO
from pathlib import Path
from typing import TextIO

import pandas as pd

from mist.app.utils.command import Command

COLS_PAF = {
    'qseqid': str,
    'qlen': int,
    'qstart': int,
    'qend': int,
    'sstrand': str,
    'sseqid': str,
    'slen': int,
    'sstart': int,
    'send': int,
    'bases_matching': int,
    'bases_total': int,
    'mq': int,
}


def create_index(path_fasta: Path) -> Path:
    """
    Creates a Minimap2 index for the input FASTA file.
    :param path_fasta: Input FASTA file
    :return: None
    """
    path_out = path_fasta.parent / f'{path_fasta.name}.mni'
    command = Command(' '.join([
        'minimap2',
        '-d', str(path_out),
        str(path_fasta)
    ]))
    command.run(path_fasta.parent, disable_logging=False)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error creating Minimap2 index: {command.stderr}')
    return path_out

def align(path_fasta: Path, path_fasta_db: Path, include_cigar: bool = False, threads: int = 1) -> pd.DataFrame:
    """
    Runs Minimap2 on the input sequence and DB.
    :param path_fasta: Input FASTA path
    :param path_fasta_db: Database FASTA path
    :param include_cigar: Include cigar string in output
    :param threads: Number of threads to use
    :return: Output results as DataFrame
    """
    command = Command(' '.join([
        'minimap2',
        str(path_fasta_db),
        str(path_fasta),
        '-O 25 -E 1', # ungapped alignment
        '-N 10 --all-chain', # include secondary alignments
        '-t', str(threads) # nb. of threads
    ] + (['-c'] if include_cigar else [])))
    command.run(Path().cwd(), disable_logging=False)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running Minimap2: {command.stderr}')
    return parse_paf(StringIO(command.stdout))


def parse_paf(file_in: Path | TextIO) -> pd.DataFrame:
    """
    Parses the input PAF file.
    :param file_in: Input path (or file-like object)
    :return: Parsed data as DataFrame
    """
    records_out = []
    offset = len(COLS_PAF)

    handle = file_in.open() if isinstance(file_in, Path) else file_in
    with handle:
        for row in handle:
            parts = row.rstrip('\n').split('\t')
            record = {col: type_(parts[i]) for i, (col, type_) in enumerate(COLS_PAF.items())}
            for tag in parts[offset:]:
                tag_name = tag.split(':', 1)[0]
                record[f"tag_{tag_name}"] = tag
            records_out.append(record)
    return pd.DataFrame(records_out)
