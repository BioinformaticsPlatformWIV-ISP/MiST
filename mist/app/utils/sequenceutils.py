import hashlib
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

REGEX_ALLELE = r'_(\d+)$'


def rev_complement(seq_in: str) -> str:
    """
    Returns the reverse complement of the input sequence.
    :param seq_in: Input sequence
    :return: Reverse complement
    """
    return str(Seq(seq_in).reverse_complement())


def hash_sequence(seq_in: str, normalize: bool = False) -> str:
    """
    Hashes the input sequence.
    :param seq_in: Input sequence
    :param normalize: If true, the hash is calculated on the canonical form (min of seq and its reverse complement)
    :return: Hashed sequence
    """
    if not normalize:
        return hashlib.sha1(seq_in.replace('-', '').encode('ascii')).hexdigest()
    seq_in = seq_in.replace('-', '')
    seq_rc = rev_complement(seq_in)
    seq_canonical = min(seq_in, seq_rc)
    return hashlib.sha1(seq_canonical.encode('ascii')).hexdigest()


def extract_hashes(path_fasta: Path, ori_by_seq_id: dict[str, str]) -> pd.DataFrame:
    """
    Extracts the hashes from the target FASTA file.
    :param path_fasta: Input FASTA file
    :param ori_by_seq_id: The sequence orientation by sequence id
    :return: DataFrame with hashes
    """
    records_out = []
    with path_fasta.open() as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if ori_by_seq_id[record.id] in ('+', None):
                seq = str(record.seq)
            else:
                seq = rev_complement(record.seq)
            records_out.append({
                'seq_id': record.id,
                'hash': hash_sequence(seq),
            })
    data_out = pd.DataFrame(records_out)
    data_out['allele_id'] = data_out['seq_id'].apply(lambda x: re.search(REGEX_ALLELE, x).group(1))
    return data_out


def merge_fasta_files(paths_fasta: list[Path], path_out: Path) -> int:
    """
    Creates a FASTA file that combines all representative alleles.
    :param paths_fasta: Input FASTA files
    :param path_out: Output path
    :return: Number of sequences
    """
    sequences_added = 0
    with open(path_out, 'w') as handle_out:
        for path_fasta in paths_fasta:
            with path_fasta.open() as handle_in:
                for seq in SeqIO.parse(handle_in, 'fasta'):
                    SeqIO.write(seq, handle_out, 'fasta')
                    sequences_added += 1
    return sequences_added


def count_sequences(path_in: Path) -> int:
    """
    Counts the number of sequences in the input FASTA file.
    :param path_in: Input path
    :return: Number of sequences
    """
    with path_in.open() as handle:
        return len(list(SeqIO.parse(handle, 'fasta')))
