from pathlib import Path

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

from mist.app.loggers.logger import logger


class InvalidLengthException(Exception):
    """
    Error that is raised when the length of the input sequence does not match any sequences in the database.
    """

    def __init__(self, len_seq: int, allowed: list[int]) -> None:
        """
        Initializes the exception.
        :param len_seq: Sequence length
        :param allowed: Allowed sequence lengths
        """
        self.len_seq = len_seq
        self.allowed = allowed


def retrieve_best_matching(dir_in: Path, seq: str, min_id: int) -> list[str]:
    """
    Retrieves the best matching allele(s) for the target sequence, only alleles with the same length are considered.
    :param dir_in: Locus directory
    :param seq: Target sequence
    :param min_id: Min. % sequence identity
    :return: Seq ids for the best matching alleles
    """
    # Parse the alleles with the same length as the target sequence (lowercase, consistent with sequence hashing)
    candidates: list[tuple[str, str]] = []
    lengths: set[int] = set()
    with (dir_in / f'{dir_in.name}.fasta').open() as handle:
        for title, seq_allele in SimpleFastaParser(handle):
            lengths.add(len(seq_allele))
            if len(seq_allele) == len(seq):
                candidates.append((title.split(None, 1)[0], seq_allele.lower()))
    logger.debug(f'Found {len(candidates):,} allele(s) matching the length of the detected sequence ({len(seq)})')
    if len(candidates) == 0:
        viable_lengths = list(lengths)
        logger.debug(
            f"Length of detected sequence ({len(seq):,}) does not match any alleles in the "
            f"database ({', '.join(str(l) for l in sorted(viable_lengths))})"
        )
        raise InvalidLengthException(len(seq), viable_lengths)

    # Count the matching positions for all candidates at once (one row per candidate)
    seq_ids, seqs = zip(*candidates)
    matrix = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8).reshape(len(seqs), len(seq))
    nb_matches = (matrix == np.frombuffer(seq.lower().encode('ascii'), dtype=np.uint8)).sum(axis=1)
    max_matches = nb_matches.max()

    # Check if the identity matches
    identity = 100 * max_matches / len(seq)
    if identity <= min_id:
        logger.debug(f'Identity ({identity:.2f}%) to best matching sequence is below threshold ({min_id}%).')
        return []
    return [seq_id for seq_id, nb in zip(seq_ids, nb_matches) if nb == max_matches]
