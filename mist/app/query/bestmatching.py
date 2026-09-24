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


class ImperfectMatchDetector:
    """
    Identifies the best matching imperfect hit.
    """

    def __init__(self, dir_in: Path, seq_length: int) -> None:
        """
        Initializes the detector.
        :param dir_in: Input directory
        :param seq_length: Length of the target sequence (only alleles with this length are retained)
        """
        self._dir_in = dir_in
        self._seq_length = seq_length
        # Parse the database sequences with the target length (lowercase)
        self._candidates: list[tuple[str, str]] = []
        self._lengths: set[int] = set()
        with (self._dir_in / f'{self._dir_in.name}.fasta').open() as handle:
            for title, seq in SimpleFastaParser(handle):
                self._lengths.add(len(seq))
                if len(seq) == seq_length:
                    self._candidates.append((title.split(None, 1)[0], seq.lower()))
        logger.debug(f'Parsed: {len(self._candidates):,} sequences with length {seq_length:,} ({dir_in.name})')

    def retrieve_best_matching(self, seq: str, min_id: int) -> list[str]:
        """
        Retrieves the best matching sequence(s) for the target sequence.
        :param seq: Target sequence
        :param min_id: Min. % sequence identity
        :return: Seq ids for the best matching sequences
        """
        if len(seq) != self._seq_length:
            raise ValueError(f'Sequence length ({len(seq):,}) does not match the detector ({self._seq_length:,})')
        candidates = self._candidates
        logger.debug(f'Found {len(candidates):,} allele(s) matching the length of the detected sequence ({len(seq)})')
        if len(candidates) == 0:
            viable_lengths = list(self._lengths)
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
