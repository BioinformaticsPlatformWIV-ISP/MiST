from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO

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

    def __init__(self, dir_in: Path) -> None:
        """
        Initializes the detector.
        :param dir_in: Input directory
        """
        self._dir_in = dir_in
        # Parse the database sequences, grouped by length (lowercase, consistent with sequence hashing)
        self._seqs_by_length: dict[int, list[tuple[str, str]]] = defaultdict(list)
        with (self._dir_in / f'{self._dir_in.name}.fasta').open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                self._seqs_by_length[len(seq)].append((seq.id, str(seq.seq).lower()))
        nb_seqs = sum(len(seqs) for seqs in self._seqs_by_length.values())
        logger.debug(f'Parsed: {nb_seqs:,} sequences ({dir_in.name})')

    def retrieve_best_matching(self, seq: str, min_id: int) -> list[str]:
        """
        Retrieves the best matching sequence(s) for the target sequence.
        :param seq: Target sequence
        :param min_id: Min. % sequence identity
        :return: Seq ids for the best matching sequences
        """
        candidates = self._seqs_by_length.get(len(seq), [])
        logger.debug(f'Found {len(candidates):,} allele(s) matching the length of the detected sequence ({len(seq)})')
        if len(candidates) == 0:
            viable_lengths = list(self._seqs_by_length.keys())
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
