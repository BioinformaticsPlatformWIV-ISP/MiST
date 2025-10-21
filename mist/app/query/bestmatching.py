from pathlib import Path

import pandas as pd
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
        # Parse the database sequences
        records_in = []
        with (self._dir_in / f'{self._dir_in.name}.fasta').open() as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                records_in.append({
                    'id': seq.id,
                    'seq': str(seq.seq),
                    'length': len(seq)
                })
        self._data_seqs_db = pd.DataFrame(records_in)
        logger.debug(f'Parsed: {len(self._data_seqs_db):,} sequences ({dir_in.name})')

    def retrieve_best_matching(self, seq: str, min_id: int) -> list[str]:
        """
        Retrieves the best matching sequence(s) for the target sequence.
        :param seq: Target sequence
        :param min_id: Min. % sequence identity
        :return: Seq ids for the best matching sequences
        """
        data_subset = self._data_seqs_db[self._data_seqs_db['length'] == len(seq)].copy()
        logger.info(f'Found {len(data_subset):,} allele(s) matching the length of the detected sequence ({len(seq)})')
        if len(data_subset) == 0:
            viable_lengths = list(self._data_seqs_db['length'].unique())
            logger.debug(
                f"Length of detected sequence ({len(seq):,}) does not match any alleles in the "
                f"database ({', '.join(str(l) for l in viable_lengths)})")
            raise InvalidLengthException(len(seq), viable_lengths)
        data_subset['nb_matches'] = data_subset['seq'].apply(lambda x: sum(c1 == c2 for c1, c2 in zip(x, seq)))
        max_matches = data_subset['nb_matches'].max()

        # Check if the identity matches
        identity = 100 * max_matches / len(seq)
        if identity <= min_id:
            logger.info(f'Identity ({identity:.2f}%) to best matching sequence is below threshold ({min_id}%).')
            return []
        return list(data_subset[data_subset['nb_matches'] == max_matches]['id'])
