from pathlib import Path

from Bio import SeqIO

from mist.app.loggers.logger import logger
from mist.app.utils import sequenceutils


class SeqHolder:
    """
    Class to store and retrieve sequences from the input FASTA file.
    """

    def __init__(self, path_fasta: Path | None = None, seqs: list[str] | None = None) -> None:
        """
        Initializes the holder.
        :param path_fasta: Input FASTA path
        :param seqs: List of sequences
        :return: None
        """
        if path_fasta is not None:
            with path_fasta.open() as handle:
                self._seq_by_id = {s.id: str(s.seq) for s in SeqIO.parse(handle, 'fasta')}
        else:
            self._seq_by_id = {str(i): seq for i, seq in enumerate(seqs)}
        logger.debug(f'Parsed {len(self._seq_by_id):,} sequence(s)')

    def get_seq(self, seq_id: str, start: int, end: int, strand: str | None = None) -> str:
        """
        Retrieves the target sequence.
        :param seq_id: Sequence id
        :param start: Sequence start position
        :param end: Sequence end position
        :param strand: Strand (+ / -)
        :return: Sequence
        """
        if strand in ('+', None):
            return self._seq_by_id[seq_id][start-1:end]
        return sequenceutils.rev_complement(self._seq_by_id[seq_id][start-1:end])
