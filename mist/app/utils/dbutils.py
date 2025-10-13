import re
from pathlib import Path

from Bio import SeqIO

from mist.app.loggers.logger import logger


def is_valid_db(dir_in: Path) -> bool:
    """
    Checks if the input directory is a valid typing database.
    :param dir_in: Input directory
    :return: True if valid
    """
    if not (dir_in / 'loci.txt').exists():
        raise FileNotFoundError("'loci.txt' file not found")
    if not (dir_in / 'loci_repr.fasta').exists():
        raise FileNotFoundError("'loci_repr.fasta' file not found")
    return True

def _get_allele_id(seq_record: SeqIO.SeqRecord, locus_name: str, requires_match: bool = False) -> str:
    """
    Returns the allele id from the input SeqRecord
    :param seq_record: Sequence record
    :param locus_name: Locus name
    :param requires_match: If True, the sequence id must match the locus name
    :return: Allele id
    """
    if re.match(r'\d+', seq_record.id):
        return f'{locus_name}_{seq_record.id}'
    if re.match(f'{locus_name}_\\d+', seq_record.id, flags=re.IGNORECASE):
        return seq_record.id
    m = re.match(r'.+_(\d+)$', seq_record.id)
    if not requires_match and (m is not None):
        return f'{locus_name}_{m.group(1)}'
    raise ValueError(f'Invalid sequencing input FASTA file: {seq_record.id}')


def reformat_fasta(fasta_in: Path, fasta_out: Path) -> None:
    """
    Reforms the input FASTA file to the standardized format required by the tool.
    :param fasta_in: Input FASTA file
    :param fasta_out: Output FASTA file
    :return: None
    """
    if not fasta_in.name.endswith('.fasta'):
        raise ValueError("Input FASTA files should have the '.fasta' extension")
    locus_name = fasta_in.name.replace('.fasta', '')
    seq_ids_match_name = True
    with open(fasta_in) as handle_in, open(fasta_out, 'w') as handle_out:
        for seq in SeqIO.parse(handle_in, 'fasta'):
            try:
                allele_id = _get_allele_id(seq, locus_name, requires_match=True)
            except ValueError:
                seq_ids_match_name = False
                allele_id = _get_allele_id(seq, locus_name, requires_match=False)
            seq_out = SeqIO.SeqRecord(
                id=allele_id,
                description='',
                seq=seq.seq
            )
            SeqIO.write(seq_out, handle_out, 'fasta')
    if not seq_ids_match_name:
        logger.warning(
            f"FASTA file '{fasta_in.name}' contains sequences that do not match the locus name '{locus_name}'.")
