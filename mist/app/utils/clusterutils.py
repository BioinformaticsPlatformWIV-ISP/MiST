import re
from pathlib import Path
from typing import TypedDict

from mist.app.loggers.logger import logger
from mist.app.utils.command import Command


class ClusterDict(TypedDict):
    """
    Holds cluster information.
    """
    name: str
    id: int
    members: list[dict[str, str]]


def cluster_fasta(path_fasta: Path, dir_out: Path, cutoff: int, threads: int = 4, debug: bool = False) -> Path:
    """
    Clusters the input FASTA file.
    :param path_fasta: Input FASTA file
    :param dir_out: Output directory
    :param cutoff: Clustering cutoff
    :param threads: Number of threads
    :param debug: If True, the command is logged
    :return: Path to the clustered FASTA file
    """
    path_out = dir_out / path_fasta.name.replace(".fasta", "-clustered.fasta")
    command = Command(
        ' '.join([
            'cd-hit-est',
            '-i', str(path_fasta),
            '-o', str(path_out),
            '-S 0',
            '-d 0',
            '-c', str(cutoff / 100),
            '-T', str(threads)])
    )
    command.run(dir_out, disable_logging=not debug)
    if not command.exit_code == 0:
        raise ValueError(f"error running clustering: {command.stderr}")
    return path_out


def parse_cluster_from_file(path_fasta: Path) -> list[ClusterDict]:
    """
    Parses a CD-HIT .clstr file and returns a list of Cluster objects.
    :param path_fasta: Path to the CD-hit clustered FASTA file
    :return: List of clusters
    """
    # Check for the file with clusters
    path_clusters = path_fasta.parent / f'{path_fasta.name}.clstr'
    if not path_clusters.exists():
        raise FileNotFoundError(f'Clusters file not found: {path_clusters}')

    clusters = []
    re_cluster = re.compile(r'\d.+, >(.+)\.\.\. (at [+-])?.*$')

    current_cluster = None
    with path_clusters.open() as handle:
        for current_line in handle.readlines():
            line = current_line.strip()
            if line.startswith('>'):
                try:
                    name = line[1:].replace(' ', '_')
                    cluster_id = int(name.split('_')[-1])
                except (ValueError, IndexError):
                    raise ValueError(f'Invalid cluster header: {line}')
                current_cluster = {'name': name, 'id': cluster_id, 'members': []}
                clusters.append(current_cluster)
            else:
                if current_cluster is None:
                    raise ValueError(f'Sequence line found before any cluster header: {line}')
                match = re_cluster.match(line)
                if not match:
                    raise ValueError(f'Invalid cluster line: {line}')
                current_cluster['members'].append({
                    'seq_id': match.group(1), 'ori': match.group(2)[-1] if match.group(2) else None})
    logger.debug(f'{len(clusters):,} clusters parsed from {path_fasta.name}')
    return clusters
