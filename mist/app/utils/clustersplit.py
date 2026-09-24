import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO

from mist.app.loggers.logger import logger
from mist.app.utils import nucmerutils


class ClusterSplit:
    """
    This class is used to split clusters into groups of loci that do not overlap (i.e., have the same exact start and
    stop position).
    """

    def __init__(self, dir_in: Path, debug: bool) -> None:
        """
        Initializes this class.
        :param dir_in: Input directory
        :param debug: If True, enables debug mode
        :return: None
        """
        with open(dir_in / 'mist_db.json') as handle:
            self._metadata = json.load(handle)
        with open(dir_in / self._metadata['fasta_full']) as handle:
            self._seq_by_id = {s.id: s for s in SeqIO.parse(handle, 'fasta')}
        with open(dir_in / self._metadata['fasta_clustered']) as handle:
            self._seqs_clust = list(SeqIO.parse(handle, 'fasta'))
        self._dir_in = dir_in
        self._dir_temp = dir_in / 'tmp'
        self._dir_temp.mkdir(exist_ok=True)
        self._locus = self._metadata['name']
        self._debug = debug

    def run(self) -> None:
        """
        Runs the cluster splitting.
        :return: None
        """
        for _, cluster_data in self._metadata['clusters'].items():
            if len(cluster_data['members']) <= 1:
                continue
            self._seqs_clust.extend(self.process_cluster(cluster_data))
        path_out_fasta_clust_ext = self._dir_in / self._metadata['fasta_clustered']
        with open(path_out_fasta_clust_ext, 'w') as handle:
            SeqIO.write(self._seqs_clust, handle, 'fasta')
        shutil.rmtree(self._dir_temp)

    def process_cluster(self, cluster_data: dict[str, Any]) -> list[SeqIO.SeqRecord]:
        """
        Processes a cluster and returns the additional sequences that should be included in the clustered FASTA file.
        :param cluster_data: Cluster data
        :return: List of sequences to add
        """
        cluster_name = cluster_data['name']

        # Create FASTA with all sequences from the cluster
        path_fasta = self.export_seqs_from_cluster(cluster_data)
        path_fasta_ref = self.export_seqs_from_cluster(cluster_data, only_first=True)

        # Run nucmer
        path_nucmer = nucmerutils.nucmer(path_fasta_ref, path_fasta, self._dir_temp, debug=self._debug)
        data_coords = nucmerutils.show_coords(path_nucmer, debug=self._debug)

        # Select one sequence for each combination of start and end offsets that differs from the representative
        offsets_by_seq_id = ClusterSplit.calculate_offsets(data_coords)
        seq_id_by_offsets = {}
        for member in cluster_data['members']:
            offsets = offsets_by_seq_id.get(member['seq_id'])
            if (offsets is not None) and (offsets != (0, 0)):
                seq_id_by_offsets.setdefault(offsets, member['seq_id'])
        if len(seq_id_by_offsets) == 0:
            return []
        seq_ids = list(seq_id_by_offsets.values())
        logger.debug(f"Splitting cluster {cluster_name} ({self._locus}) into {len(seq_ids) + 1} groups")
        logger.debug(f"Novel representative(s): {seq_ids}")
        return [self._seq_by_id[seq_id] for seq_id in seq_ids]

    @staticmethod
    def calculate_offsets(data_coords: pd.DataFrame) -> dict[str, tuple[int, int]]:
        """
        Calculates the start and end offset of each sequence compared to the representative.
        :param data_coords: Parsed show-coords output
        :return: Offsets (start, end) by sequence id
        """
        data = data_coords.copy()

        # Use the coordinates on the reverse complement for alignments on the reverse strand
        is_rev = data['[S2]'] > data['[E2]']
        data['s2'] = data['[S2]'].where(~is_rev, data['[LEN Q]'] - data['[S2]'] + 1)
        data['e2'] = data['[E2]'].where(~is_rev, data['[LEN Q]'] - data['[E2]'] + 1)
        data['offset_start'] = data['s2'] - data['[S1]']
        data['offset_end'] = (data['[LEN Q]'] - data['e2']) - (data['[LEN R]'] - data['[E1]'])

        # Only retain alignments consistent with the longest alignment of each sequence
        data_main = data.sort_values('[LEN 2]', ascending=False).drop_duplicates('[TAG Q]')
        main = data_main.set_index('[TAG Q]').loc[data['[TAG Q]']].set_axis(data.index)
        is_before = (data['[S1]'] < main['[S1]']) & (data['s2'] < main['s2'])
        is_after = (data['[E1]'] > main['[E1]']) & (data['e2'] > main['e2'])
        is_main = data.index.isin(data_main.index)
        data = data[(is_rev == (main['[S2]'] > main['[E2]'])) & (is_before | is_after | is_main)]

        # Select the outermost alignments (ties: longest alignment)
        data_first = data.sort_values(['s2', '[LEN 2]'], ascending=[True, False]).drop_duplicates('[TAG Q]')
        data_last = data.sort_values(['e2', '[LEN 2]'], ascending=[False, False]).drop_duplicates('[TAG Q]')
        offset_end_by_seq_id = dict(zip(data_last['[TAG Q]'], data_last['offset_end']))
        return {
            str(seq_id): (int(offset_start), int(offset_end_by_seq_id[seq_id]))
            for seq_id, offset_start in zip(data_first['[TAG Q]'], data_first['offset_start'])
        }

    def export_seqs_from_cluster(self, cluster_data: dict[str, Any], only_first: bool = False) -> Path:
        """
        Export the sequences from the given cluster.
        :param cluster_data: Cluster data
        :param only_first: If True, only the first is selected
        :return: FASTA file with sequences
        """
        # Create output path
        if only_first:
            basename = f"_repr-{cluster_data['name']}"
        else:
            basename = f"_all-{cluster_data['name']}"
        path_out = Path(self._dir_temp, f'{basename}.fasta')

        # Export the FASTA file
        with path_out.open('w') as handle:
            if only_first:
                SeqIO.write(self._seq_by_id[cluster_data['members'][0]['seq_id']], handle, 'fasta')
            else:
                SeqIO.write([self._seq_by_id[member['seq_id']] for member in cluster_data['members']], handle, 'fasta')
        return path_out
