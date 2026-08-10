import concurrent.futures
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from mist.app import model
from mist.app.loggers.logger import logger
from mist.app.utils import alleleutils, tsvutils

NAME_ALLELE_CODES = 'allele_codes.npy'
NAME_PROFILE_IDS = 'profile_ids.npy'
NAME_META = 'index_meta.pkl'


class ProfileIndex:
    """
    Class to query profiles files with a given allele combination, using an inverted index for speed.

    Groups profile indices by (locus, allele) at construction time, so a query only has to look up the profiles that
    share a detected allele at each locus, instead of comparing every profile against every queried locus.

    Alleles are stored internally as compact int codes rather than strings.

    Both the allele codes and the index itself are backed by on-disk numpy memmaps rather than held fully in RAM
    (which would be infeasible for huge schemes).
    """

    CHUNK_SIZE = 25_000
    THREADS = os.cpu_count() or 1

    def __init__(
        self,
        path_profiles: Path | None = None,
        loci: list[str] | None = None,
        dir_out: Path | None = None,
        dir_index: Path | None = None,
    ) -> None:
        """
        Initializes the profile index, either by parsing a profiles file and building the index from scratch
        (written directly to `dir_out`), or by reloading a previously built one from `dir_index`.
        :param path_profiles: Path to profiles file (builds the index from scratch, together with `loci`/`dir_out`)
        :param loci: List of loci (builds the index from scratch, together with `path_profiles`/`dir_out`)
        :param dir_out: Directory to build the index into (builds the index from scratch)
        :param dir_index: Directory containing a previously built index, reloaded instead of rebuilt
        :return: None
        """
        if dir_index is not None:
            dir_data = dir_index
            with open(dir_data / NAME_META, 'rb') as handle:
                meta = pickle.load(handle)
            self._loci: list[str] = meta['loci']
            self._names: list[str] = meta['names']
            self._metadata: list[list[tuple[str, str]]] = meta['metadata']
            self._buckets: dict[str, dict[int, tuple[int, int]]] = meta['buckets']
            logger.debug(f'Loaded index: {len(self._names):,} profiles ({dir_data})')
        else:
            dir_data = dir_out
            dir_data.mkdir(parents=True, exist_ok=True)
            self._loci, self._names, self._metadata, counts_by_locus_code = self._parse_profiles(
                path_profiles, set(loci), dir_data
            )
            self._buckets = self._scatter_profile_ids(dir_data, self._loci, counts_by_locus_code)
            with open(dir_data / NAME_META, 'wb') as handle:
                pickle.dump(
                    {
                        'loci': self._loci,
                        'names': self._names,
                        'metadata': self._metadata,
                        'buckets': self._buckets
                    }, handle,
                )
            logger.info(f'Profile index built: {len(self._names):,} profiles ({dir_data})')

        self._locus_to_col: dict[str, int] = {locus: i for i, locus in enumerate(self._loci)}
        self._allele_codes = np.load(dir_data / NAME_ALLELE_CODES, mmap_mode='r')
        self._profile_ids = np.load(dir_data / NAME_PROFILE_IDS, mmap_mode='r')

    def _parse_profiles(
        self, path: Path, locus_names: set[str], dir_out: Path
    ) -> tuple[list[str], list[str], list[list[tuple[str, str]]], list[dict[int, int]]]:
        """
        Parses the profile file in chunks of `CHUNK_SIZE` rows, writing the encoded allele codes directly into an
        `allele_codes.npy` memmap in `dir_out` as they are read. Real cgMLST profiles files can be too large to hold as
        a string DataFrame (or even as plain Python ints) all at once.
        :param path: Path to the TSV file
        :param locus_names: Locus names
        :param dir_out: Output directory to write `allele_codes.npy` into
        :return: Locus names (in file column order), profile names, metadata, and a tally of how many times each
            allele code occurs at each locus (one dict per locus) - used to size the buckets afterward
        """
        nb_rows = tsvutils.count_rows(path)
        logger.debug(f'Nb rows: {nb_rows:,}')

        cols_header = tsvutils.parse_header(path)
        cols_metadata: list[str] = [col for col in cols_header if col not in locus_names]
        logger.debug(f'Metadata columns: {cols_metadata}')
        cols_alleles: list[str] = [col for col in cols_header if col in locus_names]
        logger.debug(f'Gene columns: {cols_alleles}')

        # Tallied to size the buckets afterward (see `_scatter_profile_ids`)
        counts_by_locus_code: list[dict[int, int]] = [{} for _ in cols_alleles]

        # A memmap needs its shape fixed at creation (requiring the nb_rows upfront)
        allele_codes = np.lib.format.open_memmap(
            dir_out / NAME_ALLELE_CODES, mode='w+', dtype=np.int32, shape=(nb_rows, len(cols_alleles))
        )

        names: list[str] = []
        metadata: list[list[tuple[str, str]]] = []

        nb_read = 0
        for chunk in pd.read_table(path, dtype=str, chunksize=self.CHUNK_SIZE):
            chunk = chunk.fillna('n/a')

            names.extend(chunk[chunk.columns[0]].tolist())
            metadata.extend(list(zip(cols_metadata, row)) for row in chunk[cols_metadata].to_numpy())

            codes_chunk = np.array(
                [[alleleutils.encode_allele(val) for val in row] for row in chunk[cols_alleles].to_numpy()],
                dtype=np.int32,
            )

            # Write this chunk into the memmap, and tally its codes per (locus, code) - vectorized via np.unique
            # rather than a per-cell loop, since this runs once per chunk over up to CHUNK_SIZE rows.
            allele_codes[nb_read : nb_read + chunk.shape[0], :] = codes_chunk
            for j, locus_counts in enumerate(counts_by_locus_code):
                codes, counts = np.unique(codes_chunk[:, j], return_counts=True)
                for code, count in zip(codes.tolist(), counts.tolist()):
                    locus_counts[code] = locus_counts.get(code, 0) + count

            nb_read += chunk.shape[0]
            logger.debug(f'Parsed {nb_read:,} / {nb_rows:,} profiles')

        allele_codes.flush()
        return cols_alleles, names, metadata, counts_by_locus_code

    @staticmethod
    def _scatter_profile_ids(
        dir_out: Path, loci: list[str], counts_by_locus_code: list[dict[int, int]]
    ) -> dict[str, dict[int, tuple[int, int]]]:
        """
        Groups profile indices by (locus, allele code) into a `profile_ids.npy` memmap in `dir_out`, so a query
        only needs to read the (start, count) slice of one locus's column that matches a detected allele, instead
        of scanning every profile. `profile_ids.npy` is stored column-major (Fortran order) so that slice is a
        contiguous disk read rather than one page per row. Reads `allele_codes.npy` back in chunks of `CHUNK_SIZE`
        rows to do the scattering, rather than holding it fully in memory again.

        Loci are scattered across a `THREADS`-sized thread pool, one batch of loci per thread rather than one
        task per locus (submitting a fine-grained task per locus would drown the actual work in thread-pool
        overhead). This is safe because each locus owns its own output column and cursor, and it's worth doing
        because `_scatter_locus` below is vectorized numpy, which releases the GIL - so threads here give real
        parallelism instead of fighting over it.
        :param dir_out: Output directory containing `allele_codes.npy` (input) and to write `profile_ids.npy` into
        :param loci: Locus names, aligned to the columns of `allele_codes.npy` and `counts_by_locus_code`
        :param counts_by_locus_code: Tally of how many times each allele code occurs at each locus (see
            `_parse_profiles`), used to size each bucket
        :return: Mapping of locus -> allele code -> (start, count) offset into that locus's column of
            `profile_ids.npy`
        """
        buckets: dict[str, dict[int, tuple[int, int]]] = {}
        cursors: list[dict[int, int]] = []
        for locus, locus_counts in zip(loci, counts_by_locus_code):
            offset = 0
            locus_bucket: dict[int, tuple[int, int]] = {}
            cursor: dict[int, int] = {}
            for code, count in locus_counts.items():
                locus_bucket[code] = (offset, count)
                cursor[code] = offset
                offset += count
            buckets[locus] = locus_bucket
            cursors.append(cursor)

        allele_codes = np.load(dir_out / NAME_ALLELE_CODES, mmap_mode='r')
        nb_rows = allele_codes.shape[0]
        profile_ids = np.lib.format.open_memmap(
            dir_out / NAME_PROFILE_IDS, mode='w+', dtype=np.int32, shape=(nb_rows, len(loci)), fortran_order=True
        )
        locus_batches = [batch.tolist() for batch in np.array_split(np.arange(len(loci)), ProfileIndex.THREADS)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=ProfileIndex.THREADS) as executor:
            for row_start in range(0, nb_rows, ProfileIndex.CHUNK_SIZE):
                row_end = min(row_start + ProfileIndex.CHUNK_SIZE, nb_rows)
                block = allele_codes[row_start:row_end, :]
                futures = [
                    executor.submit(ProfileIndex._scatter_locus_batch, batch, block, cursors, profile_ids, row_start)
                    for batch in locus_batches
                    if len(batch) > 0
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()
                logger.debug(f'Bucketed {row_end:,} / {nb_rows:,} profiles')
        profile_ids.flush()
        return buckets

    @staticmethod
    def _scatter_locus_batch(
        locus_indices: list[int],
        block: np.ndarray,
        cursors: list[dict[int, int]],
        profile_ids: np.memmap,
        row_start: int,
    ) -> None:
        """
        Scatters a batch of locus columns of one row-chunk into their bucket positions in `profile_ids`. Meant
        to be handed a batch of loci per thread (see `_scatter_profile_ids`), not called once per locus.
        :param locus_indices: Column indices (loci) to scatter, within this chunk/`profile_ids`
        :param block: This chunk's rows of `allele_codes` (rows x loci)
        :param cursors: Per-locus next-free-slot-per-code, updated in place
        :param profile_ids: Output memmap (rows x loci) to scatter row indices into
        :param row_start: Global row index of the first row in `block`
        :return: None
        """
        for j in locus_indices:
            ProfileIndex._scatter_locus(block[:, j], cursors[j], profile_ids[:, j], row_start)

    @staticmethod
    def _scatter_locus(codes: np.ndarray, cursor: dict[int, int], column: np.ndarray, row_start: int) -> None:
        """
        Scatters one locus column of one row-chunk into its bucket positions. Vectorized rather than looping
        over every row: a single sort groups rows by code, each row's position within its own group is computed
        for the whole chunk at once, and the whole chunk is written to `column` in one indexed assignment - the
        only per-row-group Python loop left is over the (few) distinct codes in this chunk, to advance the
        cursor.
        :param codes: This chunk's allele codes for one locus
        :param cursor: Next-free-slot-per-code for this locus, updated in place
        :param column: This locus's full column of `profile_ids`, written into at the computed positions
        :param row_start: Global row index of the first row in `codes`
        :return: None
        """
        # `kind='stable'` keeps rows within the same code group in their original (file) order.
        order = np.argsort(codes, kind='stable')
        global_rows_sorted = (row_start + order).astype(np.int32)
        sorted_codes = codes[order]
        unique_codes, group_starts, group_counts = np.unique(sorted_codes, return_index=True, return_counts=True)

        # Rank of each row within its own code-group (0-indexed), computed for the whole chunk at once.
        rank_within_group = np.arange(len(sorted_codes)) - np.repeat(group_starts, group_counts)
        base_per_code = np.fromiter(
            (cursor[int(code)] for code in unique_codes), dtype=np.int64, count=len(unique_codes)
        )
        dest_positions = np.repeat(base_per_code, group_counts) + rank_within_group

        for code, count in zip(unique_codes.tolist(), group_counts.tolist()):
            cursor[code] += count
        column[dest_positions] = global_rows_sorted

    def _build_profile(self, idx: int) -> model.Profile:
        """
        Builds a `Profile` object for the given row index, decoding its allele codes back into strings.
        :param idx: Row index
        :return: Profile
        """
        return model.Profile(
            name=self._names[idx],
            alleles={
                locus: alleleutils.decode_allele(int(code))
                for locus, code in zip(self._loci, self._allele_codes[idx, :])
            },
            metadata=self._metadata[idx],
        )

    def query(self, result_by_locus: dict[str, model.QueryResult | None]) -> tuple[list[model.Profile], int]:
        """
        Queries the index using the detected alleles.

        Common alleles (especially the wildcard, which every locus checks unconditionally) can match tens of
        thousands of profiles in a real scheme, so per-profile counting is done with `np.bincount` over all
        matched indices at once, rather than incrementing a `Counter` per profile - the latter turned a fast
        bucket lookup into a slow, unvectorized Python loop over every match.
        :param result_by_locus: Detected allele(s) by locus
        :return: Best matching profiles, nb. of matching loci
        """
        matched_indices: list[np.ndarray] = []
        for locus, res in result_by_locus.items():
            bucket = self._buckets[locus]
            col = self._locus_to_col[locus]
            # Profiles with a wildcard at this locus always match, regardless of the query
            for code in alleleutils.candidate_codes(res) | {alleleutils.CODE_WILDCARD}:
                if code not in bucket:
                    continue
                start, count = bucket[code]
                matched_indices.append(self._profile_ids[start : start + count, col])

        if not matched_indices:
            # No profile matched anything (incl. wildcards) at any queried locus - nothing useful to report.
            return [], -1 if len(self._names) == 0 else 0

        counts = np.bincount(np.concatenate(matched_indices), minlength=len(self._names))
        best_matches = int(counts.max())
        best_indices = np.flatnonzero(counts == best_matches)
        return [self._build_profile(int(idx)) for idx in best_indices], best_matches
