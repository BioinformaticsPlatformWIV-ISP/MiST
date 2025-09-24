#!/usr/bin/env python
from datetime import datetime
from importlib.resources import files
from pathlib import Path
from typing import Callable

import click

from mist.app.dbs import DOWNLOADERS
from mist.app.loggers.logger import initialize_logging, logger
from mist.app.query.allelequeryminimap import MultiStrategy
from mist.scripts.mistcaller import MistCaller
from mist.scripts.mistdists import MistDists
from mist.scripts.mistdownload import MistDownload
from mist.scripts.mistindex import MistIndex
from mist.scripts.mistlist import MistList
from mist.version import __version__


def _common_options(func: Callable) -> Callable:
    """
    Defines the common options.
    :param func: Function to decorate
    :return: Decorated function
    """
    func = click.option("--log", type=Path, help="Save log to this file")(func)
    func = click.option("--debug", is_flag=True, help="Enable debug mode")(func)
    return func


@click.group(help='MiST: Minimap2-inferred Sequence Typing.')
@click.version_option(__version__, prog_name="MiST", message="%(prog)s %(version)s")
def cli() -> None:
    """
    MiST: Minimap2-inferred Sequence Typing.
    """
    pass

@cli.command(name='index')
@click.argument("fasta", nargs=-1, type=click.Path(exists=True))
@click.option("-l", "--fasta-list", type=click.Path(path_type=Path, exists=True), help="List with input FASTA path(s)")
@click.option("-p", "--profiles", type=click.Path(path_type=Path, exists=True), help="TSV file with profiles")
@click.option("-o", "--output", type=click.Path(path_type=Path), required=True, help="Output directory")
@click.option("-c", "--cutoff", type=int, default=95, show_default=True, help="Clustering cutoff")
@click.option("-t", "--threads", type=int, default=1, show_default=True, help="Number of threads to use")
@_common_options
def index_(
        fasta: list[Path],
        fasta_list: Path,
        profiles: Path | None,
        output: Path,
        cutoff: int,
        threads: int,
        debug: bool,
        log: Path) -> None:
    """
    Creates a MIST index for a set of input FASTA files.
    """
    initialize_logging(log_path=log, debug=debug)

    # Collect FASTA files
    if len(fasta) > 0:
        paths_fasta = [Path(f) for f in fasta]
    elif fasta_list is not None:
        with open(fasta_list) as handle:
            paths_fasta = [Path(line.strip()) for line in handle.readlines() if len(line.strip()) > 0]
    else:
        raise click.UsageError("No FASTA file or FASTA list provided.")
    if len(paths_fasta) < 1:
        raise click.UsageError("No input FASTA file(s) provided.")

    # Run the indexer
    indexer = MistIndex(
        paths_fasta=paths_fasta,
        path_profiles=profiles,
        cutoff=cutoff,
        debug=debug
    )
    indexer.create_index(dir_out=output.expanduser().resolve(), threads=threads)


@cli.command()
@click.option("-f", "--fasta", type=click.Path(exists=True, path_type=Path, dir_okay=False), required=True, help="Input FASTA path")
@click.option("-d", "--db", type=click.Path(exists=True, path_type=Path), required=True, help="Database path")
@click.option("-o", "--out-json", type=click.Path(), default=None, help="JSON output file")
@click.option("-t", "--threads", type=int, default=1, show_default=True, help="Number of threads to use")
@click.option("--out-tsv", type=click.Path(path_type=Path, writable=True), help="TSV output file")
@click.option("--out-dir", type=click.Path(path_type=Path), help="Output directory")
@click.option("--export-novel", is_flag=True, help="Create FASTA files for (potential) novel alleles")
@click.option("--keep-minimap2", is_flag=True, help="Store the minimap2 output")
@click.option("--min-id-novel", type=int, default=99, show_default=True, help="Minimum % identity for novel alleles")
@click.option(
    "-m", "--multi", type=click.Choice([s.value for s in MultiStrategy]), default=MultiStrategy.ALL.value,
    show_default=True, help="Strategy to handle multiple perfect hits")
@click.option("--loci", help="Limit to these loci, provided as comma seperated string (e.g., 'abcZ,fumC')")
@_common_options
def call(
        fasta: Path,
        db: Path,
        out_json: Path,
        out_tsv: Path,
        out_dir: Path,
        export_novel: bool,
        keep_minimap2: bool,
        threads: int,
        min_id_novel: int,
        multi: str,
        loci: str | None,
        debug: bool,
        log: Path) -> None:
    """
    Calls alleles from a FASTA file.
    """
    # Setup
    initialize_logging(log_path=log, debug=debug)
    t0 = datetime.now()

    # Validate args
    if (export_novel or keep_minimap2) and (out_dir is None):
        raise click.UsageError("Output directory ('--out-dir') must be specified when exporting FASTA / Minimap2 files")

    # Create the output directory
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    # Update the JSON output path
    if (out_json is None) and (out_dir is not None):
        out_json_ = Path(out_dir, 'mist.json')
    elif out_json is not None:
        out_json_ = out_json
    else:
        out_json_ = Path('mist.json')

    # Call alleles
    caller = MistCaller(
        dir_db=db,
        multi=multi,
        keep_minimap2=keep_minimap2,
        loci=loci.split(',') if loci else None,
        export_novel=export_novel,
        min_id_novel=min_id_novel
    )
    caller.call_alleles(
        path_fasta=fasta,
        out_json=out_json_,
        out_tsv=out_tsv,
        out_dir=out_dir,
        threads=threads
    )
    logger.info("Make sure to cite the corresponding database when using this in your research")
    path_citation = files('mist').joinpath('resources/citation.txt')
    logger.info(f'Please cite: {path_citation.read_text()}')
    logger.info(f"Processing time: {(datetime.now() - t0).total_seconds():.2f} seconds")


@cli.command()
@click.option("--url", required=True, help="URL to download from.")
@click.option("-o", "--output", type=click.Path(path_type=Path), default=Path("mist_download"), show_default=True, help="Output directory")
@click.option("-p", "--include-profiles", is_flag=True, help="Download the profiles")
@click.option("-d", "--downloader", type=click.Choice(list(DOWNLOADERS.keys())), required=True, help="Downloader")
@click.option("--dir-tokens", type=click.Path(path_type=Path), default=Path.cwd() / ".bigsdb_tokens", show_default=True, help="Directory with access tokens")
@click.option("--key-name", default="PubMLST", show_default=True, help="Key name")
@click.option("--site", default="PubMLST", show_default=True, help="Site")
@_common_options
def download(
        url: str,
        output: Path,
        include_profiles: bool,
        downloader: str,
        dir_tokens: Path,
        key_name: str,
        site: str,
        debug: bool,
        log: Path) -> None:
    """
    Download (cg)MLST schemes from various sources.
    """
    initialize_logging(log_path=log, debug=debug)
    downloader = MistDownload(
        url=url,
        output=output,
        include_profiles=include_profiles,
        downloader=downloader,
        dir_tokens=dir_tokens,
        key_name=key_name,
        site=site
    )
    downloader.run()

@cli.command()
@click.argument('inputs', nargs=-1)
@click.option("-d", "--out-dists", type=click.Path(path_type=Path), default=Path("distances.tsv"), show_default=True, help="Distance matrix output")
@click.option("-m", "--out-matrix", type=click.Path(path_type=Path), default=Path("allele_matrix.tsv"), show_default=True, help="Allele matrix output")
@click.option("-l", "--min-perc-loci", type=int, default=90, show_default=True, help="Minimum percentage of loci that should be present in a dataset")
@click.option("-s", "--min-perc-samples", type=int, default=90, show_default=True, help="Minimum percentage of datasets where loci should be present")
@_common_options
def dists(
        inputs: tuple[Path],
        out_dists: Path,
        out_matrix: Path,
        min_perc_loci: int,
        min_perc_samples: int,
        debug: bool,
        log: Path) -> None:
    """
    Builds distance and allele matrices from MiST output files.
    """
    initialize_logging(log_path=log, debug=debug)
    if len(inputs) == 0:
        raise click.UsageError("No input files provided.")

    mist_dists = MistDists(
        inputs=[Path(x) for x in inputs],
        out_matrix=out_matrix,
        out_dists=out_dists,
        min_perc_loci=min_perc_loci,
        min_perc_samples=min_perc_samples
    )
    mist_dists.run()

@cli.command(name='list')
@click.option('-d', '--downloader', type=click.Choice(list(DOWNLOADERS.keys())), required=True, help="Downloader")
@click.option('-h', '--host', help="Host name (required only if multiple hosts are available for this downloader).")
@click.option('--db', type=str, help="BIGSdb only: database name")
def list_(downloader: str, host: str | None, db: str) -> None:
    """
    Lists the available schemes that can be downloaded.
    """
    # Validate the arguments
    if downloader in ('bigsdb', 'bigsdb_auth') and host is None:
        raise click.UsageError("--host option is required for this downloader. Valid options: 'pubmlst' or 'pasteur'.")

    # Retrieve available schemes
    mist_list = MistList(
        downloader=downloader,
        host=host,
        db=db
    )
    mist_list.print_available_schemes()


if __name__ == "__main__":
    cli()
