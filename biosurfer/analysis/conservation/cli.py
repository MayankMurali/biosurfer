"""
CLI wiring for the conservation/AlphaMissense module. Registered via
``cli.add_command`` in ``biosurfer.py``, same pattern as ``illustrate`` and
``analyze_nmd``.
"""
from pathlib import Path

import click

from biosurfer.core.database import Database
from biosurfer.analysis.conservation.predict import DEFAULT_PHYLOP_HG38_URL, open_bigwig
from biosurfer.analysis.conservation.report import run_conservation_report


@click.command('analyze_conservation')
@click.option('-v', '--verbose', is_flag=True, help='Print verbose messages')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', required=True, help='Target gene to analyze (e.g. PPARG)')
@click.option(
    '--bigwig', is_flag=False, flag_value=DEFAULT_PHYLOP_HG38_URL,
    help=(
        'Path or URL to a phyloP/phastCons bigWig conservation track. '
        f'Pass with no value to use the UCSC hg38 100-way phyloP track ({DEFAULT_PHYLOP_HG38_URL}); '
        'omit entirely to skip conservation scoring. '
        'Requires the optional pyBigWig dependency (pip install biosurfer[conservation]).'
    ),
)
@click.option('--no-conservation', is_flag=True, help='Skip conservation scoring even if --bigwig is not given (default behavior).')
@click.option(
    '--alphamissense', default=None, type=click.Path(exists=True, dir_okay=False),
    help=(
        'Path to a bgzipped + tabix-indexed AlphaMissense TSV '
        '(AlphaMissense_hg38.tsv.gz from github.com/google-deepmind/alphamissense, '
        "indexed with `tabix -s1 -b2 -e2 -c'#' AlphaMissense_hg38.tsv.gz`). Omit to skip."
    ),
)
@click.option('-o', '--output', type=click.Path(file_okay=False, writable=True, path_type=Path), required=True, help='Directory for output tables')
def run_analyze_conservation(verbose: bool, db_name: str, gene: str, bigwig: str, no_conservation: bool, alphamissense: str, output: Path):
    """
    Scores a gene's alternative-vs-anchor isoform differences (pblocks/cblocks)
    for nucleotide conservation (phyloP/phastCons) and/or AlphaMissense
    missense pathogenicity, over the genomic span of each block.
    """
    use_conservation = bool(bigwig) and not no_conservation
    if not use_conservation and not alphamissense:
        raise click.UsageError('Provide --bigwig and/or --alphamissense -- at least one data source is required.')

    if verbose:
        click.echo(f'Initializing database: {db_name}...')
    db = Database(db_name)
    session = db.get_session()

    bw = None
    try:
        if use_conservation:
            bw_path = bigwig if bigwig else DEFAULT_PHYLOP_HG38_URL
            if verbose:
                click.echo(f'Opening conservation track: {bw_path}...')
            bw = open_bigwig(bw_path)
        run_conservation_report(session, gene, output, bw=bw, alphamissense_tabix=alphamissense)
    finally:
        if bw is not None:
            bw.close()
