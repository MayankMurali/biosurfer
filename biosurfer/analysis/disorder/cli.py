"""CLI wiring for the disorder module, registered via cli.add_command like the others."""
from pathlib import Path

import click

from biosurfer.core.database import Database
from biosurfer.analysis.disorder.report import run_disorder_report


@click.command('analyze_disorder')
@click.option('-v', '--verbose', is_flag=True, help='Print verbose messages')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', required=True, help='Target gene to analyze (e.g. PPARG)')
@click.option('-o', '--output', type=click.Path(file_okay=False, writable=True, path_type=Path), required=True, help='Directory for output tables')
def run_analyze_disorder(verbose: bool, db_name: str, gene: str, output: Path):
    """
    Scores intrinsic disorder (metapredict V3) for every alternative
    isoform vs. the gene's anchor isoform, over each protein/codon
    alignment block, to see whether alternative splicing disrupts
    disordered regions. Complementary to (not a replacement for) the
    existing MobiDB-lite IDR annotations already loaded via
    `load_feature_mappings`.

    Requires the optional 'metapredict' dependency (pip install biosurfer[disorder]).
    """
    if verbose:
        click.echo(f'Initializing database: {db_name}...')
    db = Database(db_name)
    session = db.get_session()
    run_disorder_report(session, gene, output)
