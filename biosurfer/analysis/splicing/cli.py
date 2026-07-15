"""CLI wiring for the splicing module, registered via cli.add_command like the others."""
from pathlib import Path

import click
import pysam

from biosurfer.core.database import Database
from biosurfer.analysis.splicing.predict import DEFAULT_CONTEXT, load_models
from biosurfer.analysis.splicing.report import run_splicing_report


@click.command('analyze_splicing')
@click.option('-v', '--verbose', is_flag=True, help='Print verbose messages')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', required=True, help='Target gene to analyze (e.g. PPARG)')
@click.option('--genome-fasta', required=True, type=click.Path(exists=True, dir_okay=False), help='Indexed reference genome FASTA (.fai companion required)')
@click.option('--context', default=DEFAULT_CONTEXT, show_default=True, help='SpliceAI receptive-field context size (bp)')
@click.option('-o', '--output', type=click.Path(file_okay=False, writable=True, path_type=Path), required=True, help='Directory for output tables')
def run_analyze_splicing(verbose: bool, db_name: str, gene: str, genome_fasta: str, context: int, output: Path):
    """
    Scores SpliceAI donor/acceptor probability at every splice junction
    unique to an alternative isoform or unique to the gene's anchor
    isoform, to compare predicted splice-site strength between them.

    Requires the optional 'spliceai'/'tensorflow' dependencies
    (pip install biosurfer[splicing]) and a samtools-faidx-indexed
    reference genome FASTA matching the database's assembly.
    """
    if verbose:
        click.echo(f'Initializing database: {db_name}...')
    db = Database(db_name)
    session = db.get_session()

    if verbose:
        click.echo('Loading SpliceAI models (first call is slow -- loads 5 Keras models)...')
    models = load_models()

    fasta = pysam.FastaFile(genome_fasta)
    try:
        run_splicing_report(session, gene, output, models, fasta)
    finally:
        fasta.close()
