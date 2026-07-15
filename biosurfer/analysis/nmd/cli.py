"""
CLI wiring for the NMD module, kept in its own file and registered via
``cli.add_command`` -- the same extension pattern ``illustrate`` already
uses in ``biosurfer.py`` -- so that ``biosurfer.py`` itself only needs a
two-line addition (import + ``cli.add_command``) to pick this up.
"""
from pathlib import Path

import click

from biosurfer.core.database import Database
from biosurfer.analysis.nmd.report import (
    run_isoform_nmd_delta_report,
    run_transcript_nmd_report,
    run_variant_nmd_report,
)


@click.command('analyze_nmd')
@click.option('-v', '--verbose', is_flag=True, help='Print verbose messages')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', required=True, help='Target gene to analyze (e.g. PPARG)')
@click.option(
    '--mode', required=True,
    type=click.Choice(['transcript', 'isoform', 'variant'], case_sensitive=False),
    help=(
        "'transcript': per-transcript NMD status for every transcript of the gene. "
        "'isoform': NMD status of every alternative isoform vs. the gene's anchor isoform. "
        "'variant': for each SNP in --vcf overlapping the gene, whether it flips a transcript's NMD status."
    ),
)
@click.option('--vcf', type=click.Path(exists=True, dir_okay=False), help="Path to VCF file (bgzipped & tabix indexed). Required for --mode=variant.")
@click.option('-o', '--output', type=click.Path(file_okay=False, writable=True, path_type=Path), required=True, help='Directory for output tables')
def run_analyze_nmd(verbose: bool, db_name: str, gene: str, mode: str, vcf: str, output: Path):
    """
    Predicts nonsense-mediated decay (NMD) susceptibility for a gene's
    transcripts, using the standard 50-nt rule (already implemented as
    ORF.nmd in biosurfer.core.models.biomolecules).
    """
    if mode == 'variant' and not vcf:
        raise click.UsageError('--vcf is required when --mode=variant')

    if verbose:
        click.echo(f'Initializing database: {db_name}...')
    db = Database(db_name)
    session = db.get_session()

    if mode == 'transcript':
        run_transcript_nmd_report(session, gene, output)
    elif mode == 'isoform':
        run_isoform_nmd_delta_report(session, gene, output)
    elif mode == 'variant':
        run_variant_nmd_report(session, gene, vcf, output)
