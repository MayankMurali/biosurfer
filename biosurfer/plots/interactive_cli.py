"""
CLI wiring for the interactive Plotly viewer. Reuses the same gene/
transcript_ids resolution logic as the existing static `plot` command
(analysis/plot_biosurfer.py's run_plot) rather than reinventing it, but as
its own new subcommand -- the existing `plot` command is completely
untouched.
"""
from pathlib import Path

import click

from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import Gene, Transcript
from biosurfer.plots.interactive import build_interactive_isoform_figure, save_interactive_plot


@click.command('plot_interactive')
@click.option('-v', '--verbose', is_flag=True, help='Print verbose messages')
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory in which to save the interactive plot')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', type=str, help='Name of gene for which to plot all isoforms; overrides TRANSCRIPT_IDS')
@click.argument('transcript_ids', nargs=-1)
def run_plot_interactive(verbose: bool, output: Path, db_name: str, gene: str, transcript_ids: tuple):
    """
    Plot isoforms from a single gene (or a list of TRANSCRIPT_IDS) as an
    interactive, hoverable/zoomable HTML figure -- a companion to the
    existing static `plot` command, not a replacement for it.
    """
    if not output:
        output = Path('.')
    db = Database(db_name)
    with db.get_session() as session:
        if gene:
            gene_obj = Gene.from_name(session, gene)
            if gene_obj is None:
                click.echo(f'Gene "{gene}" not found in database')
                return
            transcripts = list(gene_obj.transcripts)
            gene_name = gene
        else:
            tx_map = Transcript.from_accessions(session, transcript_ids)
            transcripts = [tx_map[tx_id] for tx_id in transcript_ids if tx_id in tx_map]
            for tx_id in transcript_ids:
                if tx_id not in tx_map:
                    click.echo(f'Transcript ID "{tx_id}" not found in database')
            if not transcripts:
                click.echo('No isoforms provided')
                return
            gene_name = transcripts[0].gene.name

        if verbose:
            click.echo(f'Building interactive plot for {len(transcripts)} isoform(s)...')
        fig = build_interactive_isoform_figure(transcripts)
        filepath = output / f'{db_name}-{gene_name}-interactive.html'
        save_interactive_plot(fig, filepath)
        click.echo(f'Saved {filepath}')
