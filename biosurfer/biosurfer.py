import csv
from pathlib import Path

import click

from biosurfer.analysis.genome_wide_alignment_analysis import run_hybrid_alignment_for_all_genes
from biosurfer.core.database import Database
from biosurfer.core.helpers import (get_ids_from_gencode_fasta,
                                    get_ids_from_lrp_fasta,
                                    get_ids_from_pacbio_fasta, skip_gencode,
                                    skip_par_y)
from biosurfer.analysis.plot_biosurfer import run_plot

from biosurfer.analysis.genetics_analyzer import analyze_nterm_risk

@click.group(chain=True)
def cli():
    """
    \b
--------------------------------------
         Welcome to Biosurfer!
--------------------------------------
    """

@cli.command("load_db")
@click.option('-v', '--verbose', is_flag=True, help="Will print verbose messages")
@click.option('-d', '--db_name', required=True, help="Database name")
@click.option('--source', type=click.Choice(['GENCODE', 'PacBio'], case_sensitive=False),  required=True, help="Source of input data")
@click.option('--gtf', required=True, type=click.Path(exists=True), help='Path to gtf file')
@click.option('--tx_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to transcript sequence fasta file')
@click.option('--tl_fasta', required=True, type=click.Path(exists=True, path_type=Path), help='Path to protein sequence fasta file')
@click.option('--sqanti', type=click.Path(exists=True, path_type=Path), help='Path to SQANTI classification tsv file (only for PacBio isoforms)')
def run_populate_database(verbose: bool, db_name: str, source: str, gtf: Path, tx_fasta: Path, tl_fasta: Path, sqanti: Path):
    """Loads transcript and protein isoform information from provided files into a Biosurfer database.
    A new database is created if the target database does not exist."""

    db = Database(db_name)
    if source == "GENCODE":
        click.echo('----- Loading database with reference ', err=True)
        db.load_gencode_gtf(gtf)
        db.load_transcript_fasta(tx_fasta, get_ids_from_gencode_fasta, skip_par_y)
        db.load_translation_fasta(tl_fasta, get_ids_from_gencode_fasta, skip_par_y)

    elif source == "PacBio":
        click.echo('----- Loading database without reference ', err=True)
        db.load_pacbio_gtf(gtf)
        db.load_transcript_fasta(tx_fasta, get_ids_from_pacbio_fasta)
        db.load_translation_fasta(tl_fasta, get_ids_from_lrp_fasta, skip_gencode)
        if sqanti:
            db.load_sqanti_classifications(sqanti)

@cli.command("hybrid_alignment")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory for output files')
@click.option('--gencode', is_flag=True, help='Also compare all GENCODE isoforms of a gene against its anchor isoform')
@click.option('--anchors', type=click.Path(exists=True, dir_okay=False, path_type=Path), help='TSV file with gene names in column 1 and anchor isoform IDs in column 2')
def run_hybrid_al(verbose, db_name, output: Path, gencode: bool, anchors: Path):
    """ This script runs hybrid alignment on the provided database. """
    click.echo('')
    click.echo('----- Running hybrid alignment: ', err=True)
    click.echo('')
    if not output:
        output = Path('.')
    if anchors:
        with open(anchors) as f:
            gene_to_anchor_tx = {gene: tx for gene, tx in csv.reader(f, delimiter='\t')}
    else:
        gene_to_anchor_tx = None
    run_hybrid_alignment_for_all_genes(db_name, output, gencode, gene_to_anchor_tx)

@cli.command("plot")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-o', '--output', type=click.Path(exists=True, file_okay=False, writable=True, path_type=Path), help='Directory in which to save plots')
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', type=str, help='Name of gene for which to plot all isoforms; overrides TRANSCRIPT_IDS')
@click.option('--snps', is_flag=True, help='Overlay genetic variants (SNPs) from the database onto the plot')  # <--- NEW OPTION
@click.argument('transcript_ids', nargs=-1)
def plot_isoforms(verbose: bool, output: Path, gene: str, db_name: str, snps: bool, transcript_ids: tuple[str]):
    """Plot isoforms from a single gene, specified by TRANSCRIPT_IDS."""
    # Pass the snps flag to the run_plot function
    run_plot(output, gene, db_name, transcript_ids, show_snps=snps)


@cli.command("analyze_nterm")
@click.option('-v', '--verbose', is_flag=True, help="Print verbose messages")
@click.option('-d', '--db_name', required=True, nargs=1, help='Database name')
@click.option('--gene', required=True, help='Target gene to analyze (e.g., PPARG)')
@click.option('--vcf', required=True, help="Path to VCF file (bgzipped & tabix indexed)")
@click.option('--gwas', required=True, help="Path to GWAS summary statistics (TSV)")
@click.option('--trait', default='GWAS', help="Trait label to record for the loaded GWAS summary statistics (e.g. T2D)")
@click.option('-o', '--output', type=click.Path(file_okay=False, writable=True, path_type=Path), help='Directory for output tables')
def analyze_nterm_risk_cmd(verbose, db_name, gene, vcf, gwas, trait, output):
    """
    Analyzes N-terminal differences for a specific gene to identify
    GWAS hits located in unique N-terminal regions (e.g. PPARG1 vs PPARG2).
    """
    if verbose:
        click.echo(f"Initializing database: {db_name}...")

    db = Database(db_name)

    # 1. Load Genetics Data into DB
    if verbose:
        click.echo(f"Loading genetics data for {gene}...")

    db.load_genetics_data(vcf_path=vcf, gwas_path=gwas, gene_name=gene, trait=trait)
    
    # 2. Run Analysis
    if verbose:
        click.echo(f"Running N-terminal risk analysis...")

    # Pass the output directory to the analyzer
    analyze_nterm_risk(db.get_session(), gene_name=gene, output_dir=output)