# %%
from pathlib import Path
from more_itertools import partition
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import Gene, Transcript
from biosurfer.plots.plotting import IsoformPlot

#%%
def run_plot(output: Path, gene: str, db_name: str, transcript_ids: tuple[str], show_snps: bool = False):
    """ Main plot function to invoke plotting for different pipelines/scripts.
    Args:
      output: Path to save the plot
      gene: Name of the gene to plot
      db_name: Name of the database to connect to
      transcript_ids: Specific transcript IDs to plot (if gene not specified)
      show_snps: Boolean flag to visualize genetic variants
    Returns:
      Nothing
    """
    if not output:
        output = Path('.')
    db = Database(db_name)
    with db.get_session() as s:
        print(f'Loading transcripts from database...')
        
        if gene:
            gene_obj = Gene.from_name(s, gene)
            if gene_obj is None:
                print(f'Gene "{gene}" not found in database')
                transcripts = dict()
                anchor = None
            else:
                transcripts = {tx.accession: tx for tx in gene_obj.transcripts}
                anchor = max(transcripts.values(), key=lambda tx: getattr(tx, 'appris', APPRIS.NONE))
            others = [tx for tx in transcripts.values() if tx is not anchor]
        else:
            transcripts: dict[str, Transcript] = {tx.accession: tx for tx in Transcript.from_accessions(s, transcript_ids).values()}
            not_found, found = partition(lambda tx_id: tx_id in transcripts, transcript_ids)
            for tx_id in not_found:
                print(f'Transcript ID "{tx_id}" not found in database')
            if transcript_ids:
                anchor = transcripts.get(transcript_ids[0], None)
            else:
                print('No isoforms provided')
                anchor = None
            others = [tx for tx in map(transcripts.get, found) if tx is not anchor]

        if anchor:
            print(f'Reference isoform: {anchor}')
            gene_name = anchor.gene.name

            alns: dict[Transcript, ProteinAlignment] = dict()
            for other in others:
                if anchor.protein is None or other.protein is None:
                    alns[other] = None
                else:
                    try:
                        alns[other] = ProteinAlignment.from_proteins(anchor.protein, other.protein)
                    except ValueError:
                        print(f'Could not plot isoform {other}')
            
            filename = f'{db_name}-{gene_name}.png'
            plot = IsoformPlot([anchor] + list(alns.keys()))
            plot.draw_all_isoforms()
            plot.draw_frameshifts()

            # --- SNP Visualization Block ---
            if show_snps:
                print(f"Retrieving variants for {gene_name}...")
                # Note: This requires the get_gene_variants method to be implemented in Database class
                variants = db.get_gene_variants(s, gene_name)
                if variants:
                    print(f"Found {len(variants)} variants. Adding to plot.")
                    plot.draw_variants(variants)
                else:
                    print("No variants found in database for this gene.")
            # -------------------------------
            
            for other, aln in alns.items():
                if aln:
                    plot.draw_protein_alignment_blocks(aln.blocks, anchor.protein, other.protein)
            
            plot.draw_legend()
            filepath = str(output/filename)
            plot.savefig(filepath)
            print(f'Saved {filepath}')


if __name__ == "__main__":
    pass