import pandas as pd
from pathlib import Path
from itertools import chain
from more_itertools import one
from sqlalchemy import select, and_

from biosurfer.core.models.biomolecules import Transcript, Protein, Gene
from biosurfer.core.models.genetics import GenomicVariant, GWASStatistic
from biosurfer.core.alignments import TranscriptAlignment, CodonAlignment, ProteinAlignment
from biosurfer.core.constants import SequenceAlignmentCategory, APPRIS
from biosurfer.core.splice_events import get_event_code

# FIX: Alias the long name to SeqAlignCat locally
SeqAlignCat = SequenceAlignmentCategory

def _query_gwas_in_protein_range(session, transcript, protein_range):
    """
    Helper: Converts a protein range (AA indices) into genomic coordinates
    and queries the GWASStatistic table.
    """
    hits = []
    
    # AA index 0 -> nucleotides 0, 1, 2 in CDS (relative to ORF start)
    start_aa = protein_range.start
    end_aa = protein_range.stop - 1 # inclusive index
    
    # Get nucleotide indices relative to the transcript start (0-based)
    # Note: These methods might return None if the range is outside CDS, so we handle safely
    try:
        tx_start_nt = transcript.protein.get_transcript_coord_from_protein_coord(start_aa)
        tx_end_nt = transcript.protein.get_transcript_coord_from_protein_coord(end_aa) + 2 
    except (ValueError, TypeError):
        return []
    
    # Iterate over Exons to handle splicing
    for exon in transcript.exons:
        # Exon bounds in transcript coordinates (0-based, derived during DB load)
        ex_tx_start = exon.transcript_start - 1
        ex_tx_stop = exon.transcript_stop - 1
        
        # Determine overlap
        overlap_start = max(tx_start_nt, ex_tx_start)
        overlap_end = min(tx_end_nt, ex_tx_stop)
        
        if overlap_start <= overlap_end:
            try:
                g_start = transcript.get_genome_coord_from_transcript_coord(overlap_start).coordinate
                g_end = transcript.get_genome_coord_from_transcript_coord(overlap_end).coordinate
                
                if g_start > g_end:
                    g_start, g_end = g_end, g_start
                
                # Query Database
                stmt = select(GWASStatistic).join(GenomicVariant).where(
                    and_(
                        GenomicVariant.chromosome == transcript.gene.chromosome_id,
                        GenomicVariant.position >= g_start,
                        GenomicVariant.position <= g_end
                    )
                )
                stats = session.execute(stmt).scalars().all()
                hits.extend(stats)
            except Exception:
                pass
                
    return hits

def generate_genetics_report(session, gene: Gene, anchor: Transcript, others: list[Transcript], output_dir: Path):
    """
    Generates a full hybrid alignment report enriched with GWAS hits.
    Saves to output_dir/{gene}_cblocks_risk.tsv.
    """
    records = []
    print(f"Generating alignment report for {len(others)} isoforms against {anchor.name}...")
    
    principal_length = anchor.protein.length

    for alternative in others:
        if not alternative.protein:
            continue
            
        try:
            alternative_length = alternative.protein.length
            
            # Perform Alignments
            tx_aln = TranscriptAlignment.from_transcripts(anchor, alternative)
            cd_aln = CodonAlignment.from_proteins(anchor.protein, alternative.protein)
            pr_aln = ProteinAlignment.from_proteins(anchor.protein, alternative.protein)
            
            pblocks = pr_aln.blocks
            
            # Iterate through protein blocks
            for p, pblock in enumerate(pblocks):
                # Skip Match blocks if you only want differences, 
                # OR keep them to show where risks are (e.g. missense).
                # Here we process all blocks to find risks everywhere.
                
                # Get associated c-blocks
                cblocks = pr_aln.pblock_to_cblocks.get(pblock, [])
                
                for c, cblock in enumerate(cblocks):
                    tblock = cd_aln.cblock_to_tblock.get(cblock)
                    events = tx_aln.block_to_events.get(tblock, ())
                    
                    # 1. Base Alignment Data
                    row = {
                        'gene': gene.name,
                        'anchor': anchor.name,
                        'other': alternative.name,
                        'pblock_number': p,
                        'pblock_category': pblock.category.name,
                        'pblock_anchor_start': pblock.anchor_range.start,
                        'pblock_anchor_stop': pblock.anchor_range.stop,
                        'pblock_other_start': pblock.other_range.start,
                        'pblock_other_stop': pblock.other_range.stop,
                        'cblock_number': c,
                        'cblock_category': cblock.category.name,
                        'events': get_event_code(events),
                        'cblock_anchor_seq': anchor.protein.sequence[cblock.anchor_range.start:cblock.anchor_range.stop],
                        'cblock_other_seq': alternative.protein.sequence[cblock.other_range.start:cblock.other_range.stop],
                    }

                    # 2. Add Genetic Risk Info
                    hits = []
                    target_tx = None
                    target_range = None
                    
                    # If Match or Deletion, the sequence exists in Anchor
                    if cblock.anchor_range:
                        target_tx = anchor
                        target_range = cblock.anchor_range
                    # If Insertion, the sequence exists in Alternative
                    elif cblock.other_range:
                        target_tx = alternative
                        target_range = cblock.other_range
                        
                    if target_tx and target_range:
                        hits = _query_gwas_in_protein_range(session, target_tx, target_range)
                    
                    row['risk_variant_count'] = len(hits)
                    if hits:
                        # Format: rsID:Pos(P-val)
                        details = "; ".join([f"{h.variant.rsid or 'NA'}:{h.variant.position}(p={h.p_value:.2e})" for h in hits])
                        row['risk_variant_details'] = details
                    else:
                        row['risk_variant_details'] = ""

                    records.append(row)
        except Exception as e:
            # print(f"Skipping {alternative.name} due to alignment error: {e}")
            continue

    if records:
        output_file = output_dir / f"{gene.name}_cblocks_risk.tsv"
        pd.DataFrame(records).to_csv(output_file, sep='\t', index=False)
        print(f"Detailed risk report saved to: {output_file}")


def analyze_nterm_risk(session, gene_name, output_dir: Path = None):
    """
    Analyzes N-terminal differences for a specific gene to identify GWAS hits.
    If output_dir is provided, saves full hybrid alignment tables with risk info.
    """
    print(f"--- Analyzing N-terminal Risk for {gene_name} ---")

    # 1. Load Gene
    gene = Gene.from_name(session, gene_name)
    if not gene:
        print(f"Gene {gene_name} not found in database.")
        return

    # 2. Identify Anchor
    transcripts = gene.transcripts
    if not transcripts:
        print("No transcripts found.")
        return

    anchor = sorted(
        [t for t in transcripts if t.protein],
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True
    )[0]
    
    print(f"Reference Isoform: {anchor.name} ({anchor.accession})")

    # 3. Standard Console Report (N-term specific)
    hits_found = False
    others = [t for t in transcripts if t != anchor and t.protein]
    
    for other in others:
        try:
            aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
            if not aln.blocks: continue
            
            first_block = aln.blocks[0]
            
            target_transcript = None
            protein_range = None
            desc = None

            if first_block.category == SeqAlignCat.DELETION:
                if first_block.anchor_range.start == 0:
                    target_transcript = anchor
                    protein_range = first_block.anchor_range
                    desc = "Unique N-term in Reference"
            elif first_block.category == SeqAlignCat.INSERTION:
                if first_block.other_range.start == 0:
                    target_transcript = other
                    protein_range = first_block.other_range
                    desc = "Unique N-term in Alternative"

            if target_transcript and protein_range:
                print(f"\nComparison: Ref vs {other.name}")
                print(f"  Pattern: {desc} (Length: {len(protein_range)} AA)")
                gwas_hits = _query_gwas_in_protein_range(session, target_transcript, protein_range)
                
                if gwas_hits:
                    hits_found = True
                    print(f"  [!] RISK DETECTED: {len(gwas_hits)} GWAS hits in this unique region!")
                    for hit in gwas_hits:
                        print(f"    -> rsID: {hit.variant.rsid or 'N/A'} | Pos: {hit.variant.position} | P-val: {hit.p_value}")
                else:
                    print(f"  No GWAS hits found in this N-terminal region.")

        except Exception as e:
            pass

    if not hits_found:
        print("\nNo N-terminal specific GWAS hits found for this gene.")
        
    # 4. Full Report Generation (if requested)
    if output_dir:
        output_dir = Path(output_dir)
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        generate_genetics_report(session, gene, anchor, others, output_dir)