from sqlalchemy import select, and_
from biosurfer.core.models.biomolecules import Transcript, Protein, Gene
from biosurfer.core.models.genetics import GenomicVariant, GWASStatistic
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import SequenceAlignmentCategory, APPRIS
from biosurfer.core.models.nonpersistent import Position

def map_variant_to_isoforms(session, variant: GenomicVariant, anchor_tx: Transcript, other_tx: Transcript):
    """
    Determines the impact of a variant on a specific isoform pair comparison.
    """
    results = {
        "variant": variant.rsid,
        "pos": variant.position,
        "anchor_impact": None,
        "other_impact": None,
        "pblock_category": None,
        "status": "Variant affects both equally"
    }

    # 1. Check if variant is inside the genomic bounds of the transcripts
    # (Simplified check; assumes variant is on correct chrom)
    
    # Helper to check exon overlap
    def get_covering_exon(tx, pos):
        return next((e for e in tx.exons if e.start <= pos <= e.stop), None)

    anchor_exon = get_covering_exon(anchor_tx, variant.position)
    other_exon = get_covering_exon(other_tx, variant.position)

    # 2. Variant Spliced Out Logic
    if anchor_exon and not other_exon:
        results["status"] = "Variant Spliced Out in Alternative"
        results["pblock_category"] = "DELETION"
    elif not anchor_exon and other_exon:
        results["status"] = "Variant Unique to Alternative"
        results["pblock_category"] = "INSERTION"
    
    # 3. Ragged Codon / Missense Logic (Future Phase)
    # If variant is in both, we would check CodonAlignment here to see 
    # if it results in different Amino Acids due to frame/split codons.

    return results

def analyze_nterm_risk(session, gene_name="PPARG"):
    """
    Specifically analyzes N-terminal differences between the Reference Isoform
    and all other isoforms to find GWAS hits located in unique N-terminal regions.
    """
    print(f"--- Analyzing N-terminal Risk for {gene_name} ---")

    # 1. Load Gene
    gene = Gene.from_name(session, gene_name)
    if not gene:
        print(f"Gene {gene_name} not found in database.")
        return

    # 2. Identify Anchor (Heuristic: APPRIS Principal or Longest Protein)
    transcripts = gene.transcripts
    if not transcripts:
        print("No transcripts found.")
        return

    # Sort by APPRIS (Principal first) then Protein Length
    anchor = sorted(
        [t for t in transcripts if t.protein],
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True
    )[0]
    
    print(f"Reference Isoform: {anchor.name} ({anchor.accession})")

    # 3. Compare against all other isoforms
    hits_found = False
    
    for other in transcripts:
        if other == anchor or not other.protein:
            continue

        try:
            # Run Hybrid Alignment
            aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
            
            # Check the very first block of the alignment
            first_block = aln.blocks[0]
            
            target_transcript = None
            protein_range = None
            desc = None

            # Case A: Deletion at start (Unique to Reference)
            if first_block.category == SequenceAlignmentCategory.DELETION:
                if first_block.anchor_range.start == 0:
                    target_transcript = anchor
                    protein_range = first_block.anchor_range
                    desc = "Unique N-term in Reference"

            # Case B: Insertion at start (Unique to Alternative)
            # This matches PPARG2 (Alternative) vs PPARG1 (Reference/Anchor) scenario
            elif first_block.category == SequenceAlignmentCategory.INSERTION:
                if first_block.other_range.start == 0:
                    target_transcript = other
                    protein_range = first_block.other_range
                    desc = "Unique N-term in Alternative"

            if target_transcript and protein_range:
                print(f"\nComparison: Ref vs {other.name}")
                print(f"  Pattern: {desc} (Length: {len(protein_range)} AA)")
                
                # 4. Map Protein Range to Genomic Coordinates to find GWAS hits
                gwas_hits = _query_gwas_in_protein_range(session, target_transcript, protein_range)
                
                if gwas_hits:
                    hits_found = True
                    print(f"  [!] RISK DETECTED: {len(gwas_hits)} GWAS hits in this unique region!")
                    for hit in gwas_hits:
                        print(f"    -> rsID: {hit.variant.rsid or 'N/A'} | Pos: {hit.variant.position} | P-val: {hit.p_value}")
                else:
                    print(f"  No GWAS hits found in this N-terminal region.")

        except Exception as e:
            print(f"  Error aligning {other.name}: {e}")

    if not hits_found:
        print("\nNo N-terminal specific GWAS hits found for this gene.")

def _query_gwas_in_protein_range(session, transcript, protein_range):
    """
    Helper: Converts a protein range (AA indices) into genomic coordinates
    and queries the GWASStatistic table.
    """
    hits = []
    
    # 1. Convert Protein Range -> Transcript Nucleotide Range
    # AA index 0 -> nucleotides 0, 1, 2 in CDS (relative to ORF start)
    start_aa = protein_range.start
    end_aa = protein_range.stop - 1 # inclusive index
    
    # Get nucleotide indices relative to the transcript start (0-based)
    tx_start_nt = transcript.protein.get_transcript_coord_from_protein_coord(start_aa)
    tx_end_nt = transcript.protein.get_transcript_coord_from_protein_coord(end_aa) + 2 # include full codon
    
    # 2. Iterate over Exons to handle splicing
    # The protein region might span multiple exons (and skip introns).
    # We must intersect the [tx_start_nt, tx_end_nt] range with every exon.
    
    for exon in transcript.exons:
        # Exon bounds in transcript coordinates (0-based, derived during DB load)
        # Note: In Biosurfer DB, exon.transcript_start/stop are 1-based.
        ex_tx_start = exon.transcript_start - 1
        ex_tx_stop = exon.transcript_stop - 1
        
        # Determine overlap
        overlap_start = max(tx_start_nt, ex_tx_start)
        overlap_end = min(tx_end_nt, ex_tx_stop)
        
        if overlap_start <= overlap_end:
            # We have an overlap in this exon. Convert to Genomic Coords.
            # Use Biosurfer's coordinate mapper
            try:
                g_start = transcript.get_genome_coord_from_transcript_coord(overlap_start).coordinate
                g_end = transcript.get_genome_coord_from_transcript_coord(overlap_end).coordinate
                
                # Handle minus strand (g_start might be > g_end)
                if g_start > g_end:
                    g_start, g_end = g_end, g_start
                
                # 3. Query Database
                stmt = select(GWASStatistic).join(GenomicVariant).where(
                    and_(
                        GenomicVariant.chromosome == transcript.gene.chromosome_id,
                        GenomicVariant.position >= g_start,
                        GenomicVariant.position <= g_end
                    )
                )
                stats = session.execute(stmt).scalars().all()
                hits.extend(stats)
                
            except Exception as e:
                # Coordinate mapping might fail at edge cases (e.g. polyA)
                pass
                
    return hits
