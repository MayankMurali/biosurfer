from biosurfer.core.models.biomolecules import Transcript, Protein
from biosurfer.core.models.genetics import GenomicVariant
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import SequenceAlignmentCategory

def map_variant_to_isoforms(session, variant: GenomicVariant, anchor_tx: Transcript, other_tx: Transcript):
    """
    Determines the impact of a variant on a specific isoform pair comparison.
    """
    results = {
        "variant": variant.rsid,
        "pos": variant.position,
        "anchor_impact": None,
        "other_impact": None,
        "pblock_category": None
    }

    # 1. Check if variant is inside the transcript bounds
    # Note: Biosurfer transcripts have genomic start/stop properties
    # We need to check exons specifically for finer granularity
    
    anchor_exon = next((e for e in anchor.exons if e.start <= variant.position <= e.stop), None)
    other_exon = next((e for e in other_tx.exons if e.start <= variant.position <= e.stop), None)

    # 2. If it falls in a deletion p-block (present in anchor, missing in other)
    if anchor_exon and not other_exon:
        results["status"] = "Variant Spliced Out in Alternative"
        # Determine specific p-block
        # (Requires instantiating ProteinAlignment)
        aln = ProteinAlignment.from_proteins(anchor_tx.protein, other_tx.protein)
        
        # Convert genomic coord to transcript coord, then protein coord
        try:
            tx_coord = anchor_tx.get_transcript_coord_from_genome_coord(
                Position(variant.chromosome, anchor_tx.strand, variant.position)
            )
            # Check if this falls into a Deletion p-block
            # (Logic required to traverse alignment blocks to find specific residue)
        except:
            pass

    return results
