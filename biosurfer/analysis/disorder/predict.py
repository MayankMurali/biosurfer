"""
Core disorder-scoring logic. Scores the *whole* protein sequence once (so
metapredict's network sees real flanking context everywhere, including
near block boundaries) and slices the cached per-residue array for
whatever protein range is being compared -- much cheaper than rescoring an
isolated fragment per block, and more accurate near boundaries.
"""
from typing import List, Optional

from biosurfer.core.models.biomolecules import Protein, Transcript

try:
    import metapredict as meta
    _METAPREDICT_AVAILABLE = True
except ImportError:
    _METAPREDICT_AVAILABLE = False


def _require_metapredict():
    if not _METAPREDICT_AVAILABLE:
        raise ImportError(
            "Disorder scoring requires the optional 'metapredict' dependency, which "
            "isn't part of the base biosurfer install. Install it with:\n"
            "    pip install biosurfer[disorder]"
        )


def score_disorder(sequence: str) -> List[float]:
    """Per-residue disorder score (0-1) for a full protein sequence, via metapredict V3."""
    _require_metapredict()
    scores = meta.predict_disorder(sequence)
    return [float(s) for s in scores]


def score_protein_range_disorder(disorder_scores: List[float], protein_range) -> dict:
    """
    Mean/max/min disorder score over a protein coordinate range (0-based
    half-open, e.g. a pblock/cblock's ``anchor_range``/``other_range``),
    sliced out of an already-computed full-sequence score array.
    """
    segment = disorder_scores[protein_range.start:protein_range.stop]
    if not segment:
        return {'mean': None, 'max': None, 'min': None}
    return {'mean': sum(segment) / len(segment), 'max': max(segment), 'min': min(segment)}


def compare_isoform_disorder(anchor: Transcript, other: Transcript) -> List[dict]:
    """
    For every protein/codon alignment block between ``other`` and the
    anchor isoform, compare mean disorder score in the anchor's sequence
    vs. the alternative's sequence over that block's range. Mirrors
    ``biosurfer.analysis.conservation.report``'s alignment-block iteration.
    """
    from biosurfer.core.alignments import CodonAlignment, ProteinAlignment, TranscriptAlignment
    from biosurfer.core.splice_events import get_event_code

    anchor_scores = score_disorder(anchor.protein.sequence)
    other_scores = score_disorder(other.protein.sequence)

    tx_aln = TranscriptAlignment.from_transcripts(anchor, other)
    cd_aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    pr_aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)

    rows = []
    for p, pblock in enumerate(pr_aln.blocks):
        cblocks = pr_aln.pblock_to_cblocks.get(pblock, [])
        for c, cblock in enumerate(cblocks):
            tblock = cd_aln.cblock_to_tblock.get(cblock)
            events = tx_aln.block_to_events.get(tblock, ())
            anchor_disorder = score_protein_range_disorder(anchor_scores, cblock.anchor_range) if cblock.anchor_range else None
            other_disorder = score_protein_range_disorder(other_scores, cblock.other_range) if cblock.other_range else None
            rows.append({
                'pblock_number': p,
                'pblock_category': pblock.category.name,
                'cblock_number': c,
                'cblock_category': cblock.category.name,
                'events': get_event_code(events),
                'anchor_disorder_mean': anchor_disorder['mean'] if anchor_disorder else None,
                'other_disorder_mean': other_disorder['mean'] if other_disorder else None,
            })
    return rows
