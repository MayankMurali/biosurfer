"""
Conservation (phyloP/phastCons) and AlphaMissense pathogenicity scoring.

Additive Stage B module. Nothing here modifies ``biosurfer.core``.

Two independent, optional data sources:

1. Nucleotide-level conservation via bigWig tracks (``pyBigWig`` -- an
   *optional* dependency, ``pip install biosurfer[conservation]``, since the
   base install must stay free of it).
2. AlphaMissense's precomputed missense pathogenicity scores, read via
   ``pysam`` -- already a hard dependency (used for VCF loading), so this
   needs no new dependency at all.
"""
from biosurfer.analysis.conservation.predict import (
    query_alphamissense,
    score_genomic_range,
    score_protein_range_conservation,
)

__all__ = [
    'query_alphamissense',
    'score_genomic_range',
    'score_protein_range_conservation',
]
