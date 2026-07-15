"""
Splice-site-strength scoring (SpliceAI) at isoform-specific junctions.

Additive Stage B module. Nothing here modifies ``biosurfer.core``.

Answers: "is this alternative splice site predicted stronger/weaker than
the reference?" -- scores the donor/acceptor probability SpliceAI assigns
at each of an isoform's exon-exon junctions, and compares junctions unique
to one isoform against the anchor.

``spliceai`` + ``tensorflow`` are *optional* dependencies
(``pip install biosurfer[splicing]``); the base install stays free of them.
Genomic reference sequence is read via ``pysam.FastaFile`` -- ``pysam`` is
already a hard biosurfer dependency (used for VCF loading), so no new
dependency is needed for that half.

Licensing note: SpliceAI's trained models are distributed under
CC BY-NC 4.0 (academic/non-commercial use only) and its source code under
the PolyForm Strict License -- fine for this project's academic use, but
flagged here since it's a real constraint on how this module's output may
be used or redistributed.
"""
from biosurfer.analysis.splicing.predict import compare_isoform_splice_sites, score_junction

__all__ = ['compare_isoform_splice_sites', 'score_junction']
