"""
Intrinsically disordered region (IDR) disruption scoring via metapredict.

Additive Stage B module. Nothing here modifies ``biosurfer.core``.

Note on scope: Biosurfer's existing ``core.data_loaders.load_feature_mappings``
already loads MobiDB-lite IDR annotations as ``ProteinFeature`` rows (same
pipeline as Pfam domains), and ``core.alignments`` already projects them
onto alternative isoforms via ``ProjectedFeature``/``altered_residues`` --
this is live, working infrastructure, not dormant like ``ORF.nmd`` was.
This module is a deliberately independent, complementary data source
(metapredict's ML-based per-residue disorder score) rather than a
duplicate of that machinery -- it doesn't touch the database at all,
comparing anchor-vs-alternative disorder scores directly over protein
alignment blocks, the same pattern used in
``biosurfer.analysis.conservation``.

``metapredict`` is an *optional* dependency (``pip install
biosurfer[disorder]``); the base install stays free of it.
"""
from biosurfer.analysis.disorder.predict import compare_isoform_disorder, score_disorder

__all__ = ['compare_isoform_disorder', 'score_disorder']
