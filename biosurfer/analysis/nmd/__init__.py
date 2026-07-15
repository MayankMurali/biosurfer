"""
Nonsense-mediated decay (NMD) analysis.

Additive Stage B module. Nothing here modifies ``biosurfer.core`` -- it only
reads from existing Transcript/ORF/Exon objects and surfaces/extends the
50-nt-rule ``ORF.nmd`` cached_property that already existed in
``biosurfer.core.models.biomolecules`` but was unused anywhere in the
codebase before this module was added.
"""
from biosurfer.analysis.nmd.predict import (
    VariantNMDEffect,
    compare_isoform_nmd,
    get_nmd_status,
    predict_variant_nmd_effect,
)

__all__ = [
    'VariantNMDEffect',
    'compare_isoform_nmd',
    'get_nmd_status',
    'predict_variant_nmd_effect',
]
