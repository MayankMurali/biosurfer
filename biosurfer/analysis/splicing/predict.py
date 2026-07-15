"""
Core SpliceAI scoring logic.

SpliceAI's "custom sequence" API (see its README) feeds a sequence of
length ``context + N`` through 5 ensembled Keras models and gets back
per-position acceptor/donor probabilities for the central ``N``
positions (the flanking ``context`` bases are consumed as receptive
field, not scored themselves). Here ``N=1`` -- every call scores exactly
one genomic position -- so the real genomic sequence in a
``context``-sized window around that position is fetched from an indexed
FASTA (real flanking sequence, not the README example's illustrative
N-padding, which would throw away the actual signal the model needs).

For a minus-strand position, the fetched (+strand) window is
reverse-complemented before scoring -- the position of interest, being
at the exact center of an odd-length window, stays at the same index
either way.
"""
from typing import List, Optional

import numpy as np
import pysam
from Bio.Seq import Seq

from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.models.nonpersistent import Junction

try:
    from keras.models import load_model
    from pkg_resources import resource_filename
    from spliceai.utils import one_hot_encode
    _SPLICEAI_AVAILABLE = True
except ImportError:
    _SPLICEAI_AVAILABLE = False

DEFAULT_CONTEXT = 10000  # SpliceAI's standard 5000bp-per-side receptive field

_models_cache = None  # lazy singleton -- loading 5 Keras models is expensive, do it once per process


def _require_spliceai():
    if not _SPLICEAI_AVAILABLE:
        raise ImportError(
            "Splice-site scoring requires the optional 'spliceai' and 'tensorflow' "
            "dependencies, which aren't part of the base biosurfer install. Install with:\n"
            "    pip install biosurfer[splicing]\n"
            "Note: SpliceAI's trained models are distributed under CC BY-NC 4.0 "
            "(academic/non-commercial use) and its source under the PolyForm Strict "
            "License -- see github.com/Illumina/SpliceAI for details."
        )


def load_models() -> list:
    """Load (and cache) SpliceAI's 5-model ensemble."""
    global _models_cache
    _require_spliceai()
    if _models_cache is None:
        paths = (resource_filename('spliceai', f'models/spliceai{x}.h5') for x in range(1, 6))
        _models_cache = [load_model(path) for path in paths]
    return _models_cache


def _fetch_window(fasta: 'pysam.FastaFile', chromosome: str, strand: Strand, position: int, context: int = DEFAULT_CONTEXT) -> str:
    """
    Fetch a ``context + 1``-bp real genomic window centered on ``position``
    (1-based), oriented along ``strand``. Pads with 'N' only for the part
    of the window that would run past the chromosome's actual ends.
    """
    half = context // 2
    fetch_start = position - half - 1  # 0-based, pysam.fetch convention
    fetch_end = position + half        # 0-based, exclusive
    chrom_len = fasta.get_reference_length(chromosome)
    clipped_start = max(0, fetch_start)
    clipped_end = min(chrom_len, fetch_end)
    seq = fasta.fetch(chromosome, clipped_start, clipped_end).upper()
    seq = ('N' * (clipped_start - fetch_start)) + seq + ('N' * (fetch_end - clipped_end))
    if strand is Strand.MINUS:
        seq = str(Seq(seq).reverse_complement())
    return seq


def score_position(models: list, fasta: 'pysam.FastaFile', chromosome: str, strand: Strand, position: int, context: int = DEFAULT_CONTEXT) -> dict:
    """
    Return SpliceAI's acceptor/donor probability at a single 1-based
    genomic ``position``, oriented along ``strand``.
    """
    _require_spliceai()
    seq = _fetch_window(fasta, chromosome, strand, position, context)
    x = one_hot_encode(seq)[None, :]
    y = np.mean([model.predict(x) for model in models], axis=0)
    # output length == 1 (the single scored position); index 0 is it
    return {'acceptor_prob': float(y[0, 0, 1]), 'donor_prob': float(y[0, 0, 2])}


def score_junction(models: list, fasta: 'pysam.FastaFile', junction: Junction, context: int = DEFAULT_CONTEXT) -> dict:
    """
    Score a splice junction's donor site (at ``junction.donor``) and
    acceptor site (at ``junction.acceptor``) -- the biologically relevant
    read-out for each side.
    """
    donor, acceptor = junction.donor, junction.acceptor
    donor_scores = score_position(models, fasta, donor.chromosome, donor.strand, donor.coordinate, context)
    acceptor_scores = score_position(models, fasta, acceptor.chromosome, acceptor.strand, acceptor.coordinate, context)
    return {'donor_prob': donor_scores['donor_prob'], 'acceptor_prob': acceptor_scores['acceptor_prob']}


def compare_isoform_splice_sites(models: list, fasta: 'pysam.FastaFile', anchor: Transcript, other: Transcript, context: int = DEFAULT_CONTEXT) -> List[dict]:
    """
    Score every junction unique to ``anchor`` or unique to ``other``
    (junctions both isoforms share are skipped -- nothing splicing-wise
    differs there). Answers: is this alternative splice site predicted
    stronger/weaker than the sites the anchor isoform actually uses?
    """
    anchor_junctions = set(anchor.junctions)
    other_junctions = set(other.junctions)
    anchor_only = anchor_junctions - other_junctions
    other_only = other_junctions - anchor_junctions

    rows = []
    for label, junctions in (('anchor_only', anchor_only), ('other_only', other_only)):
        for junction in junctions:
            scores = score_junction(models, fasta, junction, context)
            rows.append({
                'label': label,
                'junction': repr(junction),
                'donor_position': junction.donor.coordinate,
                'acceptor_position': junction.acceptor.coordinate,
                'donor_prob': scores['donor_prob'],
                'acceptor_prob': scores['acceptor_prob'],
            })
    return rows
