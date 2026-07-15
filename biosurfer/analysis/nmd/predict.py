"""
Core NMD prediction logic.

Biosurfer's ``ORF`` class already implements the standard 50-nt rule for
nonsense-mediated decay (NMD) susceptibility:

    ORF.nmd  ==  (last_junction - transcript_stop) >= 50

where ``last_junction`` is the transcript-coordinate position of the start of
the last exon (i.e. the last exon-exon junction) and ``transcript_stop`` is
the last nucleotide of the stop codon, both in transcript coordinates. That
one-line formula already handles the two standard NMD-escape cases correctly
without special-casing them:

- single-exon transcripts have no junction, so ``last_junction`` collapses to
  0 and the difference is always negative -> ``nmd`` is ``False``;
- a stop codon in the last exon means ``transcript_stop >= last_junction``,
  so the difference is <= 0 -> ``nmd`` is ``False``.

This property existed in the codebase already but was never called from
anywhere else (confirmed via repo-wide grep before writing this module) --
this module surfaces it at three levels of increasing specificity:

1. ``get_nmd_status`` -- per-transcript lookup.
2. ``compare_isoform_nmd`` -- anchor-vs-alternative isoform delta, reusing
   the anchor/other vocabulary already established by
   ``biosurfer.analysis.genetics_analyzer``.
3. ``predict_variant_nmd_effect`` -- re-evaluates the same rule after
   applying a single-nucleotide substitution, to answer "does this specific
   variant flip this transcript's predicted NMD fate?"
"""
from typing import NamedTuple, Optional

from Bio.Seq import Seq

from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import Transcript


class VariantNMDEffect(NamedTuple):
    """Result of applying a single variant to a transcript and re-checking NMD status."""
    applicable: bool
    reason: str
    reference_nmd: Optional[bool]
    variant_nmd: Optional[bool]
    reference_stop_tx_coord: Optional[int]
    variant_stop_tx_coord: Optional[int]

    @property
    def delta(self) -> str:
        if not self.applicable or self.reference_nmd is None or self.variant_nmd is None:
            return 'not_applicable'
        if self.reference_nmd == self.variant_nmd:
            return 'unchanged'
        return 'triggers_nmd' if self.variant_nmd else 'rescues_nmd'


def get_nmd_status(transcript: Transcript) -> Optional[bool]:
    """
    Return ``transcript``'s primary ORF's NMD status, or ``None`` if the
    transcript has no annotated ORF (e.g. a non-coding transcript).
    """
    orf = transcript.primary_orf
    return orf.nmd if orf else None


def compare_isoform_nmd(anchor: Transcript, other: Transcript) -> dict:
    """
    Compare NMD susceptibility between an anchor (reference) isoform and an
    alternative isoform of the same gene -- the isoform-switch question:
    does this alternative splicing event newly trigger, or rescue from,
    NMD relative to the anchor?

    Returns a dict with both transcripts' NMD status and a ``delta`` label:
    ``'triggers_nmd'`` (anchor isn't NMD-susceptible, other is),
    ``'rescues_nmd'`` (anchor is, other isn't), ``'unchanged'`` (both the
    same), or ``'not_applicable'`` (either transcript has no ORF).
    """
    anchor_nmd = get_nmd_status(anchor)
    other_nmd = get_nmd_status(other)
    if anchor_nmd is None or other_nmd is None:
        delta = 'not_applicable'
    elif anchor_nmd == other_nmd:
        delta = 'unchanged'
    elif other_nmd:
        delta = 'triggers_nmd'
    else:
        delta = 'rescues_nmd'
    return {
        'anchor': anchor.name,
        'other': other.name,
        'anchor_nmd': anchor_nmd,
        'other_nmd': other_nmd,
        'delta': delta,
    }


def predict_variant_nmd_effect(transcript: Transcript, chromosome: str, position: int, ref: str, alt: str) -> VariantNMDEffect:
    """
    Predict whether a single-nucleotide substitution shifts ``transcript``'s
    NMD status, by re-translating its ORF from the variant-containing
    sequence and re-evaluating the same 50-nt-from-last-junction rule that
    ``ORF.nmd`` uses, against the new stop codon position.

    ``chromosome``/``position``/``ref``/``alt`` follow standard VCF
    convention: 1-based genomic coordinate, alleles given relative to the
    forward (+) genomic strand -- ``biosurfer``'s own ``GenomicVariant``
    model and ``Database.load_genetics_data`` use the same convention.

    SNPs only in this first version. Indels shift the reading frame in ways
    that need careful transcript-coordinate bookkeeping this function
    doesn't attempt yet -- they're explicitly rejected (``applicable=False``)
    rather than silently mishandled.
    """
    orf = transcript.primary_orf
    if orf is None:
        return VariantNMDEffect(False, 'transcript has no annotated ORF', None, None, None, None)
    if len(ref) != 1 or len(alt) != 1:
        return VariantNMDEffect(False, 'only single-nucleotide substitutions are supported', None, None, None, None)

    nt = transcript.get_nucleotide_from_coordinate(position)
    if nt is None:
        return VariantNMDEffect(False, 'variant position not covered by this transcript', None, None, None, None)
    tx_index = nt.position - 1  # 0-based index into transcript.sequence

    # transcript.sequence is mRNA-sense (5'->3'); for a minus-strand
    # transcript that's the reverse complement of the +strand genomic
    # sequence the VCF ref/alt alleles are given against, so complement the
    # alleles before comparing/substituting. If this assumption is ever
    # wrong for some edge case, the mismatch check below fails loud instead
    # of silently producing a wrong prediction.
    if transcript.strand is Strand.MINUS:
        expected_base = str(Seq(ref).complement())
        variant_base = str(Seq(alt).complement())
    else:
        expected_base = ref
        variant_base = alt

    tx_seq = transcript.sequence
    observed_base = tx_seq[tx_index]
    if observed_base.upper() != expected_base.upper():
        return VariantNMDEffect(
            False,
            f'reference allele mismatch: transcript has {observed_base!r} at this position, expected {expected_base!r}',
            None, None, None, None,
        )

    mutated_seq = tx_seq[:tx_index] + variant_base + tx_seq[tx_index + 1:]

    last_junction = transcript.exons[-1].transcript_start
    reference_nmd = orf.nmd
    reference_stop_tx_coord = orf.transcript_stop

    variant_stop_tx_coord = _find_first_inframe_stop(mutated_seq, orf.transcript_start)
    if variant_stop_tx_coord is None:
        # ORF runs off the end of the transcript without hitting a stop --
        # the 50-nt rule doesn't have a well-defined "transcript_stop" to
        # compare against here, so report this rather than guessing.
        return VariantNMDEffect(
            True, 'variant ORF has no in-frame stop codon before transcript end',
            reference_nmd, None, reference_stop_tx_coord, None,
        )

    variant_nmd = (last_junction - variant_stop_tx_coord) >= 50
    return VariantNMDEffect(
        True, 'ok', reference_nmd, variant_nmd, reference_stop_tx_coord, variant_stop_tx_coord,
    )


def _find_first_inframe_stop(sequence: str, orf_start_1based: int) -> Optional[int]:
    """
    Scan codons from ``orf_start_1based`` (1-based transcript coordinate)
    and return the 1-based transcript coordinate of the last nucleotide of
    the first in-frame stop codon, or ``None`` if translation runs off the
    end of ``sequence`` without hitting one.
    """
    cds = sequence[orf_start_1based - 1:]
    in_frame_len = len(cds) - len(cds) % 3
    protein = str(Seq(cds[:in_frame_len]).translate())
    stop_index = protein.find('*')
    if stop_index == -1:
        return None
    # stop_index is a 0-based amino acid index; that codon occupies
    # nucleotides [stop_index*3, stop_index*3 + 3) of `cds`
    return orf_start_1based + (stop_index + 1) * 3 - 1
