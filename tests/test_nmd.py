"""
Tests for biosurfer.analysis.nmd. Uses plain, hand-constructed ORM objects
(same pattern as the `transcript` hypothesis strategy in test_models.py)
rather than the toy-GENCODE database fixture, since the module under test
only ever reads Transcript/ORF/Exon attributes -- no DB round-trip needed.

Coordinate/genomic-mapping choices here are deliberately derived from
construction parameters (lengths) rather than hardcoded, and cross-checked
by hand in code comments, to avoid silently-wrong magic numbers.
"""
from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import ORF, Exon, Gene, Protein, Transcript

from biosurfer.analysis.nmd.predict import compare_isoform_nmd, get_nmd_status, predict_variant_nmd_effect


def _make_gene():
    return Gene(name='TEST', chromosome_id='chr1')


def _make_transcript(name, strand, exon_lengths, genomic_starts, orf_tx_start=None, orf_tx_stop=None, sequence=None):
    """
    Build a single-ORF (or non-coding) transcript from a list of exon
    lengths (in transcript-coordinate order) and matching genomic start
    coordinates (one per exon; exon.stop is derived from length).
    """
    assert len(exon_lengths) == len(genomic_starts)
    total_len = sum(exon_lengths)
    if sequence is None:
        sequence = 'A' * total_len
    assert len(sequence) == total_len

    exons = []
    tx_pos = 1
    for length, g_start in zip(exon_lengths, genomic_starts):
        exons.append(Exon(
            transcript_start=tx_pos,
            transcript_stop=tx_pos + length - 1,
            start=g_start,
            stop=g_start + length - 1,
        ))
        tx_pos += length

    orfs = []
    if orf_tx_start is not None:
        orf_len = orf_tx_stop - orf_tx_start + 1
        protein = Protein(sequence='X' * ((orf_len // 3) - 1))
        orfs = [ORF(transcript_start=orf_tx_start, transcript_stop=orf_tx_stop, protein=protein)]

    return Transcript(
        name=name,
        gene=_make_gene(),
        strand=strand,
        sequence=sequence,
        exons=exons,
        orfs=orfs,
    )


def test_get_nmd_status_single_exon_transcript_is_never_nmd():
    # single exon -> no junction -> never NMD, regardless of stop position
    tx = _make_transcript('TX1', Strand.PLUS, exon_lengths=[90], genomic_starts=[1000],
                           orf_tx_start=11, orf_tx_stop=40)
    assert get_nmd_status(tx) is False


def test_get_nmd_status_stop_in_last_exon_is_not_nmd():
    # last_junction = 61 (exon2 tx_start); stop at tx 88-90 -> 61-90 = -29 < 50
    tx = _make_transcript('TX2', Strand.PLUS, exon_lengths=[60, 30], genomic_starts=[1000, 2000],
                           orf_tx_start=11, orf_tx_stop=90)
    assert get_nmd_status(tx) is False


def test_get_nmd_status_ptc_far_upstream_of_last_junction_is_nmd():
    # last_junction = 61; stop at tx 8 -> 61-8 = 53 >= 50
    tx = _make_transcript('TX3', Strand.PLUS, exon_lengths=[60, 30], genomic_starts=[1000, 2000],
                           orf_tx_start=3, orf_tx_stop=8)
    assert get_nmd_status(tx) is True


def test_get_nmd_status_non_coding_transcript_is_none():
    tx = _make_transcript('TX4', Strand.PLUS, exon_lengths=[50], genomic_starts=[1000])
    assert get_nmd_status(tx) is None


def test_compare_isoform_nmd_labels():
    anchor_not_nmd = _make_transcript('ANCHOR', Strand.PLUS, exon_lengths=[60, 30], genomic_starts=[1000, 2000],
                                       orf_tx_start=11, orf_tx_stop=90)
    other_nmd = _make_transcript('OTHER', Strand.PLUS, exon_lengths=[60, 30], genomic_starts=[1000, 2000],
                                  orf_tx_start=3, orf_tx_stop=8)
    result = compare_isoform_nmd(anchor_not_nmd, other_nmd)
    assert result['anchor_nmd'] is False
    assert result['other_nmd'] is True
    assert result['delta'] == 'triggers_nmd'

    # reciprocal
    result2 = compare_isoform_nmd(other_nmd, anchor_not_nmd)
    assert result2['delta'] == 'rescues_nmd'

    # unchanged
    result3 = compare_isoform_nmd(anchor_not_nmd, anchor_not_nmd)
    assert result3['delta'] == 'unchanged'

    # not_applicable: other has no ORF
    non_coding = _make_transcript('NC', Strand.PLUS, exon_lengths=[50], genomic_starts=[1000])
    result4 = compare_isoform_nmd(anchor_not_nmd, non_coding)
    assert result4['delta'] == 'not_applicable'


def _variant_test_transcript():
    """
    utr5(2nt) + ORF[ATG, CAA, 50x AAC filler, TAA](159nt) + utr3(4nt) = 165nt,
    split into exon1 (tx 1-60) and exon2 (tx 61-165).

    - Reference stop (real ORF stop, TAA at the end): tx_stop = 2 + 159 = 161.
      last_junction(61) - 161 = -100 < 50 -> reference NOT NMD.
    - Codon 2 (CAA, tx 6-8) is the mutation target: C->T at tx position 6
      turns CAA (Gln) into TAA (stop). New stop tx_coord = 8.
      last_junction(61) - 8 = 53 >= 50 -> mutant IS NMD-triggering.
    - tx position 6 is the first base of exon1 (genomic start=1000), so its
      genomic coordinate is 1000 + (6-1) = 1005.
    """
    utr5 = 'AA'
    orf_seq = 'ATG' + 'CAA' + 'AAC' * 50 + 'TAA'
    utr3 = 'GGGG'
    sequence = utr5 + orf_seq + utr3
    orf_tx_start = len(utr5) + 1  # 3
    orf_tx_stop = len(utr5) + len(orf_seq)  # 161
    total_len = len(sequence)  # 165
    exon1_len = 60
    exon2_len = total_len - exon1_len  # 105
    tx = _make_transcript(
        'VARTX', Strand.PLUS,
        exon_lengths=[exon1_len, exon2_len], genomic_starts=[1000, 2000],
        orf_tx_start=orf_tx_start, orf_tx_stop=orf_tx_stop, sequence=sequence,
    )
    return tx


def test_predict_variant_nmd_effect_triggers_nmd():
    tx = _variant_test_transcript()
    assert get_nmd_status(tx) is False  # sanity check on the unmutated reference

    effect = predict_variant_nmd_effect(tx, 'chr1', 1005, 'C', 'T')
    assert effect.applicable is True
    assert effect.reference_nmd is False
    assert effect.variant_nmd is True
    assert effect.delta == 'triggers_nmd'
    assert effect.variant_stop_tx_coord == 8


def test_predict_variant_nmd_effect_ref_allele_mismatch():
    tx = _variant_test_transcript()
    # tx position 6 is actually 'C' (first base of codon CAA); assert wrong ref is rejected
    effect = predict_variant_nmd_effect(tx, 'chr1', 1005, 'G', 'T')
    assert effect.applicable is False
    assert 'mismatch' in effect.reason


def test_predict_variant_nmd_effect_rejects_indels():
    tx = _variant_test_transcript()
    effect = predict_variant_nmd_effect(tx, 'chr1', 1005, 'C', 'TT')
    assert effect.applicable is False
    assert 'single-nucleotide' in effect.reason


def test_predict_variant_nmd_effect_no_orf():
    tx = _make_transcript('NC2', Strand.PLUS, exon_lengths=[50], genomic_starts=[1000])
    effect = predict_variant_nmd_effect(tx, 'chr1', 1010, 'A', 'T')
    assert effect.applicable is False
    assert 'no annotated ORF' in effect.reason
