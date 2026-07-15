"""
Tests for biosurfer.plots.interactive. Uses the same hand-built
Transcript/Exon/ORF/Protein fixture pattern as test_nmd.py -- no DB
round-trip needed. `plotly` is already a hard biosurfer dependency, so
these should run wherever biosurfer itself runs (unlike the optional-
dependency modules' tests).
"""
import pytest

from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import ORF, Exon, Gene, Protein, Transcript

from biosurfer.plots.interactive import build_interactive_isoform_figure


def _make_transcript(name, strand, exon_specs, orf_tx_start=None, orf_tx_stop=None):
    """
    exon_specs: list of (tx_start, tx_stop, g_start, g_stop). When an ORF
    range is given, the nucleotide sequence there is real codons (ATG +
    GCT*n + TAA) matching the placeholder protein sequence exactly, so
    Protein.residues' amino-acid-to-nucleotide linking (which translates
    and pattern-matches the sequence) actually succeeds -- codon
    coordinates are needed for the start/stop-codon markers.
    """
    exons = [Exon(transcript_start=a, transcript_stop=b, start=c, stop=d) for a, b, c, d in exon_specs]
    total_len = exon_specs[-1][1]
    orfs = []
    if orf_tx_start is not None:
        orf_len = orf_tx_stop - orf_tx_start + 1
        assert orf_len % 3 == 0 and orf_len >= 6
        n_ala = orf_len // 3 - 2
        orf_seq = 'ATG' + 'GCT' * n_ala + 'TAA'
        protein = Protein(sequence='M' + 'A' * n_ala + '*')
        orfs = [ORF(transcript_start=orf_tx_start, transcript_stop=orf_tx_stop, protein=protein)]
        sequence = ('A' * (orf_tx_start - 1)) + orf_seq + ('A' * (total_len - orf_tx_stop))
    else:
        sequence = 'A' * total_len
    return Transcript(
        name=name, gene=Gene(name='TEST', chromosome_id='chr1'),
        strand=strand, sequence=sequence, exons=exons, orfs=orfs,
    )


def test_build_figure_raises_on_empty_list():
    with pytest.raises(ValueError, match='No transcripts'):
        build_interactive_isoform_figure([])


def test_build_figure_raises_on_mixed_strands():
    plus_tx = _make_transcript('PLUS', Strand.PLUS, [(1, 30, 1000, 1029)])
    minus_tx = _make_transcript('MINUS', Strand.MINUS, [(1, 30, 2000, 2029)])
    with pytest.raises(ValueError, match='different strands'):
        build_interactive_isoform_figure([plus_tx, minus_tx])


def test_build_figure_coding_transcript_has_cds_and_utr_and_codon_markers():
    # exon1 tx1-60 g1000-1059; ORF tx11-58 (well within exon1, leaves 5' and 3' UTR)
    tx = _make_transcript('TX1', Strand.PLUS, [(1, 60, 1000, 1059)], orf_tx_start=11, orf_tx_stop=58)
    fig = build_interactive_isoform_figure([tx])

    hover_texts = [trace.text for trace in fig.data if getattr(trace, 'text', None)]
    assert any('(CDS)' in t for t in hover_texts)
    assert any('(UTR)' in t for t in hover_texts)
    assert any('start codon' in t for t in hover_texts)
    # last residue's amino acid needs to actually be a stop codon to get a stop marker;
    # our placeholder sequence ends in '*' (AminoAcid.STOP), so it should appear too
    assert any('stop codon' in t for t in hover_texts)


def test_build_figure_non_coding_transcript_has_no_codon_markers():
    tx = _make_transcript('NC', Strand.PLUS, [(1, 40, 1000, 1039)])
    fig = build_interactive_isoform_figure([tx])
    hover_texts = [trace.text for trace in fig.data if getattr(trace, 'text', None)]
    assert all('codon' not in t for t in hover_texts)
    assert any('(UTR)' in t for t in hover_texts)
    assert not any('(CDS)' in t for t in hover_texts)


def test_build_figure_multi_isoform_yaxis_labels():
    tx1 = _make_transcript('TX1', Strand.PLUS, [(1, 30, 1000, 1029)])
    tx2 = _make_transcript('TX2', Strand.PLUS, [(1, 30, 1000, 1029)])
    fig = build_interactive_isoform_figure([tx1, tx2])
    assert list(fig.layout.yaxis.ticktext) == ['TX1', 'TX2']
