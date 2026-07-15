"""
Tests for biosurfer.analysis.splicing.

Keras/tensorflow/spliceai aren't installed in this environment, so
`score_position`/`load_models` are exercised via monkeypatching rather than
real inference -- these tests check the orchestration logic (window
fetching/reverse-complement, model-output indexing/averaging, junction
set-difference labeling), not SpliceAI's actual predictions. `_fetch_window`
is tested against a real pysam-indexed FASTA (pysam is a hard dependency).
"""
import pytest
import pysam

from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import Exon, Gene, Transcript

import biosurfer.analysis.splicing.predict as predict_module
from biosurfer.analysis.splicing.predict import (
    _fetch_window,
    compare_isoform_splice_sites,
    score_position,
)


def _write_fasta(tmp_path, name, sequence):
    fasta_path = tmp_path / f'{name}.fa'
    with open(fasta_path, 'w') as f:
        f.write(f'>{name}\n{sequence}\n')
    pysam.faidx(str(fasta_path))
    return pysam.FastaFile(str(fasta_path))


def test_fetch_window_plus_strand_exact_slice(tmp_path):
    seq = ''.join('ACGT'[i % 4] for i in range(200))  # 200bp deterministic sequence
    fasta = _write_fasta(tmp_path, 'chrT', seq)
    window = _fetch_window(fasta, 'chrT', Strand.PLUS, position=100, context=10)
    # context=10 -> 11bp window (5 each side + the position itself), centered on 1-based position 100
    assert len(window) == 11
    assert window == seq[100 - 1 - 5:100 - 1 + 6]
    fasta.close()


def test_fetch_window_minus_strand_is_reverse_complemented(tmp_path):
    seq = ''.join('ACGT'[i % 4] for i in range(200))
    fasta = _write_fasta(tmp_path, 'chrT', seq)
    plus = _fetch_window(fasta, 'chrT', Strand.PLUS, position=100, context=10)
    minus = _fetch_window(fasta, 'chrT', Strand.MINUS, position=100, context=10)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    assert minus == ''.join(complement[b] for b in reversed(plus))
    fasta.close()


def test_fetch_window_pads_with_n_past_contig_end(tmp_path):
    seq = 'ACGT' * 5  # 20bp
    fasta = _write_fasta(tmp_path, 'chrT', seq)
    window = _fetch_window(fasta, 'chrT', Strand.PLUS, position=2, context=10)
    # position=2, half=5 -> would need bases at 1-based positions -3..7; only 1..7 exist
    assert len(window) == 11
    assert window.startswith('NNNN')  # 4bp of padding for the out-of-bounds left side
    fasta.close()


class _FakeModel:
    def __init__(self, acceptor_val, donor_val):
        self.acceptor_val = acceptor_val
        self.donor_val = donor_val

    def predict(self, x):
        import numpy as np
        # shape (1, 1, 3): [neither, acceptor, donor]
        return np.array([[[1 - self.acceptor_val - self.donor_val, self.acceptor_val, self.donor_val]]])


def test_score_position_averages_model_ensemble(tmp_path, monkeypatch):
    seq = 'A' * 50
    fasta = _write_fasta(tmp_path, 'chrT', seq)
    monkeypatch.setattr(predict_module, '_SPLICEAI_AVAILABLE', True)
    monkeypatch.setattr(predict_module, 'one_hot_encode', lambda s: __import__('numpy').zeros((len(s), 4)), raising=False)
    models = [_FakeModel(0.2, 0.1), _FakeModel(0.4, 0.3)]
    result = score_position(models, fasta, 'chrT', Strand.PLUS, position=25, context=10)
    assert result['acceptor_prob'] == pytest.approx(0.3)  # mean(0.2, 0.4)
    assert result['donor_prob'] == pytest.approx(0.2)     # mean(0.1, 0.3)
    fasta.close()


def test_require_spliceai_raises_clear_error(monkeypatch):
    monkeypatch.setattr(predict_module, '_SPLICEAI_AVAILABLE', False)
    with pytest.raises(ImportError, match=r'biosurfer\[splicing\]'):
        predict_module.load_models()


def _make_transcript(name, strand, exon_specs):
    """exon_specs: list of (tx_start, tx_stop, g_start, g_stop)"""
    exons = [Exon(transcript_start=a, transcript_stop=b, start=c, stop=d) for a, b, c, d in exon_specs]
    total_len = exon_specs[-1][1]
    return Transcript(
        name=name, gene=Gene(name='TEST', chromosome_id='chr1'),
        strand=strand, sequence='A' * total_len, exons=exons, orfs=[],
    )


def test_compare_isoform_splice_sites_labels_unique_junctions(monkeypatch):
    anchor = _make_transcript('ANCHOR', Strand.PLUS, [
        (1, 30, 1000, 1029),
        (31, 60, 2000, 2029),
        (61, 90, 3000, 3029),
    ])
    # exon-skipping isoform: same first exon, jumps straight to anchor's third exon
    other = _make_transcript('OTHER', Strand.PLUS, [
        (1, 30, 1000, 1029),
        (31, 60, 3000, 3029),
    ])

    assert len(anchor.junctions) == 2
    assert len(other.junctions) == 1

    calls = []

    def fake_score_junction(models, fasta, junction, context=predict_module.DEFAULT_CONTEXT):
        calls.append(junction)
        return {'donor_prob': 0.5, 'acceptor_prob': 0.5}

    monkeypatch.setattr(predict_module, 'score_junction', fake_score_junction)
    rows = compare_isoform_splice_sites(models=[], fasta=None, anchor=anchor, other=other)

    labels = {row['label'] for row in rows}
    assert labels == {'anchor_only', 'other_only'}
    assert sum(1 for row in rows if row['label'] == 'anchor_only') == 2
    assert sum(1 for row in rows if row['label'] == 'other_only') == 1
    assert len(calls) == 3
