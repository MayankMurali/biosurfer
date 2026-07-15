"""
Tests for biosurfer.analysis.disorder.

`metapredict` isn't installed in this environment, so `compare_isoform_disorder`
is exercised against the real toy-GENCODE `session` fixture (CRYBG2 has
several real coding isoforms with genuine protein alignments) while
monkeypatching `score_disorder` with a deterministic stand-in -- this
tests the real alignment-block iteration/row-construction logic without
needing the optional dependency installed.
"""
import pytest

from biosurfer.core.models.biomolecules import Gene

import biosurfer.analysis.disorder.predict as predict_module
from biosurfer.analysis.disorder.predict import score_protein_range_disorder


def test_score_protein_range_disorder_basic_stats():
    scores = [0.1, 0.2, 0.9, 0.9, 0.3]
    result = score_protein_range_disorder(scores, range(1, 4))  # indices 1,2,3
    assert result['mean'] == pytest.approx((0.2 + 0.9 + 0.9) / 3)
    assert result['max'] == 0.9
    assert result['min'] == 0.2


def test_score_protein_range_disorder_empty_range_returns_none():
    result = score_protein_range_disorder([0.1, 0.2, 0.3], range(5, 5))
    assert result == {'mean': None, 'max': None, 'min': None}


def test_require_metapredict_raises_clear_error(monkeypatch):
    monkeypatch.setattr(predict_module, '_METAPREDICT_AVAILABLE', False)
    with pytest.raises(ImportError, match=r'biosurfer\[disorder\]'):
        predict_module.score_disorder('ACDEFG')


def test_compare_isoform_disorder_against_real_toy_gene(session, monkeypatch):
    gene = Gene.from_name(session, 'CRYBG2')
    coding = [t for t in gene.transcripts if t.protein]
    assert len(coding) >= 2, 'toy CRYBG2 fixture should have at least 2 coding isoforms'
    anchor, other = coding[0], coding[1]

    # deterministic stand-in: every residue scores 0.5, sized to each protein's real length
    def fake_score_disorder(sequence):
        return [0.5] * len(sequence)

    monkeypatch.setattr(predict_module, 'score_disorder', fake_score_disorder)

    rows = predict_module.compare_isoform_disorder(anchor, other)
    assert len(rows) > 0
    for row in rows:
        if row['anchor_disorder_mean'] is not None:
            assert row['anchor_disorder_mean'] == pytest.approx(0.5)
        if row['other_disorder_mean'] is not None:
            assert row['other_disorder_mean'] == pytest.approx(0.5)
        assert 'pblock_category' in row and 'cblock_category' in row
