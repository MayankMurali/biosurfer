"""
Tests for biosurfer.analysis.conservation.

Conservation (bigWig) tests use a duck-typed fake in place of a real
pyBigWig handle -- `score_genomic_range`/`score_protein_range_conservation`
never touch the `pyBigWig` module directly (only `open_bigwig` does), so
they're testable without the optional dependency installed. AlphaMissense
tests build a real tiny bgzipped+tabix-indexed TSV via pysam (already a
hard biosurfer dependency) rather than mocking it.
"""
import gzip

import pysam
import pytest

from biosurfer.core.constants import Strand
from biosurfer.core.models.biomolecules import ORF, Exon, Gene, Protein, Transcript

from biosurfer.analysis.conservation.predict import (
    open_bigwig,
    query_alphamissense,
    score_genomic_range,
    score_protein_range_conservation,
)


class FakeBigWig:
    """Duck-typed stand-in for a pyBigWig handle: differentiates scores by
    genomic region so multi-exon aggregation can be exercised without a
    real bigWig file."""

    def stats(self, chrom, start, stop, type='mean'):
        vals = {'mean': 0.2, 'max': 0.5, 'min': 0.0} if start < 1500 else {'mean': 0.8, 'max': 1.0, 'min': 0.6}
        return [vals[type]]


def _make_transcript_with_orf(exon_lengths, genomic_starts, orf_tx_start, orf_tx_stop):
    total_len = sum(exon_lengths)
    exons = []
    tx_pos = 1
    for length, g_start in zip(exon_lengths, genomic_starts):
        exons.append(Exon(
            transcript_start=tx_pos, transcript_stop=tx_pos + length - 1,
            start=g_start, stop=g_start + length - 1,
        ))
        tx_pos += length

    orf_len = orf_tx_stop - orf_tx_start + 1
    protein = Protein(sequence='X' * (orf_len // 3))
    orf = ORF(transcript_start=orf_tx_start, transcript_stop=orf_tx_stop, protein=protein)

    return Transcript(
        name='TX', gene=Gene(name='TEST', chromosome_id='chr1'),
        strand=Strand.PLUS, sequence='A' * total_len, exons=exons, orfs=[orf],
    )


def test_score_genomic_range_swaps_reversed_bounds():
    bw = FakeBigWig()
    forward = score_genomic_range(bw, 'chr1', 1010, 1024)
    reversed_ = score_genomic_range(bw, 'chr1', 1024, 1010)
    assert forward == reversed_


def test_score_protein_range_conservation_single_exon():
    # exon1 tx1-60 genomic 1000-1059; ORF tx11-118 (spans both exons, but
    # this protein range stays within exon1)
    tx = _make_transcript_with_orf([60, 60], [1000, 2000], orf_tx_start=11, orf_tx_stop=118)
    bw = FakeBigWig()
    result = score_protein_range_conservation(bw, tx, range(0, 5))
    assert result['error'] is None
    assert result['mean'] == 0.2
    assert result['max'] == 0.5
    assert result['min'] == 0.0


def test_score_protein_range_conservation_spans_two_exons():
    tx = _make_transcript_with_orf([60, 60], [1000, 2000], orf_tx_start=11, orf_tx_stop=118)
    bw = FakeBigWig()
    # aa 15-19 (0-based half-open [15,20)) straddles the exon1/exon2 boundary
    result = score_protein_range_conservation(bw, tx, range(15, 20))
    assert result['error'] is None
    assert result['mean'] == pytest.approx(0.5)  # mean of 0.2 and 0.8
    assert result['max'] == 1.0
    assert result['min'] == 0.0


def test_score_protein_range_conservation_no_protein_returns_error():
    # non-coding transcript
    non_coding = Transcript(
        name='NC', gene=Gene(name='TEST', chromosome_id='chr1'),
        strand=Strand.PLUS, sequence='A' * 50,
        exons=[Exon(transcript_start=1, transcript_stop=50, start=1000, stop=1049)],
        orfs=[],
    )
    result = score_protein_range_conservation(FakeBigWig(), non_coding, range(0, 5))
    assert result['error'] is not None
    assert result['mean'] is None


def test_open_bigwig_raises_clear_import_error_when_pybigwig_missing(monkeypatch):
    import biosurfer.analysis.conservation.predict as predict_module
    monkeypatch.setattr(predict_module, 'pyBigWig', None)
    with pytest.raises(ImportError, match='biosurfer\\[conservation\\]'):
        open_bigwig('/some/path.bw')


def _write_alphamissense_tsv(path, rows):
    """rows: list of (chrom, pos, ref, alt, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class)"""
    plain_path = str(path) + '.plain'
    with open(plain_path, 'w') as f:
        f.write('#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')
    pysam.tabix_compress(plain_path, str(path), force=True)
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, meta_char='#', force=True)


def test_query_alphamissense_returns_overlapping_rows(tmp_path):
    tsv_path = tmp_path / 'AlphaMissense_test.tsv.gz'
    _write_alphamissense_tsv(tsv_path, [
        ('chr1', 1005, 'C', 'T', 'hg38', 'P12345', 'ENST00000001', 'Q2*', 0.95, 'likely_pathogenic'),
        ('chr1', 1050, 'G', 'A', 'hg38', 'P12345', 'ENST00000001', 'A17T', 0.10, 'likely_benign'),
        ('chr1', 5000, 'A', 'C', 'hg38', 'P99999', 'ENST00000002', 'K3Q', 0.50, 'ambiguous'),
    ])
    hits = query_alphamissense(str(tsv_path), 'chr1', 1000, 1024)
    assert len(hits) == 1
    assert hits[0]['position'] == 1005
    assert hits[0]['am_pathogenicity'] == pytest.approx(0.95)
    assert hits[0]['am_class'] == 'likely_pathogenic'


def test_query_alphamissense_no_hits_in_range(tmp_path):
    tsv_path = tmp_path / 'AlphaMissense_test2.tsv.gz'
    _write_alphamissense_tsv(tsv_path, [
        ('chr1', 5000, 'A', 'C', 'hg38', 'P99999', 'ENST00000002', 'K3Q', 0.50, 'ambiguous'),
    ])
    hits = query_alphamissense(str(tsv_path), 'chr1', 1000, 1024)
    assert hits == []


def test_query_alphamissense_unknown_contig_returns_empty(tmp_path):
    tsv_path = tmp_path / 'AlphaMissense_test3.tsv.gz'
    _write_alphamissense_tsv(tsv_path, [
        ('chr1', 1005, 'C', 'T', 'hg38', 'P12345', 'ENST00000001', 'Q2*', 0.95, 'likely_pathogenic'),
    ])
    hits = query_alphamissense(str(tsv_path), 'chr99', 1, 100)
    assert hits == []


def test_query_alphamissense_missing_file_raises_clear_error(tmp_path):
    with pytest.raises(ValueError, match='tabix'):
        query_alphamissense(str(tmp_path / 'does_not_exist.tsv.gz'), 'chr1', 1, 100)
