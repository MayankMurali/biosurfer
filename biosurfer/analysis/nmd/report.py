"""
Gene-level orchestration and TSV report generation for the three NMD
report modes exposed by ``biosurfer analyze_nmd``. Structurally mirrors
``biosurfer.analysis.genetics_analyzer`` (anchor selection, pandas
DataFrame -> TSV) rather than inventing a new reporting convention.
"""
from pathlib import Path

import pandas as pd
import pysam

from biosurfer.core.constants import APPRIS
from biosurfer.core.models.biomolecules import Gene

from biosurfer.analysis.nmd.predict import (
    compare_isoform_nmd,
    get_nmd_status,
    predict_variant_nmd_effect,
)


def _select_anchor(transcripts):
    """Same anchor-selection heuristic as genetics_analyzer.analyze_nterm_risk:
    prefer the APPRIS principal isoform, then the longest protein."""
    coding = [t for t in transcripts if t.protein]
    if not coding:
        return None
    return sorted(
        coding,
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True,
    )[0]


def run_transcript_nmd_report(session, gene_name: str, output_dir: Path) -> Path:
    """Mode 1: per-transcript NMD status for every transcript of a gene."""
    gene = Gene.from_name(session, gene_name)
    if not gene:
        print(f'Gene {gene_name} not found in database.')
        return None

    rows = []
    for transcript in gene.transcripts:
        rows.append({
            'gene': gene.name,
            'transcript': transcript.name,
            'accession': transcript.accession,
            'has_orf': transcript.primary_orf is not None,
            'nmd': get_nmd_status(transcript),
        })

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_nmd_transcripts.tsv'
    pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
    print(f'Per-transcript NMD report saved to: {output_file}')
    return output_file


def run_isoform_nmd_delta_report(session, gene_name: str, output_dir: Path) -> Path:
    """Mode 2: NMD status of every alternative isoform vs. the gene's anchor."""
    gene = Gene.from_name(session, gene_name)
    if not gene:
        print(f'Gene {gene_name} not found in database.')
        return None

    anchor = _select_anchor(gene.transcripts)
    if anchor is None:
        print(f'No coding transcripts found for {gene_name}.')
        return None
    print(f'Reference isoform: {anchor.name} ({anchor.accession})')

    others = [t for t in gene.transcripts if t is not anchor and t.protein]
    rows = [compare_isoform_nmd(anchor, other) for other in others]
    for row in rows:
        row['gene'] = gene.name

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_nmd_isoform_delta.tsv'
    pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
    print(f'Isoform NMD delta report saved to: {output_file}')
    return output_file


def run_variant_nmd_report(session, gene_name: str, vcf_path: str, output_dir: Path) -> Path:
    """
    Mode 3: for every SNP in ``vcf_path`` that falls within the gene's
    genomic span, predict whether it flips the NMD status of each of the
    gene's coding transcripts.

    Reads the VCF directly (via pysam, same library ``Database.load_genetics_data``
    already depends on) without writing anything to the Biosurfer database --
    this mode is read-only and doesn't touch ``core/`` or the schema at all.
    """
    gene = Gene.from_name(session, gene_name)
    if not gene:
        print(f'Gene {gene_name} not found in database.')
        return None

    coding_transcripts = [t for t in gene.transcripts if t.protein]
    if not coding_transcripts:
        print(f'No coding transcripts found for {gene_name}.')
        return None

    vcf = pysam.VariantFile(vcf_path)
    rows = []
    for record in vcf.fetch(gene.chromosome_id, gene.start, gene.stop):
        ref = record.ref
        for alt in record.alts or ():
            if len(ref) != 1 or len(alt) != 1:
                continue  # indels: out of scope for this first version, skip rather than misreport
            for transcript in coding_transcripts:
                effect = predict_variant_nmd_effect(transcript, gene.chromosome_id, record.pos, ref, alt)
                rows.append({
                    'gene': gene.name,
                    'transcript': transcript.name,
                    'chromosome': gene.chromosome_id,
                    'position': record.pos,
                    'ref': ref,
                    'alt': alt,
                    'rsid': record.id,
                    'applicable': effect.applicable,
                    'reason': effect.reason,
                    'reference_nmd': effect.reference_nmd,
                    'variant_nmd': effect.variant_nmd,
                    'delta': effect.delta,
                })

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_nmd_variant_effects.tsv'
    pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
    print(f'Variant NMD effect report saved to: {output_file}')
    return output_file
