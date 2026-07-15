"""
Gene-level orchestration and TSV report generation for
``biosurfer analyze_splicing``. Same anchor-selection heuristic as
``genetics_analyzer.py``/``analysis.nmd.report``.
"""
from pathlib import Path

import pandas as pd

from biosurfer.core.constants import APPRIS
from biosurfer.core.models.biomolecules import Gene

from biosurfer.analysis.splicing.predict import compare_isoform_splice_sites


def _select_anchor(transcripts):
    coding = [t for t in transcripts if t.protein]
    if not coding:
        return None
    return sorted(
        coding,
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True,
    )[0]


def run_splicing_report(session, gene_name: str, output_dir: Path, models, fasta) -> Path:
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
    records = []
    for other in others:
        rows = compare_isoform_splice_sites(models, fasta, anchor, other)
        for row in rows:
            row['gene'] = gene.name
            row['anchor'] = anchor.name
            row['other'] = other.name
        records.extend(rows)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_splice_site_strength.tsv'
    pd.DataFrame(records).to_csv(output_file, sep='\t', index=False)
    print(f'Splice-site-strength report saved to: {output_file}')
    return output_file
