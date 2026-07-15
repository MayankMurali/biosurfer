"""Gene-level orchestration and TSV report generation for `biosurfer analyze_disorder`."""
from pathlib import Path

import pandas as pd

from biosurfer.core.constants import APPRIS
from biosurfer.core.models.biomolecules import Gene

from biosurfer.analysis.disorder.predict import compare_isoform_disorder


def _select_anchor(transcripts):
    coding = [t for t in transcripts if t.protein]
    if not coding:
        return None
    return sorted(
        coding,
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True,
    )[0]


def run_disorder_report(session, gene_name: str, output_dir: Path) -> Path:
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
        try:
            rows = compare_isoform_disorder(anchor, other)
        except Exception as e:
            print(f'Skipping {other.name}: {e}')
            continue
        for row in rows:
            row['gene'] = gene.name
            row['anchor'] = anchor.name
            row['other'] = other.name
        records.extend(rows)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_disorder.tsv'
    pd.DataFrame(records).to_csv(output_file, sep='\t', index=False)
    print(f'Disorder report saved to: {output_file}')
    return output_file
