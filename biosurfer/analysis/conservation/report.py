"""
Gene-level orchestration and TSV report generation for
``biosurfer analyze_conservation``. Structurally mirrors
``biosurfer.analysis.genetics_analyzer.generate_genetics_report`` --
same anchor-selection heuristic, same hybrid-alignment-block iteration,
same DataFrame -> TSV pattern -- just with conservation/AlphaMissense
columns instead of GWAS columns.
"""
from pathlib import Path
from typing import Optional

import pandas as pd

from biosurfer.core.alignments import CodonAlignment, ProteinAlignment, TranscriptAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.models.biomolecules import Gene
from biosurfer.core.splice_events import get_event_code

from biosurfer.analysis.conservation.predict import (
    query_alphamissense_for_protein_range,
    score_protein_range_conservation,
)


def _select_anchor(transcripts):
    coding = [t for t in transcripts if t.protein]
    if not coding:
        return None
    return sorted(
        coding,
        key=lambda t: (t.appris == APPRIS.PRINCIPAL, t.protein.length),
        reverse=True,
    )[0]


def run_conservation_report(
    session, gene_name: str, output_dir: Path,
    bw=None, alphamissense_tabix: Optional[str] = None,
) -> Path:
    """
    For every alternative isoform's protein/codon alignment blocks against
    the gene's anchor isoform, report conservation (if ``bw`` is an open
    pyBigWig handle) and/or AlphaMissense pathogenicity (if
    ``alphamissense_tabix`` is given) over each block's genomic span.

    At least one of ``bw``/``alphamissense_tabix`` should be provided, but
    neither is required by this function itself -- the CLI layer enforces
    that at least one was requested.
    """
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
            tx_aln = TranscriptAlignment.from_transcripts(anchor, other)
            cd_aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
            pr_aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
        except Exception:
            continue

        for p, pblock in enumerate(pr_aln.blocks):
            cblocks = pr_aln.pblock_to_cblocks.get(pblock, [])
            for c, cblock in enumerate(cblocks):
                tblock = cd_aln.cblock_to_tblock.get(cblock)
                events = tx_aln.block_to_events.get(tblock, ())

                target_tx, target_range = None, None
                if cblock.anchor_range:
                    target_tx, target_range = anchor, cblock.anchor_range
                elif cblock.other_range:
                    target_tx, target_range = other, cblock.other_range

                row = {
                    'gene': gene.name,
                    'anchor': anchor.name,
                    'other': other.name,
                    'pblock_number': p,
                    'pblock_category': pblock.category.name,
                    'cblock_number': c,
                    'cblock_category': cblock.category.name,
                    'events': get_event_code(events),
                }

                if target_tx and target_range:
                    if bw is not None:
                        cons = score_protein_range_conservation(bw, target_tx, target_range)
                        row['conservation_mean'] = cons['mean']
                        row['conservation_max'] = cons['max']
                        row['conservation_min'] = cons['min']
                    if alphamissense_tabix:
                        am_hits = query_alphamissense_for_protein_range(alphamissense_tabix, target_tx, target_range)
                        row['alphamissense_variant_count'] = len(am_hits)
                        if am_hits:
                            row['alphamissense_max_pathogenicity'] = max(h['am_pathogenicity'] for h in am_hits)
                        else:
                            row['alphamissense_max_pathogenicity'] = None

                records.append(row)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{gene.name}_conservation.tsv'
    pd.DataFrame(records).to_csv(output_file, sep='\t', index=False)
    print(f'Conservation report saved to: {output_file}')
    return output_file
