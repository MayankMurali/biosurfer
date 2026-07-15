"""
Interactive Plotly isoform-structure viewer, alongside the existing static
matplotlib ``IsoformPlot`` (``biosurfer.plots.plotting``). Purely additive:
reads the same ``Transcript``/``Exon``/``ORF`` data and renders a new kind
of output (a zoomable/pannable/hoverable HTML figure) instead of a static
PNG. Doesn't modify ``plots/plotting.py`` or ``plots/canvas.py`` at all.

v1 scope: exon/CDS/UTR structure with start/stop-codon markers and
per-exon hover tooltips, one row per isoform -- the same information
``IsoformPlot.draw_isoform`` shows, re-rendered as an interactive figure.
Domain/pblock/cblock overlays aren't included yet (flagged as a natural
follow-up, not silently dropped).

``plotly`` is already a base biosurfer dependency (``setup.py``
``install_requires``), so this introduces no new dependency at all.
"""
from typing import List

import plotly.graph_objects as go

from biosurfer.core.constants import AminoAcid, Strand
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.plots.plotting import TRANSCRIPT_COLORS


def _exon_rectangle(exon_start: int, exon_stop: int, track: int, height: float, color: str, hover_text: str) -> go.Scatter:
    y0, y1 = track - height / 2, track + height / 2
    return go.Scatter(
        x=[exon_start, exon_stop, exon_stop, exon_start, exon_start],
        y=[y0, y0, y1, y1, y0],
        fill='toself', fillcolor=color, line=dict(color='black', width=1),
        mode='lines', hoverinfo='text', text=hover_text, showlegend=False,
    )


def build_interactive_isoform_figure(transcripts: List['Transcript']) -> go.Figure:
    """
    Build an interactive Plotly figure showing exon/CDS/UTR structure for
    each transcript in ``transcripts``, one row per isoform, hoverable.
    """
    transcripts = [tx for tx in transcripts if tx is not None]
    if not transcripts:
        raise ValueError('No transcripts to plot')

    strands = {tx.strand for tx in transcripts}
    if len(strands) > 1:
        raise ValueError("Can't plot isoforms from different strands")
    strand = next(iter(strands))

    fig = go.Figure()
    height = 0.4

    for track, tx in enumerate(transcripts):
        color_dark, color_light = TRANSCRIPT_COLORS.get(type(tx), TRANSCRIPT_COLORS[None])

        fig.add_trace(go.Scatter(
            x=[tx.start, tx.stop], y=[track, track],
            mode='lines', line=dict(color='gray', width=1.5),
            hoverinfo='skip', showlegend=False,
        ))

        orf = tx.primary_orf
        for exon in tx.exons:
            cds_start = cds_stop = None
            if orf:
                overlap_start = max(exon.start, orf.start)
                overlap_stop = min(exon.stop, orf.stop)
                if overlap_start <= overlap_stop:
                    cds_start, cds_stop = overlap_start, overlap_stop

            hover = f'{tx.name} exon {exon.position}<br>{exon.start}-{exon.stop}'
            if cds_start is not None:
                if exon.start < cds_start:
                    fig.add_trace(_exon_rectangle(exon.start, cds_start, track, height, color_light, hover + ' (UTR)'))
                if exon.stop > cds_stop:
                    fig.add_trace(_exon_rectangle(cds_stop, exon.stop, track, height, color_light, hover + ' (UTR)'))
                fig.add_trace(_exon_rectangle(cds_start, cds_stop, track, height, color_dark, hover + ' (CDS)'))
            else:
                fig.add_trace(_exon_rectangle(exon.start, exon.stop, track, height, color_light, hover + ' (UTR)'))

        if orf and orf.protein and orf.protein.residues:
            first_res, last_res = orf.protein.residues[0], orf.protein.residues[-1]
            if first_res.amino_acid is AminoAcid.MET:
                fig.add_trace(go.Scatter(
                    x=[first_res.codon[0].coordinate], y=[track], mode='markers',
                    marker=dict(color='lime', symbol='line-ns', size=12, line=dict(width=2)),
                    hoverinfo='text', text=f'{tx.name} start codon', showlegend=False,
                ))
            if last_res.amino_acid is AminoAcid.STOP:
                fig.add_trace(go.Scatter(
                    x=[last_res.codon[2].coordinate], y=[track], mode='markers',
                    marker=dict(color='red', symbol='line-ns', size=12, line=dict(width=2)),
                    hoverinfo='text', text=f'{tx.name} stop codon', showlegend=False,
                ))

    fig.update_yaxes(
        tickmode='array', tickvals=list(range(len(transcripts))),
        ticktext=[tx.name for tx in transcripts], autorange='reversed',
    )
    fig.update_xaxes(
        title='Genomic position',
        autorange='reversed' if strand is Strand.MINUS else True,
    )
    fig.update_layout(showlegend=False, height=max(200, 60 * len(transcripts)))
    return fig


def save_interactive_plot(fig: go.Figure, output_path) -> None:
    fig.write_html(str(output_path))
