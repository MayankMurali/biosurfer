# functions to create different visualizations of isoforms/clones/domains/muts
from itertools import groupby, tee
from operator import attrgetter
from typing import (TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple)
from warnings import filterwarnings, warn

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from biosurfer.core.alignments import CodonAlignment, ProjectedFeature
from biosurfer.core.algorithms import generate_subtracks
from biosurfer.core.constants import (FRAMESHIFT, SPLIT_CODON, AminoAcid,
                                      CodonAlignmentCategory, FeatureType,
                                      SequenceAlignmentCategory, Strand)
from biosurfer.core.collections_utils import ExceptionLogger
from biosurfer.core.models.biomolecules import (GencodeTranscript,
                                                PacBioTranscript, Transcript)
from biosurfer.core.splice_events import (AcceptorSpliceEvent,
                                          DonorSpliceEvent, ExonBypassEvent,
                                          ExonSpliceEvent, IntronSpliceEvent)
from biosurfer.plots.canvas import PlotCanvas, TableColumn
from brokenaxes import BrokenAxes
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from more_itertools import first, last, only

if TYPE_CHECKING:
    from biosurfer.core.alignments import ProteinAlignmentBlock, CodonAlignmentBlock
    from biosurfer.core.models.biomolecules import Protein

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

TRANSCRIPT_SOURCE = {
    'gencodetranscript': 'GENCODE',
    'pacbiotranscript': 'PacBio'
}
def get_transcript_source(tx: 'Transcript'):
    return TRANSCRIPT_SOURCE.get(tx.type, '')

# colors for different transcript types
TRANSCRIPT_COLORS = {
    None: ('#404040', '#777777'),
    GencodeTranscript: ('#343553', '#5D5E7C'),
    PacBioTranscript: ('#61374D', '#91677D')
}

# colors for transcript events
EVENT_COLORS = {
    IntronSpliceEvent: '#e69138',
    DonorSpliceEvent: '#6aa84f',
    AcceptorSpliceEvent: '#674ea7',
    ExonSpliceEvent: '#3d85c6',
    ExonBypassEvent: '#bebebe',
}

# alpha values for different absolute reading frames
ABS_FRAME_ALPHA = {0: 1.0, 1: 0.45, 2: 0.15}

# hatching styles for different relative frameshifts
REL_FRAME_STYLE = {
    CodonAlignmentCategory.FRAME_AHEAD: '////',
    CodonAlignmentCategory.FRAME_BEHIND: 'xxxx'
}

# colors for CodonAlignmentBlocks
CBLOCK_COLORS = {
    CodonAlignmentCategory.TRANSLATED: '#9bf3ff',
    CodonAlignmentCategory.INSERTION: '#05e0ff',
    CodonAlignmentCategory.FRAME_AHEAD: '#fff099',
    CodonAlignmentCategory.FRAME_BEHIND: '#ffd700',
    CodonAlignmentCategory.UNTRANSLATED: '#ff99ce',
    CodonAlignmentCategory.DELETION: '#ff0082',
    CodonAlignmentCategory.EDGE: '#8270c1',
    CodonAlignmentCategory.COMPLEX: '#aaaaaa'
}

PBLOCK_COLORS = {
    SequenceAlignmentCategory.DELETION: '#FF0082',
    SequenceAlignmentCategory.INSERTION: '#05E0FF',
    SequenceAlignmentCategory.SUBSTITUTION: '#FFD700'
}

FEATURE_COLORS = {
    'MobiDB': '#AAAAAA'
}


class IsoformPlot(PlotCanvas):
    """Encapsulates methods for drawing one or more isoforms aligned to the same genomic x-axis."""
    def __init__(self, transcripts: Iterable['Transcript'], columns: Dict[str, TableColumn] = None, **kwargs):
        self.transcripts: List['Transcript'] = list(transcripts)  # list of orf objects to be drawn
        gene = {tx.gene for tx in filter(None, self.transcripts)}
        if len(gene) > 1:
            raise ValueError(f'Found isoforms from multiple genes: {", ".join(g.name for g in gene)}')
        strand = only(
            {tx.strand for tx in filter(None, self.transcripts)},
            too_long = ValueError("Can't plot isoforms from different strands")
        )
        self._columns: Dict[str, TableColumn] = {'Source': get_transcript_source} | (columns if columns else dict())
        super().__init__(strand=strand, num_tracks=len(self.transcripts), **kwargs)
        self.reset_xlims()

    def reset_xlims(self):
        """Set xlims automatically based on exons in isoforms."""
        space = self.opts.intron_spacing
        self.xlims = tuple((exon.start - space, exon.stop + space) for tx in filter(None, self.transcripts) for exon in tx.exons)

    def draw_isoform(self, tx: 'Transcript', track: int):
        """Plot a single isoform in the given track."""
        start, stop = tx.start, tx.stop
        align_start, align_stop = 'right', 'left'
        if self.strand is Strand.MINUS:
            align_start, align_stop = align_stop, align_start

        # plot intron line
        self.draw_region(
            track,
            start = start,
            stop = stop,
            type = 'line',
            linewidth = 1.5,
            color = 'gray',
            zorder = 1.5
        )

        # plot exons
        utr_kwargs = {
            'type': 'rect',
            'edgecolor': 'k',
            'facecolor': TRANSCRIPT_COLORS[type(tx)][1],
            'height': 0.5*self.opts.max_track_width,
            'zorder': 1.5
        }
        cds_kwargs = {
            'type': 'rect',
            'edgecolor': 'k',
            'facecolor': TRANSCRIPT_COLORS[type(tx)][0],
            'zorder': 1.5
        }
        if tx.orfs:
            orf = tx.primary_orf
            if orf.utr5:
                for exon in orf.utr5.exons:
                    if self.strand is Strand.PLUS:
                        start = exon.start
                        stop = min(exon.stop, orf.start)
                    elif self.strand is Strand.MINUS:
                        start = max(exon.start, orf.stop)
                        stop = exon.stop
                    self.draw_region(track, start=start, stop=stop, **utr_kwargs)
            for exon in orf.exons:
                start = max(exon.start, orf.start)
                stop = min(exon.stop, orf.stop)
                self.draw_region(track, start=start, stop=stop, **cds_kwargs)
            if orf.utr3:
                for exon in orf.utr3.exons:
                    if self.strand is Strand.PLUS:
                        start = max(exon.start, orf.stop)
                        stop = exon.stop
                    elif self.strand is Strand.MINUS:
                        start = exon.start
                        stop = min(exon.stop, orf.start)
                    self.draw_region(track, start=start, stop=stop, **utr_kwargs)
        else:
            for exon in tx.exons:
                self.draw_region(track, start=exon.start, stop=exon.stop, **utr_kwargs)

        for exon in tx.exons:
            # label every 5th exon in anchor isoform for easier navigation
            if track == 0 and exon.position % 5 == 0:
                self.draw_text((exon.start + exon.stop)//2, track - self.opts.max_track_width, f'E{exon.position}', ha='center', va='baseline')

        for orf in tx.orfs:
            first_res = orf.protein.residues[0]
            last_res = orf.protein.residues[-1]
            if first_res.amino_acid is AminoAcid.MET:
                start_codon = first_res.codon[0].coordinate
                self.draw_point(track, start_codon, color='lime')
            if last_res.amino_acid is AminoAcid.STOP:
                stop_codon = last_res.codon[2].coordinate
                self.draw_point(track, stop_codon, color='red')

        if hasattr(tx, 'start_nf') and tx.start_nf:
            self.draw_text(tx.start if self.strand is Strand.PLUS else tx.stop, track, '! ', ha='right', va='center', weight='bold', color='r')
        if hasattr(tx, 'end_nf') and tx.end_nf:
            self.draw_text(tx.stop if self.strand is Strand.PLUS else tx.start, track, ' !', ha='left', va='center', weight='bold', color='r')

    def draw_all_isoforms(self, subplot_spec = None):
        """Plot all isoforms."""
        R = len(self.transcripts)
        C = len(self._columns)
        self.fig = plt.figure()
        self._bax = BrokenAxes(fig=self.fig, xlims=self.xlims, ylims=((R-0.5, -0.5),), wspace=0, d=0.008, subplot_spec=subplot_spec)
        self._handles['Intron'] = mlines.Line2D([], [], linewidth=1.5, color='gray')
        self._handles['Exon (translated)'] = mpatches.Patch(facecolor=TRANSCRIPT_COLORS[None][0], edgecolor='k')
        self._handles['Exon (untranslated)'] = mpatches.Patch(facecolor=TRANSCRIPT_COLORS[None][1], edgecolor='k')
        self._handles['Start codon'] = mlines.Line2D([], [], linestyle='None', color='lime', marker='|', markersize=10, markeredgewidth=1)
        self._handles['Stop codon'] = mlines.Line2D([], [], linestyle='None', color='red', marker='|', markersize=10, markeredgewidth=1)


        for i, tx in enumerate(self.transcripts):
            with ExceptionLogger(f'Error plotting {tx}'):
                if tx:
                    self.draw_isoform(tx, i)

        # plot genomic region label
        # gene = self.transcripts[0].gene
        # start, end = self.xlims[0][0], self.xlims[-1][1]
        # self._bax.set_title(f'{gene.chromosome}({self.strand}):{start}-{end}')

        # hide y axis spine
        left_subaxes = self._bax.axs[0]
        left_subaxes.spines['left'].set_visible(False)
        left_subaxes.set_yticks([])

        # plot table
        # https://stackoverflow.com/a/57169705
        table = self._bax.big_ax.table(
            rowLabels = [getattr(tx, 'name', '') for tx in self.transcripts],
            colLabels = list(self._columns.keys()),
            cellText = [[f(tx) if tx else '' for f in self._columns.values()] for tx in self.transcripts],
            cellLoc = 'center',
            edges = 'open',
            bbox = (-0.1*C, 0.0, 0.1*C, (R+1)/R)
        )
        # table.auto_set_font_size(False)
        # table.set_fontsize(10)

        # rotate x axis tick labels for better readability
        for subaxes in self._bax.axs:
            subaxes.xaxis.set_major_formatter('{x:.0f}')
            for label in subaxes.get_xticklabels():
                label.set_va('top')
                label.set_rotation(90)
                label.set_size(8)

    def draw_frameshifts(self, anchor: Optional['Transcript'] = None, hatch_color='white'):
        """Plot relative frameshifts on all isoforms. Uses first isoform as the anchor by default."""
        self._handles['Frame +1'] = mpatches.Patch(facecolor='k', edgecolor='w', hatch=REL_FRAME_STYLE[CodonAlignmentCategory.FRAME_AHEAD])
        self._handles['Frame +2'] = mpatches.Patch(facecolor='k', edgecolor='w', hatch=REL_FRAME_STYLE[CodonAlignmentCategory.FRAME_BEHIND])

        if anchor is None:
            anchor = next(filter(None, self.transcripts))
        if not anchor or not anchor.protein:
            warn(
                'Cannot draw frameshifts without an anchor ORF'
            )
            return
        for i, other in enumerate(self.transcripts):
            if not other or not other.protein or other is anchor:
                continue
            aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
            for block in filter(lambda block: block.category in FRAMESHIFT, aln.blocks):
                for exons, residues in groupby(other.protein.residues[block.other_range.start:block.other_range.stop], key=attrgetter('exons')):
                    if len(exons) > 1:
                        continue
                    r1, r2 = tee(residues, 2)
                    start = first(r1).codon[1].coordinate
                    stop = last(r2).codon[1].coordinate
                    self.draw_region(
                        track = i,
                        start = start,
                        stop = stop,
                        facecolor = 'none',
                        edgecolor = hatch_color,
                        linewidth = 0.0,
                        zorder = 1.9,
                        hatch = REL_FRAME_STYLE[block.category]
                    )

    def draw_codon_alignment_blocks(self, cd_aln: 'CodonAlignment', alpha: float = 0.5):
        for category, color in CBLOCK_COLORS.items():
            label = category.name.capitalize().replace('_', ' ')
            if label not in self._handles:
                self._handles[label] = mpatches.Patch(facecolor=color)
        height = 0.25*self.opts.max_track_width
        track = self.transcripts.index(cd_aln.other.transcript)
        for block in filter(lambda block: block.category is not CodonAlignmentCategory.MATCH, cd_aln.blocks):
            if block.other_range:
                start = cd_aln.other.residues[block.other_range[0]].codon[1].coordinate
                stop = cd_aln.other.residues[block.other_range[-1]].codon[1].coordinate
            else:
                start = cd_aln.anchor.residues[block.anchor_range[0]].codon[1].coordinate
                stop = cd_aln.anchor.residues[block.anchor_range[-1]].codon[1].coordinate
            if block.category in SPLIT_CODON:
                self.draw_point(  # TODO: fix
                    track,
                    start,
                    height = height,
                    type = 'lollipop',
                    marker = '.',
                    color = CBLOCK_COLORS[block.category],
                    zorder = 1.9,
                    alpha = alpha
                )
            else:
                self.draw_region(
                    track,
                    start,
                    stop,
                    y_offset = -0.5*self.opts.max_track_width,
                    height = -height,
                    facecolor = CBLOCK_COLORS[block.category],
                    alpha = alpha
                )

    def draw_protein_alignment_blocks(self, pblocks: Iterable['ProteinAlignmentBlock'], anchor: 'Protein', other: 'Protein', alpha: float = 1.0):
        for category, color in PBLOCK_COLORS.items():
            label = category.name.capitalize().replace('_', ' ')
            if label not in self._handles:
                self._handles[label] = mpatches.Patch(facecolor=color)
        self._handles['Ragged 5\' end'] = mlines.Line2D([], [], linestyle='None', color='#999999', marker='<', markersize=8, markeredgewidth=1)
        self._handles['Ragged 3\' end'] = mlines.Line2D([], [], linestyle='None', color='#999999', marker='>', markersize=8, markeredgewidth=1)

        for pblock in filter(lambda block: block.category is not SequenceAlignmentCategory.MATCH, pblocks):
            anchor_start, anchor_stop, other_start, other_stop = None, None, None, None
            if pblock.category is not SequenceAlignmentCategory.INSERTION:
                anchor_start = anchor.transcript.get_genome_coord_from_transcript_coord(
                    anchor.get_transcript_coord_from_protein_coord(pblock.anchor_range[0]) + 1
                ).coordinate
                anchor_stop = anchor.transcript.get_genome_coord_from_transcript_coord(
                    anchor.get_transcript_coord_from_protein_coord(pblock.anchor_range[-1]) + 1
                ).coordinate
            if pblock.category is not SequenceAlignmentCategory.DELETION:
                other_start = other.transcript.get_genome_coord_from_transcript_coord(
                    other.get_transcript_coord_from_protein_coord(pblock.other_range[0]) + 1
                ).coordinate
                other_stop = other.transcript.get_genome_coord_from_transcript_coord(
                    other.get_transcript_coord_from_protein_coord(pblock.other_range[-1]) + 1
                ).coordinate

            other_track = self.transcripts.index(other.transcript)
            lollipop_direction = 1 if pblock.category is SequenceAlignmentCategory.INSERTION else -1

            if pblock.ragged5:
                self.draw_point(
                    other_track,
                    pos = anchor_start,
                    ylims = (lollipop_direction*0.75*self.opts.max_track_width, 0),
                    linewidth = 0,
                    marker = '<',
                    markersize = 6,
                    color = PBLOCK_COLORS[pblock.category],
                    zorder = 1.9
                )
            if pblock.ragged3:
                self.draw_point(
                    other_track,
                    pos = anchor_stop,
                    ylims = (lollipop_direction*0.75*self.opts.max_track_width, 0),
                    linewidth = 0,
                    marker = '>',
                    markersize = 6,
                    color = PBLOCK_COLORS[pblock.category],
                    zorder = 1.9
                )
            self.draw_region(
                other_track,
                start = anchor_start,
                stop = anchor_stop,
                y_offset = -1.0*self.opts.max_track_width,
                height = 0.5*self.opts.max_track_width,
                edgecolor = 'none',
                facecolor = PBLOCK_COLORS[pblock.category],
                alpha = alpha
            )
            self.draw_region(
                other_track,
                start = other_start,
                stop = other_stop,
                y_offset = 0.5*self.opts.max_track_width,
                height = 0.5*self.opts.max_track_width,
                edgecolor = 'none',
                facecolor = PBLOCK_COLORS[pblock.category],
                alpha = alpha
            )

    def draw_features(self):
        h = self.opts.max_track_width
        feature_names = sorted({feature.name for tx in filter(None, self.transcripts) if tx.protein for feature in tx.protein.features if feature.type is not FeatureType.IDR})
        cmap = sns.color_palette('pastel', len(feature_names))
        colors = dict(zip(feature_names, cmap))
        colors.update(FEATURE_COLORS)
        self._handles.update({name: mpatches.Patch(facecolor=color) for name, color in colors.items()})
        for track, tx in enumerate(self.transcripts):
            if not tx or not tx.protein:
                continue
            features = tx.protein.features
            if not features:
                continue
            subtracks, n_subtracks = generate_subtracks(
                ((feature.protein_start, feature.protein_stop) for feature in features),
                (feature.name for feature in features)
            )
            for feature in features:
                subtrack = subtracks[feature.name]
                color = colors[feature.name]
                if feature.reference:
                    subfeatures = groupby(feature.residues, key=lambda res: (False, res.primary_exon))
                    # n_subtracks_temp = n_subtracks
                else:
                    subfeatures = groupby(feature.residues, key=lambda res: (res in feature.altered_residues, res.primary_exon))
                    # n_subtracks_temp = 2*n_subtracks
                for (altered, _), subfeature in subfeatures:
                    subfeature = list(subfeature)
                    start = subfeature[0].codon[1].coordinate
                    stop = subfeature[-1].codon[1].coordinate
                    self.draw_region(
                        track,
                        start = start,
                        stop = stop,
                        y_offset = (-0.5 + subtrack/n_subtracks)*h,
                        height = h/n_subtracks,
                        edgecolor = 'none',
                        facecolor = color,
                        alpha = 0.5 if altered else 1.0,
                        zorder = 1.8,
                        label = feature.name
                    )
                # draw box behind entire feature
                self.draw_region(
                    track,
                    start = feature.residues[0].codon[1].coordinate,
                    stop = feature.residues[-1].codon[1].coordinate,
                    y_offset = (-0.5 + subtrack/n_subtracks)*h,
                    height = h/n_subtracks,
                    edgecolor = 'none',
                    facecolor = color,
                    alpha = 0.5,
                    zorder = 1.4
                )

    # --- NEW METHOD (Corrected) ---
    def draw_variants(self, variants: Iterable[Any], color='red', marker='v'):
        """
        Draws markers for genomic variants (SNPs) on the isoform tracks.
        Only plots a variant on a track if it falls within an exon of that isoform.
        """
        # Add entry to legend
        if 'Genetic Variant' not in self._handles:
            self._handles['Genetic Variant'] = mlines.Line2D(
                [], [], color=color, marker=marker, linestyle='None',
                markersize=8, label='Genetic Variant',
                markeredgecolor='black', markeredgewidth=0.5
            )

        for track_idx, tx in enumerate(self.transcripts):
            if not tx:
                continue

            for variant in variants:
                # Check if the variant position overlaps with this transcript (using simple exon bounds check)
                is_exonic = any(exon.start <= variant.position <= exon.stop for exon in tx.exons)

                if is_exonic:
                    # Draw the point on the specific track
                    self.draw_point(
                        track=track_idx,
                        pos=variant.position,
                        marker=marker,
                        color=color,
                        markersize=8,
                        zorder=3.0,
                        markeredgecolor='black',  # Fixed kwarg here
                        linewidth=0.5
                    )
