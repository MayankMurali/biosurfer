# Generic broken-axes drawing primitives, with no isoform/biosurfer-domain knowledge.
# IsoformPlot (plots/plotting.py) subclasses PlotCanvas and adds the domain-specific
# draw_* layers on top of these primitives.
from copy import copy
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Dict, Iterable, Optional, Tuple, Union
from warnings import warn

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from biosurfer.core.constants import Strand
from intervaltree import IntervalTree

if TYPE_CHECKING:
    from biosurfer.core.models.biomolecules import Transcript
    from brokenaxes import BrokenAxes
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

StartStop = Tuple[int, int]

TableColumn = Callable[['Transcript'], str]


@dataclass
class IsoformPlotOptions:
    """Bundles various options for adjusting plots made by IsoformPlot."""
    intron_spacing: int = 30  # number of bases to show in each intron
    track_spacing: float = 2  # ratio of space between tracks to max track width
    subtle_splicing_threshold: int = 20  # maximum difference (in bases) between exon boundaries to display subtle splicing

    @property
    def max_track_width(self) -> float:
        return 1/(self.track_spacing + 1)

    @max_track_width.setter
    def max_track_width(self, width: float):
        self.track_spacing = (1 - width)/width


class PlotCanvas:
    """Generic broken-axes drawing surface: tracks, subaxes, and low-level artist
    primitives (points, regions, text, legend). No knowledge of transcripts/isoforms --
    subclasses provide domain-specific draw_* methods built on top of these."""
    def __init__(self, strand: 'Strand', num_tracks: int, **kwargs):
        self.strand = strand
        self.num_tracks = num_tracks
        self.fig: Optional['Figure'] = None
        self._bax: Optional['BrokenAxes'] = None
        self.opts = IsoformPlotOptions(**kwargs)
        # keep track of artists for legend
        self._handles = dict()

    # Internally, PlotCanvas stores _subaxes, which maps each genomic region to the subaxes that plots the region's features.
    # The xlims property provides a simple interface to allow users to control which genomic regions are plotted.
    @property
    def xlims(self) -> Tuple[StartStop]:
        """Coordinates of the genomic regions to be plotted, as a tuple of (start, end) tuples."""
        return self._xlims

    @xlims.setter
    def xlims(self, xlims: Iterable[StartStop]):
        xregions = IntervalTree.from_tuples((min(start, stop), max(start, stop)+1) for start, stop in xlims)  # condense xlims into single IntervalTree object
        xregions.merge_equals()
        xregions.merge_overlaps()
        xregions.merge_neighbors()
        xregions = sorted(xregions.all_intervals)
        if self.strand is Strand.MINUS:
            xregions.reverse()
        self._subaxes = IntervalTree.from_tuples((start, end, i) for i, (start, end, _) in enumerate(xregions))
        self._xlims = tuple((start, end-1) if self.strand is Strand.PLUS else (end-1, start) for start, end , _ in xregions)

    # This method speeds up plotting by allowing PlotCanvas to add artists only to the subaxes where they are needed.
    def _get_subaxes(self, xcoords: Union[int, StartStop]) -> Tuple['Axes']:
        """For a specific coordinate or range of coordinates, retrieve corresponding subaxes."""
        if isinstance(xcoords, tuple):
            if xcoords[0] > xcoords[1]:
                xcoords = (xcoords[1], xcoords[0])
            xcoords = slice(*xcoords)
        subax_ids = [interval[-1] for interval in self._subaxes[xcoords]]
        if not subax_ids:
            raise ValueError(f"{xcoords} is not within plot's xlims")
        return tuple(self._bax.axs[id] for id in subax_ids)

    def draw_point(self, track: int, pos: int,
                    ylims: tuple[float, float] = None,
                    marker='', linewidth=1, **kwargs):
        """Draw a feature at a specific point. Appears as a vertical line with an optional marker."""
        if ylims is None:
            ylims = -0.5*self.opts.max_track_width, 0.5*self.opts.max_track_width
        artist = mlines.Line2D(
            xdata = (pos, pos),
            ydata = (track + ylims[0], track + ylims[1]),
            linewidth = linewidth,
            marker = marker,
            markevery = 2,
            **kwargs
        )

        try:
            subaxes = self._get_subaxes(pos)[0]
        except ValueError as e:
            warn(str(e))
        else:
            subaxes.add_artist(artist)
        return artist

    def draw_region(self, track: int, start: int, stop: int,
                    y_offset: Optional[float] = None,
                    height: Optional[float] = None,
                    type='rect', **kwargs):
        """Draw a feature that spans a region. Appearance types are rectangle and line."""
        if start == stop:
            return
        # TODO: make type an enum?
        if type == 'rect':
            if height is None:
                height = self.opts.max_track_width
            if y_offset is None:
                y_offset = -0.5*height
            artist = mpatches.Rectangle(
                xy = (start, track + y_offset),
                width = stop - start,
                height = height,
                **kwargs
            )
        elif type == 'line':
            if y_offset is None:
                y_offset = 0
            artist = mlines.Line2D(
                xdata = (start, stop),
                ydata = (track + y_offset, track + y_offset),
                **kwargs
            )
        else:
            raise ValueError(f'Region type "{type}" is not defined')

        subaxes = self._get_subaxes((start, stop))
        for ax in subaxes:
            ax.add_artist(copy(artist))
        return artist

    def draw_background_rect(self, start: int, stop: int,
                            track_first: int = None, track_last: int = None,
                            padding: float = None, **kwargs):
        """Draw a rectangle in the background of the plot."""
        if start == stop:
            return
        if track_first is None:
            track_first = 0
        if track_last is None:
            track_last = self.num_tracks - 1
        if padding is None:
            padding = self.opts.max_track_width
        top = track_first - padding
        bottom = track_last + padding
        artist = mpatches.Rectangle(
            xy = (start, top),
            width = stop - start,
            height = bottom - top,
            zorder = 0.5,
            **kwargs
        )

        subaxes = self._get_subaxes((start, stop))
        for ax in subaxes:
            ax.add_artist(copy(artist))
        return artist

    def draw_text(self, x: int, y: float, text: str, **kwargs):
        """Draw text at a specific location. x-coordinate is genomic, y-coordinate is w/ respect to tracks (0-indexed).
        Ex: x=20000, y=2 will center text on track 2 at position 20,000."""
        # TODO: make this use Axes.annotate instead
        # we can't know how much horizontal space text will take up ahead of time
        # so text is plotted using BrokenAxes' big_ax, since it spans the entire x-axis
        big_ax = self._bax.big_ax
        try:
            subaxes = self._get_subaxes(x)[0]  # grab coord transform from correct subaxes
        except ValueError as e:
            warn(str(e))
        else:
            big_ax.text(x, y, text, transform=subaxes.transData, **kwargs)

    def draw_legend(self, only_labels: Optional[Iterable[str]] = None, except_labels: Optional[Iterable[str]] = None, **kwargs):
        if only_labels and except_labels:
            raise ValueError('Cannot set both "only_labels" and "except_labels"')
        elif only_labels:
            labels = [label for label in self._handles if label in only_labels]
        elif except_labels:
            labels = [label for label in self._handles if label not in except_labels]
        else:
            labels = list(self._handles.keys())
        handles = [self._handles[label] for label in labels]
        self.fig.legend(
            handles = handles,
            labels = labels,
            # ncol = 1,
            # loc = 'center left',
            # mode = 'expand',
            # bbox_to_anchor = (1.05, 0.5),
            **kwargs
        )

    def savefig(self, fig_path):
        self.fig.set_size_inches(20, 0.8 + 0.4*self.num_tracks)
        plt.figure(self.fig)
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=700, bbox_inches='tight')
