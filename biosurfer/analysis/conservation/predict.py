"""
Core scoring logic for conservation (phyloP/phastCons) and AlphaMissense
pathogenicity, both queried over genomic ranges.

Coordinate conversion (protein range -> transcript range -> genomic ranges)
follows the same pattern already established in
``biosurfer.analysis.genetics_analyzer._query_gwas_in_protein_range`` --
duplicated locally rather than imported, so this module stays independently
testable/isolated per the project's Stage B rules (a new subpackage should
not require edits to, or introduce a dependency on the internals of,
another analysis module).
"""
from typing import List, Optional

import pysam

from biosurfer.core.models.biomolecules import Transcript

try:
    import pyBigWig
except ImportError:
    pyBigWig = None

# UCSC hg38, 100-way vertebrate alignment. bigWig supports random-access
# range queries over HTTP -- this does NOT download the (multi-GB)
# genome-wide file, only the byte ranges covering the queried region.
DEFAULT_PHYLOP_HG38_URL = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw'
DEFAULT_PHASTCONS_HG38_URL = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw'

# AlphaMissense_hg38.tsv.gz column layout (google-deepmind/alphamissense):
# CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
_AM_COLUMNS = ('chromosome', 'position', 'ref', 'alt', 'genome', 'uniprot_id',
               'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class')


def _require_pybigwig():
    if pyBigWig is None:
        raise ImportError(
            "Conservation scoring requires the optional 'pyBigWig' dependency, which "
            "isn't part of the base biosurfer install. Install it with:\n"
            "    pip install biosurfer[conservation]"
        )


def open_bigwig(path_or_url: str):
    """
    Open a local or remote (http/https/ftp) bigWig file. The caller is
    responsible for calling ``.close()`` on the result when done -- pass
    the same open handle into repeated ``score_genomic_range``/
    ``score_protein_range_conservation`` calls rather than reopening it per
    query (each open of a remote URL does a small handshake).
    """
    _require_pybigwig()
    return pyBigWig.open(path_or_url)


def score_genomic_range(bw, chromosome: str, start: int, stop: int) -> dict:
    """
    Mean/min/max conservation score over a 1-based, inclusive genomic range
    ``[start, stop]`` (biosurfer's usual coordinate convention). bigWig
    itself is 0-based half-open, so ``start`` is shifted here -- callers
    should always pass 1-based coordinates.

    ``bw`` just needs a pyBigWig-compatible ``.stats(chrom, start, end,
    type=...)`` method -- this function doesn't touch the ``pyBigWig``
    module itself (only ``open_bigwig`` does), so it works with a real
    handle or a duck-typed test double regardless of whether the optional
    ``pyBigWig`` dependency is actually installed.
    """
    if start > stop:
        start, stop = stop, start
    try:
        mean = bw.stats(chromosome, start - 1, stop, type='mean')[0]
        maximum = bw.stats(chromosome, start - 1, stop, type='max')[0]
        minimum = bw.stats(chromosome, start - 1, stop, type='min')[0]
    except RuntimeError as e:
        # pyBigWig raises RuntimeError for an unknown chromosome or an
        # out-of-bounds range rather than returning an empty result.
        return {'mean': None, 'max': None, 'min': None, 'error': str(e)}
    return {'mean': mean, 'max': maximum, 'min': minimum, 'error': None}


def _protein_range_to_genomic_ranges(transcript: Transcript, protein_range) -> List[tuple]:
    """
    Convert a protein coordinate range (e.g. a pblock/cblock's anchor_range
    or other_range) to a list of (genomic_start, genomic_stop) tuples, one
    per exon it spans -- mirrors
    ``genetics_analyzer._query_gwas_in_protein_range``'s conversion.
    """
    protein = transcript.protein
    if protein is None:
        return []
    start_aa = protein_range.start
    end_aa = protein_range.stop - 1  # inclusive index
    try:
        tx_start_nt = protein.get_transcript_coord_from_protein_coord(start_aa)
        tx_end_nt = protein.get_transcript_coord_from_protein_coord(end_aa) + 2
    except (ValueError, TypeError):
        return []

    genomic_ranges = []
    for exon in transcript.exons:
        ex_tx_start = exon.transcript_start - 1
        ex_tx_stop = exon.transcript_stop - 1
        overlap_start = max(tx_start_nt, ex_tx_start)
        overlap_end = min(tx_end_nt, ex_tx_stop)
        if overlap_start <= overlap_end:
            try:
                g_start = transcript.get_genome_coord_from_transcript_coord(overlap_start).coordinate
                g_end = transcript.get_genome_coord_from_transcript_coord(overlap_end).coordinate
            except Exception:
                continue
            if g_start > g_end:
                g_start, g_end = g_end, g_start
            genomic_ranges.append((g_start, g_end))
    return genomic_ranges


def score_protein_range_conservation(bw, transcript: Transcript, protein_range) -> dict:
    """
    Aggregate conservation score (mean of means, max of maxes, min of
    mins across exons) over the genomic span of a protein coordinate range.
    """
    genomic_ranges = _protein_range_to_genomic_ranges(transcript, protein_range)
    if not genomic_ranges:
        return {'mean': None, 'max': None, 'min': None, 'error': 'could not resolve genomic coordinates'}

    chromosome = transcript.gene.chromosome_id
    scores = [score_genomic_range(bw, chromosome, g_start, g_stop) for g_start, g_stop in genomic_ranges]
    means = [s['mean'] for s in scores if s['mean'] is not None]
    maxes = [s['max'] for s in scores if s['max'] is not None]
    mins = [s['min'] for s in scores if s['min'] is not None]
    return {
        'mean': (sum(means) / len(means)) if means else None,
        'max': max(maxes) if maxes else None,
        'min': min(mins) if mins else None,
        'error': None,
    }


def query_alphamissense(tabix_path: str, chromosome: str, start: int, stop: int) -> List[dict]:
    """
    Query AlphaMissense pathogenicity scores overlapping a 1-based
    inclusive genomic range, from a bgzipped + tabix-indexed AlphaMissense
    TSV (``AlphaMissense_hg38.tsv.gz`` from
    https://github.com/google-deepmind/alphamissense). AlphaMissense ships
    the TSV bgzipped but not tabix-indexed -- index it once with:

        tabix -s1 -b2 -e2 -c'#' AlphaMissense_hg38.tsv.gz

    Uses ``pysam`` (already a hard biosurfer dependency, used for VCF
    loading) -- no new dependency needed for this data source.
    """
    try:
        tabix = pysam.TabixFile(tabix_path)
    except (OSError, ValueError) as e:
        raise ValueError(
            f"Could not open {tabix_path!r} as a tabix-indexed file. AlphaMissense's "
            "downloaded TSV needs to be bgzipped and tabix-indexed first:\n"
            "    tabix -s1 -b2 -e2 -c'#' AlphaMissense_hg38.tsv.gz"
        ) from e

    rows = []
    try:
        for line in tabix.fetch(chromosome, max(0, start - 1), stop):
            fields = line.split('\t')
            if len(fields) < len(_AM_COLUMNS):
                continue
            row = dict(zip(_AM_COLUMNS, fields))
            row['position'] = int(row['position'])
            row['am_pathogenicity'] = float(row['am_pathogenicity'])
            rows.append(row)
    except ValueError:
        # chromosome not present in this tabix file's index at all
        return []
    finally:
        tabix.close()
    return rows


def query_alphamissense_for_protein_range(tabix_path: str, transcript: Transcript, protein_range) -> List[dict]:
    """Same protein-range -> genomic-range conversion as conservation scoring, applied to AlphaMissense."""
    genomic_ranges = _protein_range_to_genomic_ranges(transcript, protein_range)
    chromosome = transcript.gene.chromosome_id
    rows = []
    for g_start, g_stop in genomic_ranges:
        rows.extend(query_alphamissense(tabix_path, chromosome, g_start, g_stop))
    return rows
