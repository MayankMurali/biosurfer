# %%
from pathlib import Path
from more_itertools import partition
from biosurfer.core.alignments import ProteinAlignment
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.helpers import (get_ids_from_gencode_fasta,
                                    get_ids_from_lrp_fasta,
                                    get_ids_from_pacbio_fasta, skip_gencode,
                                    skip_par_y)
from biosurfer.core.models.biomolecules import Gene, Transcript, Chromosome
from biosurfer.plots.plotting import IsoformPlot

from copy import copy
from dataclasses import dataclass
from itertools import chain, groupby, islice, tee
from operator import attrgetter, sub
from typing import (TYPE_CHECKING, Any, Callable, Collection, Dict, Iterable,
                    List, Literal, Optional, Set, Tuple, Union)
from warnings import filterwarnings, warn

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import Align
from biosurfer.core.alignments import CodonAlignment, ProjectedFeature
from biosurfer.core.constants import (FRAMESHIFT, SPLIT_CODON, AminoAcid,
                                      CodonAlignmentCategory, FeatureType,
                                      SequenceAlignmentCategory, Strand)
from biosurfer.core.helpers import (ExceptionLogger, Interval, IntervalTree,
                                    get_interval_overlap_graph)
from biosurfer.core.models.biomolecules import (GencodeTranscript,
                                                PacBioTranscript, Transcript)
from biosurfer.core.splice_events import (AcceptorSpliceEvent,
                                          DonorSpliceEvent, ExonBypassEvent,
                                          ExonSpliceEvent, IntronSpliceEvent)
from brokenaxes import BrokenAxes
from graph_tool import Graph
from graph_tool.topology import sequential_vertex_coloring
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from more_itertools import first, last, only
import pandas as pd
import csv
if TYPE_CHECKING:
    from biosurfer.core.alignments import ProteinAlignmentBlock, CodonAlignmentBlock
    from biosurfer.core.models.biomolecules import Protein
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

def get_frameshift_positions(s, chr, writer):
    """ 
    Function to frameshift postions 
    """
    chr_obj = Chromosome.from_name(s, chr)
    gene_list = list(chr_obj.genes)
    i = iter(gene_list)
    for gene in gene_list:
        # print("\n Gene :", gene)
        if gene:
            try:
                gene_obj = Gene.from_name(s, str(gene))
            except:
                continue
            if gene_obj is None:
                #print(f'Gene "{gene}" not found in database')
                transcripts = dict()
                anchor = None
            else:
                transcripts = {tx.accession: tx for tx in gene_obj.transcripts}
                try:
                    anchor = max(transcripts.values(), key=lambda tx: getattr(tx, 'appris', APPRIS.NONE))
                except:
                    anchor=None
                    print("Error with :", gene)
            others = [tx for tx in transcripts.values() if tx is not anchor]
        else:
            transcripts: dict[str, Transcript] = {tx.accession: tx for tx in Transcript.from_accessions(s, transcript_ids).values()}
            not_found, found = partition(lambda tx_id: tx_id in transcripts, transcript_ids)
            for tx_id in not_found:
                print(f'Transcript ID "{tx_id}" not found in database')
            if transcript_ids:
                anchor = transcripts.get(transcript_ids[0], None)
            else:
                #print('No isoforms provided')
                anchor = None
            others = [tx for tx in map(transcripts.get, found) if tx is not anchor]

        if anchor:
            # print(f'Reference isoform: {anchor}')
            
            gene = anchor.gene.name

            alns: dict[Transcript, ProteinAlignment] = dict()
            for other in others:
                if anchor.protein is None or other.protein is None:
                    alns[other] = None
                else:
                    try:
                        alns[other] = ProteinAlignment.from_proteins(anchor.protein, other.protein)
                    except:
                        # print(f'Could not plot isoform {other}')
                        continue
            

            #plot = IsoformPlot([anchor] + list(alns.keys()))
            transcripts: Iterable['Transcript']= [anchor] + list(alns.keys())
            transcripts: List['Transcript'] = list(transcripts)  # list of orf objects to be drawn
            gene = {tx.gene for tx in filter(None, transcripts)}
            if len(gene) > 1:
                raise ValueError(f'Found isoforms from multiple genes: {", ".join(g.name for g in gene)}')
            strand = only(
                {tx.strand for tx in filter(None, transcripts)},
                too_long = ValueError("Can't plot isoforms from different strands")
            )
            strand: Strand = strand

            # plot.draw_all_isoforms()

            for i, tx in enumerate(transcripts):
                with ExceptionLogger(f'Error plotting {tx}'):
                    if tx:
                        #self.draw_isoform(tx, i)
                        
                        start, stop = tx.start, tx.stop
                        align_start, align_stop = 'right', 'left'
                        if strand is Strand.MINUS:
                            align_start, align_stop = align_stop, align_start
                        transcript_start = start
                        transcript_stop = stop
                        isoform_name = tx
                        isoform_accession = tx.accession
                        transcript_orf_detail = tx.orfs
                        if tx.orfs:
                            orf = tx.primary_orf
                            if orf.utr5:
                                for exon in orf.utr5.exons:
                                    if strand is Strand.PLUS:
                                        start = exon.start
                                        stop = min(exon.stop, orf.start)
                                    elif strand is Strand.MINUS:
                                        start = max(exon.start, orf.stop)
                                        stop = exon.stop
                                UTR5_start = start
                                UTR5_stop = stop
                            for exon in orf.exons:
                                start = max(exon.start, orf.start)
                                stop = min(exon.stop, orf.stop)
                                exon_start = start
                                exon_stop = stop
                                
                            if orf.utr3:
                                for exon in orf.utr3.exons:
                                    if strand is Strand.PLUS:
                                        start = max(exon.start, orf.stop)
                                        stop = exon.stop
                                    elif strand is Strand.MINUS:
                                        start = exon.start
                                        stop = min(exon.stop, orf.start)
                                UTR3_start = start
                                UTR3_stop = stop
                        else:
                            for exon in tx.exons:
                                exon_start = start
                                exon_stop = stop
                        
                        for orf in tx.orfs:
                            first_res = orf.protein.residues[0]
                            last_res = orf.protein.residues[-1]
                            if first_res.amino_acid is AminoAcid.MET:
                                try:
                                    start_codon = first_res.codon[0].coordinate
                                except ValueError:
                                    print ("AttributeError avoided")
                            if last_res.amino_acid is AminoAcid.STOP:
                                stop_codon = last_res.codon[2].coordinate

            # plot.draw_frameshifts()  
            if anchor is None:
                anchor = next(filter(None, transcripts))
            if not anchor or not anchor.protein:
                warn(
                    'Cannot draw frameshifts without an anchor ORF'
                )
            for i, other in enumerate(transcripts):
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
                        frameshift_start = stop
                        frameshift_stop = start    
                        row = (gene,other,other.accession,frameshift_start,frameshift_stop)
                        writer.writerow(row)                        
                        

if __name__ == "__main__":
    db = Database('gencode_v44_comp')
    s = db.get_session() 

    chrs = [f'chr{i}' for i in list(range(23, 23)) + ['X']]    
    for chr in chrs:
        with open(f'{chr}_frameshift_positions.csv', 'wt') as f:
            writer = csv.writer(f)
            writer.writerow(('gene','isoform_name','accession','frameshift_start','frameshift_stop'))
            get_frameshift_positions(s, chr, writer)