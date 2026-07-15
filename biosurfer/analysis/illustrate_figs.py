# %%
from pathlib import Path
import colorsys
import matplotlib as mpl
import matplotlib.colors as mc
import matplotlib.font_manager as fm
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency
from itertools import combinations
from seaborn import color_palette
from re import M
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib.patches import Patch

def run_illustrate_analysis(pblock_table:Path, output: Path):
    """ Main plot function to invoke plotting for different pipelines/scripts.
    Args:
      pblock_table (Path): Path to the pblocks.tsv file.
      output (Path): Directory to save the illustrations.
    Returns:
      Nothing
    """
    
    #####################################################
    ## Setting configurations for illustrating figures ##
    #####################################################

    font = {
        'family': 'sans-serif',
        'sans-serif': ['Arial'],
        'weight': 'normal',
        'size': 16
    }
    mpl.rc('font', **font)

    # from https://stackoverflow.com/a/49601444
    def adjust_lightness(color, amount=0.5):
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        cnew = colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
        return mc.to_hex(cnew)

    PBLOCK_COLORS = {
        'DELETION': '#f800c0',
        'INSERTION': '#00c0f8',
        'SUBSTITUTION': '#f8c000',
    }

    PBLOCK_COLORS['SUBSTITUTION (reference)'] = adjust_lightness(PBLOCK_COLORS['SUBSTITUTION'], 1)
    PBLOCK_COLORS['SUBSTITUTION (alternative)'] = adjust_lightness(PBLOCK_COLORS['SUBSTITUTION'], 1)

    SECTION_COLORS = {
        'N-terminal': color_palette('pastel')[2],
        'Internal': color_palette('pastel')[7],
        'C-terminal': color_palette('pastel')[3],
        'Full-length': 'none',
    }

    NTERM_CLASSES = {
        'MUTUALLY_EXCLUSIVE': 'Mutually exclusive starts (MSX)',
        'DOWNSTREAM_SHARED': 'Shared downstream start (SDS)',
        'UPSTREAM_SHARED': 'Shared upstream start (SUS)',
        'MUTUALLY_SHARED': 'Mutually shared starts (MSS)'
    }
    NTERM_COLORS = dict(zip(
        NTERM_CLASSES.values(),
        color_palette('viridis_r', n_colors=len(NTERM_CLASSES)+1)[:-1]
    ))

    SPLICE_EVENT_COLORS = {
        'Intron': '#EBA85F',
        'Single exon': '#649FD2',
        'Alt. donor': '#86BB6F',
        'Alt. acceptor': '#A26FBB',
        'Compound': '#888888',
        'Frameshift': '#F7D76E',
    }

    CTERM_CLASSES = {
        'SPLICING': 'Splice-driven',
        'FRAMESHIFT': 'Frameshift-driven',
    }
    cterm_splice_palette = color_palette('RdPu_r', n_colors=6)
    cterm_frameshift_palette = color_palette('YlOrRd_r', n_colors=5)
    CTERM_PALETTE = [cterm_splice_palette[0], cterm_frameshift_palette[0]]

    pblocks = pd.read_csv(pblock_table, sep='\t')

    #######################################
    ## Genome-wide summary illustrations ##
    #######################################
    gw_output = output / 'gw_summary_plots'   
    gw_output.mkdir(exist_ok=True)
    # %% Fig2 panel A: Number of altered isoforms per gene vs number of genes
    fig = plt.figure(figsize=(4, 2.4))
    bins = list(range(1, 11)) + [100]
    ax = sns.histplot(
        x=pd.cut(
            pblocks.groupby('anchor')['other'].nunique(),
            bins=bins,
            right=False,
            labels=[str(x) for x in bins[:-2]] + [f'{bins[-2]}+'],
        ),
        shrink=0.75,
        color='#888888',
        edgecolor='k',
        alpha=1,
    )
    ax.set_xlabel('Number of alternative isoforms\nper gene')
    ax.set_ylabel('Number of genes')
    ax.set_ylim(0, 5000)
    ### output
    fig.savefig(gw_output / 'alternative-isoforms-per-gene.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    pblocks.groupby('anchor')['other'].nunique().to_frame(name='count').to_csv(gw_output / 'alternative-isoforms-per-gene-table.tsv', sep='\t')
    # %% Fig2 panel B: Number of observed pblocks per alternative protein isoforms
    fig = plt.figure(figsize=(4, 2.4))
    ax = sns.histplot(
        x=pd.cut(
            pblocks.groupby(['anchor', 'other']).size(),
            bins=[1, 2, 3, 4, 5, 14],
            right=False,
            labels=['1', '2', '3', '4', '5+']
        ),
        shrink=0.75,
        color='#888888',
        edgecolor='k',
        alpha=1,
    )
    ax.set_xlabel('Number of altered regions\nper isoform')
    ax.set_ylabel('Number of alternative\nprotein isoforms')

    ### output
    fig.savefig(gw_output / 'altered-regions-per-isoform.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    pblocks.groupby(['anchor', 'other']).size().to_frame(name='num_alt_regions').to_csv(gw_output / 'altered-regions-per-isoform-table.tsv', sep='\t')
    # %% Fig2 panel C: Distribution of lengths of the insertion, deletion and substitution affected regions for proteins 
    aa_loss = pblocks[pblocks['pblock_category'].isin({'DELETION', 'SUBSTITUTION'})].reset_index()[['anchor', 'other', 'pblock_category', 'aa_loss']]
    aa_loss['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (reference)', inplace=True)
    aa_loss.rename(columns={'aa_loss': 'length'}, inplace=True)
    aa_gain = pblocks[pblocks['pblock_category'].isin({'INSERTION', 'SUBSTITUTION'})].reset_index()[['anchor', 'other', 'pblock_category', 'aa_gain']]
    aa_gain['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (alternative)', inplace=True)
    aa_gain.rename(columns={'aa_gain': 'length'}, inplace=True)
    affected_lengths = pd.concat([aa_loss, aa_gain])

    binwidth = 50
    xmax = 600
    xtick = 200

    fig = plt.figure(figsize=(5, 2))
    data = affected_lengths[affected_lengths['pblock_category'] != 'SUBSTITUTION (alternative)']
    ax = sns.histplot(
        data=data,
        x='length',
        binwidth=binwidth,
        binrange=(0, xmax),
        stat='count',
        color='#808080',
        alpha=1,
    )
    ax.set_xlabel('Length of altered region (amino acids)')
    ax.set_ylabel('Number of\naltered regions')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.vlines(data['length'].median(), *ax.get_ylim(), color='#b0b0b0', linestyle='-', linewidth=1)

    ### output
    fig.savefig(gw_output / 'altered-region-affected-lengths.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    affected_lengths[affected_lengths['pblock_category'] != 'SUBSTITUTION (alternative)'].to_csv(gw_output / 'altered-region-affected-lengths-table.tsv', sep='\t')
    # %% Fig2 panel D: Distribution of the length of altered protein regions across the annotated proteome
    facets = sns.displot(
        data=affected_lengths,
        x='length',
        binwidth=binwidth,
        binrange=(0, xmax),
        stat='count',
        row='pblock_category',
        hue='pblock_category',
        palette=PBLOCK_COLORS,
        row_order=('DELETION', 'INSERTION', 'SUBSTITUTION (reference)', 'SUBSTITUTION (alternative)'),
        legend=False,
        alpha=1,
        height=2,
        aspect=2.5
    )
    facets.set_xlabels('Length of altered region (amino acids)')
    facets.set_ylabels('Number of\naltered regions')
    for category, ax in facets.axes_dict.items():
        ax.set_title(category.capitalize())
        ax.set_xticks(range(0, xmax + 1, xtick))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
        ax.vlines(affected_lengths[affected_lengths['pblock_category'] == category]['length'].median(), *ax.get_ylim(), color='#808080', linestyle='-', linewidth=1)

    ### output
    facets.fig.savefig(gw_output / 'altered-region-affected-lengths-categories.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    affected_lengths.to_csv(gw_output / 'altered-region-affected-lengths-categories-table.tsv', sep='\t')
    # %% Fig2 panel I =: Substitution scatter plot 
    plt.figure(figsize=(4.8, 3.6))
    ax = sns.histplot(
        data=pblocks[pblocks['pblock_category'] == 'SUBSTITUTION'],
        x='aa_gain',
        y='aa_loss',
        binwidth=binwidth / 2,
        stat='count',
        color=PBLOCK_COLORS['SUBSTITUTION'],
        legend=False,
        cbar=True,
        cbar_kws={
            'label': 'Number of regions',
        },
        alpha=1,
    )
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, xmax)
    ax.set_xticks(range(0, xmax + 1, xtick))
    ax.set_yticks(range(0, xmax + 1, xtick))
    ax.set_xlabel('Length of substitution region \nin alternative isoform (AA)')
    ax.set_ylabel('Length of substitution region \nin reference isoform (AA)')
    ### output
    plt.savefig(gw_output / 'substitution-reference-alternative-lengths.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    pblocks.query("pblock_category == 'SUBSTITUTION'")[['anchor', 'other', 'pblock_category', 'aa_gain', 'aa_loss']].to_csv(gw_output / 'substitution-reference-alternative-lengths-table.tsv', sep='\t')
    # %% Fig2 panel D: Pie chart
    category_counts = pblocks['pblock_category'].value_counts()
    total_pblocks = category_counts.sum()
    fig, ax = plt.subplots()
    wedges, texts, autotexts = plt.pie(
        category_counts,
        colors=category_counts.index.map(PBLOCK_COLORS),
        wedgeprops={'width': 0.4},
        startangle=180,
        counterclock=False,
        autopct=lambda x: f'{np.round(total_pblocks * x / 100):.0f}\n({x:.0f}%)',
        pctdistance=1.3,
    )
    for i, wedge in enumerate(wedges):
        wedge.set_edgecolor('k')
    ### output
    fig.savefig(gw_output / 'altered-region-category-donut.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    pblocks['pblock_category'].value_counts().to_csv(gw_output / 'altered-region-category-donut-table.tsv', sep='\t')
    # %%
    def get_section(nterm, cterm):
        if nterm and cterm:
            return 'Full-length'
        elif nterm:
            return 'N-terminal'
        elif cterm:
            return 'C-terminal'
        else:
            return 'Internal'

    pblocks['protein_section'] = list(map(get_section, ~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()))
    pblock_sections = pblocks['protein_section'].value_counts()

    fig, ax = plt.subplots(figsize=(6, 1))
    left = 0
    for section, color in SECTION_COLORS.items():
        val = pblock_sections[section]
        label = f'{val:g}\n({100 * val / pblock_sections.sum():0.1f}%)'
        if section == 'Full-length':
            left += 5000
            label_type = 'edge'
            padding = 5
        else:
            label_type = 'center'
            padding = 0
        bar = plt.barh(
            [0],
            val,
            left=left,
            color=color,
            edgecolor='k',
            label=section,
        )
        plt.bar_label(bar, labels=[label], label_type=label_type, padding=padding)
        left = left + pblock_sections[section]
    ax.legend(loc='upper left', bbox_to_anchor=(0, 0, 1, -0.1), ncols=2, frameon=False)
    plt.axis('off')
    ### output
    fig.savefig(gw_output / 'protein-section-counts.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    with open(gw_output / 'protein-section-counts-table.tsv', 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=SECTION_COLORS.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(SECTION_COLORS)
    # %%
    ##################################
    ## N-term summary illustrations ##
    ##################################
    nterm_output = output / 'nterm_summary_plots'   
    nterm_output.mkdir(exist_ok=True)
    nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != 'ALTERNATIVE_ORF') & (pblocks['cterm'].isna())].copy()
    nterm_pblocks['nterm'].replace(NTERM_CLASSES, inplace=True)
    nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: eval(x).intersection('BbPp')).astype(bool)

    # %% Fig3 panel A (both Alt TSS and 5' UTR AS)
    tss_fig = plt.figure(figsize=(5, 4))
    ax = sns.countplot(
        data=nterm_pblocks,
        y='nterm',
        order=NTERM_COLORS.keys(),
        palette=NTERM_COLORS,
        edgecolor='k',
        saturation=1,
    )
    sns.countplot(
        ax=ax,
        data=nterm_pblocks[nterm_pblocks['altTSS']],
        y='nterm',
        order=NTERM_COLORS.keys(),
        palette=NTERM_COLORS,
        edgecolor='k',
        fill=False,
        hatch='//',
    )
    ax.legend(
        loc=(0, 1),
        frameon=False,
        handles=[Patch(facecolor='w', edgecolor='k', hatch='///'), Patch(facecolor='w', edgecolor='k')],
        labels=['Alternative transcription start site', '5\' UTR alternative splicing'],
    )
    ax.set_xlabel('Number of alternative isoforms')
    ax.set_ylabel(None)
    plt.savefig(nterm_output / 'nterm-counts-all_mechanism.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    nterm_pblocks.query("nterm in ['Mutually exclusive starts (MSX)', 'Shared downstream start (SDS)']")[['anchor', 'other', 'nterm', 'altTSS']].to_csv(nterm_output / 'nterm-counts-all_mechanism.tsv', sep='\t')
    # %% Fig3 panel C: MXS vs SDS scatterplot 
    font = {
        'family': 'sans-serif',
        'sans-serif': ['Arial'],
        'weight': 'normal',
        'size': 10
    }
    mpl.rc('font', **font)

    # Filter the dataframe for 'Mutually exclusive starts (MXS)' and 'Shared downstream start (SDS)'
    msx_data = nterm_pblocks[nterm_pblocks['nterm'] == 'Mutually exclusive starts (MSX)']
    sds_data = nterm_pblocks[nterm_pblocks['nterm'] == 'Shared downstream start (SDS)']
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5.5, 5.5))
    msx_color = (0.565498, 0.84243, 0.262877)
    sds_color = (0.20803, 0.718701, 0.472873)

    sns.scatterplot(data=msx_data, x='aa_loss', y='aa_gain', marker='.', ax=axes[0], alpha=0.2,
                    color=msx_color)
    axes[0].set_title('Mutually exclusive starts (MSX)', fontsize=11)
    axes[0].set_xlabel('Reference \n(amino acids)', fontsize=10)
    axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=10)
    axes[0].set_xlim(0, 2000)
    axes[0].set_ylim(0, 2000)
    axes[0].set_aspect('equal')
    axes[0].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=sds_data, x='aa_loss', y='aa_gain', marker='.', ax=axes[1], alpha=0.2, color=sds_color)
    axes[1].set_title('Shared downstream start (SDS)', fontsize=11)
    axes[1].set_xlabel('Reference \n(amino acids)', fontsize=10)
    axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=10)
    axes[1].set_xlim(0, 2000)
    axes[1].set_ylim(0, 2000)
    axes[1].set_aspect('equal')
    axes[1].grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    # Save plot
    plt.savefig(nterm_output / 'nterm-rel-length-change_scatterplot.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    nterm_pblocks.query("nterm in ['Mutually exclusive starts (MSX)', 'Shared downstream start (SDS)']")[['anchor', 'other', 'aa_loss', 'aa_gain']].to_csv(nterm_output / 'nterm_mechanism_affected_len.tsv', sep='\t')
    # %%
    ############################################
    ## Internal region  summary illustrations ##
    ############################################
    internal_output = output / 'internal_summary_plots'   
    internal_output.mkdir(exist_ok=True)
    internal_pblocks = (
        pblocks[pblocks['internal']].
        drop(columns=[col for col in pblocks.columns if 'start' in col or 'stop' in col]).
        copy()
    )
    # convert string repr back to Python object
    internal_pblocks['tblock_events'] = internal_pblocks['tblock_events'].map(eval)
    internal_pblocks['events'] = internal_pblocks['events'].map(eval)
    internal_subcats = pd.DataFrame(
        {
            'Frameshift': internal_pblocks['frameshift'],
            'Intron': internal_pblocks['tblock_events'].isin({('I',), ('i',)}),
            'Alt. donor': internal_pblocks['tblock_events'].isin({('D',), ('d',)}),
            'Alt. acceptor': internal_pblocks['tblock_events'].isin({('A',), ('a',)}),
            'Single exon': internal_pblocks['tblock_events'].isin({('E',), ('e',)}),
            'Compound': [True for _ in internal_pblocks.index]
        }
    )
    subcat_order = ('Single exon', 'Alt. acceptor', 'Alt. donor', 'Intron', 'Compound', 'Frameshift')
    internal_pblocks['splice_event'] = internal_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(subcat_order, ordered=True))
    # %% Fig4 panel A: Internal splicing events frequencies
    internal_pblocks_fig = plt.figure(figsize=(4.6, 3.8))
    ax = sns.countplot(
        data=internal_pblocks.sort_values('pblock_category', ascending=True),
        y='splice_event',
        dodge=True,
        hue='pblock_category',
        palette=PBLOCK_COLORS,
        saturation=1,
        edgecolor='k',
    )
    plt.legend(loc='center right', labels=['Deletions', 'Insertions', 'Substitutions'])
    ax.set_xlabel('Number of altered internal regions')
    ax.set_ylabel(None)
    internal_pblocks_fig.savefig(internal_output / 'internal-events.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    internal_pblocks[['splice_event', 'pblock_category']].to_csv(internal_output / 'internal-events-table.tsv', sep='\t')
    # %% Fig4 panel C: Proportion of each internal protein region that are ragged codons
    internal_pblocks_ragged_fig = plt.figure(figsize=(4.6, 3.8))
    ax = sns.countplot(
        data=internal_pblocks.sort_values('pblock_category', ascending=True),
        y='splice_event',
        palette=SPLICE_EVENT_COLORS,
        saturation=1,
        edgecolor='k',
    )
    sns.countplot(
        ax=ax,
        data=internal_pblocks[internal_pblocks['split_codons']].sort_values('pblock_category', ascending=True),
        y='splice_event',
        fill=False,
        edgecolor='k',
        hatch='///',
    )
    plt.gca()
    ax.set_xlabel('Number of altered internal regions')
    ax.set_ylabel(None)
    internal_pblocks_ragged_fig.savefig(internal_output / 'internal-events-ragged.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    internal_pblocks[['splice_event', 'split_codons']].to_csv(internal_output / 'internal-events-ragged-table.tsv', sep='\t')
    
    alpha = 0.01
    ragged_contingency = pd.crosstab(internal_pblocks['split_codons'], internal_pblocks['splice_event'])
    chi2, p_all, dof, expected = chi2_contingency(ragged_contingency)

    ps = dict()
    for event1, event2 in combinations(internal_subcats.columns, 2):
        sub_contingency = ragged_contingency[[event1, event2]]
        _, ps[event1, event2], _, _ = chi2_contingency(sub_contingency)

    ps_sig = {k: p for k, p in ps.items() if p < alpha/len(ps)}
    ps_insig = {k: p for k, p in ps.items() if k not in ps_sig}

    # %%
    nagnag_pblocks = internal_pblocks[(internal_pblocks['splice_event'] == 'Alt. acceptor') & (internal_pblocks['length_change'].abs() == 1)]

    # %% Fig4 panel B: Frequency of compound splicing events
    internal_compound_pblocks = internal_pblocks[internal_pblocks['splice_event'] == 'Compound'].copy()

    internal_compound_subcats = pd.DataFrame(
        {
            'Multi-exon skipping': internal_compound_pblocks['events'] == frozenset('e'),
            'Exon skipping + \nalt. donor/acceptor': internal_compound_pblocks['events'].isin({
                frozenset(sorted('de')),
                frozenset(sorted('De')),
                frozenset(sorted('ea')),
                frozenset(sorted('eA')),
                frozenset(sorted('dea')),
                frozenset(sorted('Dea')),
                frozenset(sorted('deA')),
                frozenset(sorted('DeA')),
            }),
            'Mutually exclusive exons': internal_compound_pblocks['tblock_events'].isin({('E', 'e'), ('e', 'E')}),
            'Multi-exon inclusion': internal_compound_pblocks['events'] == frozenset('E'),
            'Alt. donor + alt. acceptor': internal_compound_pblocks['events'].isin({
                frozenset(sorted('ad')),
                frozenset(sorted('Ad')),
                frozenset(sorted('aD')),
                frozenset(sorted('AD')),
            }),
            'Exon inclusion + \nalt. donor/acceptor': internal_compound_pblocks['events'].isin({
                frozenset(sorted('dE')),
                frozenset(sorted('DE')),
                frozenset(sorted('Ea')),
                frozenset(sorted('EA')),
                frozenset(sorted('dEa')),
                frozenset(sorted('DEa')),
                frozenset(sorted('dEA')),
                frozenset(sorted('DEA')),
            }),
            'Other': [True for _ in internal_compound_pblocks.index]
        }
    )
    internal_compound_pblocks['compound_subcat'] = internal_compound_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(internal_compound_subcats.columns, ordered=True))

    internal_pblocks_compound_fig = plt.figure(figsize=(3, 3))
    ax = sns.countplot(
            data=internal_compound_pblocks,
            y='compound_subcat',
            palette='Greys_r',
            saturation=1,
            edgecolor='k',
    )
    ax.set_xlabel('Number of altered\ninternal regions'),
    ax.set_ylabel(None)
    internal_pblocks_compound_fig.savefig(internal_output / 'internal-compound-events.png', dpi=500, facecolor=None, bbox_inches='tight')
    # Output source data
    internal_compound_pblocks[['anchor', 'other', 'compound_subcat']].to_csv(internal_output / 'internal-compound-events-table.tsv', sep='\t')
    #%%
    ##################################
    ## C-term summary illustrations ##
    ##################################
    cterm_output = output / 'cterm_summary_plots'   
    cterm_output.mkdir(exist_ok=True)
    
    cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'].isna()) & (pblocks['cterm'] != "ALTERNATIVE_ORF") & (pblocks['cterm'] != "UNKNOWN")].copy()
    cterm_pblocks['cterm'] = cterm_pblocks['cterm'].map(CTERM_CLASSES).astype('category')
    # Changed string to set for intersection
    cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: set(x).intersection('BbPp')).astype(bool)
    
    #%% Fig5 panel A: Frequency of splice-driven and frameshift-driven C-terminal events
    cterm_fig = plt.figure(figsize=(3.8, 2))
    ax = sns.countplot(
        data = cterm_pblocks,
        y = 'cterm',
        order = CTERM_CLASSES.values(),
        palette = CTERM_PALETTE,
        saturation = 1,
        linewidth = 1,
        edgecolor = 'k',
    )
    ax.set_xlabel('Number of alternative isoforms')
    ax.set_ylabel('')
    plt.savefig(cterm_output/'cterm-class-counts.png', dpi=500, facecolor=None, bbox_inches='tight')
    #Output source data
    cterm_pblocks.query("cterm in ['Splice-driven', 'Frameshift-driven']")[['anchor','other','cterm']].to_csv(cterm_output/'cterm-class-counts-table.tsv', sep='\t')
    # %% Fig5 panel B: Frequency of splice-driven patterns
    cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
    single_ATE = (cterm_pblocks['cterm'] == 'Splice-driven') & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})
    cterm_splice_subcats = pd.DataFrame(
        {
            'Exon extension introduces termination': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
            'Alternative terminal exon(s)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
            'Poison exon inclusion': cterm_pblocks['up_stop_events'] == 'E',
            'Other': [True for _ in cterm_pblocks.index]
            #'Alternative last exon in UTR': cterm_pblocks['cblocks'].apply(lambda x: 'TRANSLATED' in x and 'DELETION' in x and 'UNTRANSLATED' not in x)
        }
    )
    cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

    cterm_splice_palette_dict = dict(zip(
        cterm_splice_subcats.columns,
        cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
    ))
    splice_subcat_order = tuple(cterm_splice_subcats.keys())

    cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
    single_ATE = (cterm_pblocks['cterm'] == 'Splice-driven') & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})

    cterm_splice_subcats = pd.DataFrame(
        {
            'Exon extension introduces \n termination (EXIT)': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
            'Alternative terminal \n exon(s) (ATE)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
            'Alternative last exon \n in UTR (ALE in UTR)': cterm_pblocks.apply(lambda row: 'TRANSLATED' in row['cblocks'] and 'DELETION' in row['cblocks'] and 'UNTRANSLATED' not in row['cblocks'] if row['cterm'] == 'Splice-driven' and row['splice_subcat'] == 'Other' else False, axis=1),
            'Poison exon inclusion': cterm_pblocks['up_stop_events'] == 'E',
            'Cut-out splice terminal \n exon (COSTE)': cterm_pblocks.apply(lambda row: 'DELETION' in row['cblocks'] and 'INSERTION' in row['cblocks'] and 'TRANSLATED' not in row['cblocks'] and 'UNTRANSLATED' not in row['cblocks'] and 'FRAME' not in row['cblocks'] and 'p' in row['tblock_events'] and row['tblock_events'].count('B') == 1 if row['cterm'] == 'Splice-driven' and row['splice_subcat'] == 'Other' else False, axis=1),
            'Other': [True for _ in cterm_pblocks.index]
        }
    )
    cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

    cterm_splice_palette_dict = dict(zip(
        cterm_splice_subcats.columns,
        cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + cterm_splice_palette[3:4] + cterm_splice_palette[4:5] +  ['#bbbbbb']
    ))
    splice_subcat_order = tuple(cterm_splice_subcats.keys())

    cterm_splice_fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    sns.countplot(
        ax = axs[0],
        data = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven'],
        y = 'splice_subcat',
        order = splice_subcat_order,
        palette = cterm_splice_palette_dict,
        saturation = 1,
        edgecolor = 'k',
        linewidth = 1,
    )
    axs[0].set_xlabel('Number of alternative isoforms')
    axs[0].set_ylabel(None)

    plt.savefig(cterm_output/'cterm-splicing-subcats.png', dpi=500, facecolor=None, bbox_inches='tight')
    #Output source data
    cterm_pblocks.assign(anchor_relative_length_change = cterm_pblocks['anchor_relative_length_change'].abs())[['anchor','other', 'splice_subcat','anchor_relative_length_change']].to_csv(cterm_output/'cterm-splicing-subcats-table.tsv', sep='\t')
    cterm_pblocks.to_csv(cterm_output / 'cterm_pblocks.tsv', sep='\t', index=False)
    
    # %% Alternative Last Exon in 3' UTR case from Splice-driven 'Other' category.
    cterm_pblock_splice = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven']
    cterm_splice_other = cterm_pblock_splice[cterm_pblock_splice['splice_subcat'] == 'Other']
    condition1 = cterm_splice_other['cblocks'].apply(lambda x: 'DELETION' in x and 'TRANSLATED' in x)
    condition2 = cterm_splice_other['cblocks'].apply(lambda x: 'UNTRANSLATED' not in x)
    cterm_aleutr = cterm_splice_other[condition1 & condition2].copy()
    cterm_aleutr.to_csv(cterm_output / 'cterm-splice-driven-ALEinUTR.tsv', sep='\t')

    # %% Cut-out splice terminal exon case from Splice-driven 'Other' category.
    condition3 = cterm_splice_other['cblocks'].apply(lambda x: 'DELETION' in x and 'INSERTION' in x)
    condition4 = cterm_splice_other['cblocks'].apply(lambda x: 'TRANSLATED' not in x and 'UNTRANSLATED' not in x and 'FRAME' not in x)
    condition5 = cterm_splice_other['tblock_events'].apply(lambda x: x.count('B') == 1 and 'p' in x)
    cterm_other_new = cterm_splice_other[condition3 & condition4 & condition5].copy()
    cterm_other_new.to_csv(cterm_output / 'cterm-splice-driven-other-NEW.tsv', sep='\t')
    # %% Fig5 panel C & D: 2D scatter plot v2 splice-driven vs frameshift-driven
    font = {
        'family': 'sans-serif',
        'sans-serif': ['Arial'],
        'weight': 'normal',
        'size': 10
    }
    mpl.rc('font', **font)
    msx_data = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven']
    sds_data = cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven']
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6, 6))
    msx_color = (0.5048212226066897, 0.00392156862745098, 0.47021914648212226)
    sds_color = (0.7885121107266436, 0.03238754325259515, 0.13656286043829297)

    sns.scatterplot(data=msx_data, x='aa_loss', y='aa_gain', marker='o', ax=axes[0], alpha=0.2,
                    color=msx_color)
    axes[0].set_title('Splice-driven', fontsize=13)
    axes[0].set_xlabel('Reference \n(amino acids)', fontsize=12)
    axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=12)
    axes[0].set_xlim(0, 3000)
    axes[0].set_ylim(0, 3000)
    axes[0].set_aspect('equal')
    axes[0].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=sds_data, x='aa_loss', y='aa_gain', marker='o', ax=axes[1], alpha=0.2, color=sds_color)
    axes[1].set_title('Frameshift-driven', fontsize=13)
    axes[1].set_xlabel('Reference \n(amino acids)', fontsize=12)
    axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=12)
    axes[1].set_xlim(0, 3000)
    axes[1].set_ylim(0, 3000)
    axes[1].set_aspect('equal')
    axes[1].grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(cterm_output / 'cterm-rel-length-change_scatterplot.png', dpi=800, facecolor=None, bbox_inches='tight')
    # Output source data
    cterm_pblocks.query("cterm in ['Splice-driven', 'Frameshift-driven']")[['anchor', 'other', 'aa_loss', 'aa_gain']].to_csv(cterm_output / 'cterm_mechanism_affected_len.tsv', sep='\t')
    
    
    # %% Supplementary Figure S5: 2D scatter plot v2 frameshift-driven subcats
    d1 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Exon extension introduces \n termination (EXIT)']
    d2 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Alternative terminal \n exon(s) (ATE)']
    d3 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Alternative last exon \n in UTR (ALE in UTR)']
    d4 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Poison exon inclusion']
    d5 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Cut-out splice terminal \n exon (COSTE)']

    fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(6, 15))
    colors = [(0.5048212226066897, 0.00392156862745098, 0.47021914648212226),
              (0.735840061514802, 0.061960784313725495, 0.5225682429834679), 
              (0.9094502114571319, 0.2894886582083814, 0.6086120722798923), 
              (0.9754555940023067, 0.5330257593233372, 0.6768935024990388), 
              (0.9859592464436755, 0.7293041138023837, 0.7404229142637447)]

    sns.scatterplot(data=d1, x='aa_loss', y='aa_gain', marker='o', ax=axes[0], alpha=0.2,
                    color=colors[0])
    axes[0].set_title('Exon extension introduces termination', fontsize=30, pad=20)
    axes[0].set_xlabel('Reference \n(amino acids)', fontsize=25)
    axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=25)
    axes[0].set_xlim(0, 3000)
    axes[0].set_ylim(0, 3000)
    axes[0].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=d2, x='aa_loss', y='aa_gain', marker='o', ax=axes[1], alpha=0.2,
                    color=colors[1])
    axes[1].set_title('Alternative terminal exon(s)', fontsize=30, pad=20)
    axes[1].set_xlabel('Reference \n(amino acids)', fontsize=25)
    axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=25)
    axes[1].set_xlim(0, 3000)
    axes[1].set_ylim(0, 3000)
    axes[1].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=d3, x='aa_loss', y='aa_gain', marker='o', ax=axes[2], alpha=0.2,
                    color=colors[2])
    axes[2].set_title('Alternative last exon in UTR', fontsize=30, pad=20)
    axes[2].set_xlabel('Reference \n(amino acids)', fontsize=25)
    axes[2].set_ylabel('Alternative \n(amino acids)', fontsize=25)
    axes[2].set_xlim(0, 3000)
    axes[2].set_ylim(0, 3000)
    axes[2].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=d4, x='aa_loss', y='aa_gain', marker='o', ax=axes[3], alpha=0.2,
                    color=colors[3])
    axes[3].set_title('Poison exon inclusion', fontsize=30, pad=20)
    axes[3].set_xlabel('Reference \n(amino acids)', fontsize=25)
    axes[3].set_ylabel('Alternative \n(amino acids)', fontsize=25)
    axes[3].set_xlim(0, 3000)
    axes[3].set_ylim(0, 3000)
    axes[3].grid(True, linestyle='--', linewidth=0.5)

    sns.scatterplot(data=d5, x='aa_loss', y='aa_gain', marker='o', ax=axes[4], alpha=0.2,
                    color=colors[4])
    axes[4].set_title('Cut-out splice terminal exon', fontsize=30, pad=20)
    axes[4].set_xlabel('Reference \n(amino acids)', fontsize=25)
    axes[4].set_ylabel('Alternative \n(amino acids)', fontsize=25)
    axes[4].set_xlim(0, 3000)
    axes[4].set_ylim(0, 3000)
    axes[4].grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(cterm_output / 'cterm-rel-splice-driven-subcat-length-change_scatterplot.png', dpi=500, facecolor=None, bbox_inches='tight')

if __name__ == "__main__":
    run_illustrate_analysis(pblock_table, output)