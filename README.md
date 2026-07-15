
# Biosurfer

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7297008.svg)](https://doi.org/10.5281/zenodo.7297008)


"Surf" the biological network, from genome to transcriptome to proteome and back to gain insights into human disease biology.

> This is a fork ([MayankMurali/biosurfer](https://github.com/MayankMurali/biosurfer))
> undergoing internal modernization with no change to the CLI, database schema, or
> install workflow below. See [`CHANGELOG.md`](CHANGELOG.md) for what's changed and why.

**Contents**

- [Installation](#installation)
- [Usage](#usage)
  - [1. Load database](#1-load-database)
  - [2. Hybrid alignment](#2-hybrid-alignment)
  - [3. Visualize protein isoforms](#3-visualize-protein-isoforms)
  - [4. Genetics & Risk Analysis](#4-genetics--risk-analysis)
  - [5. Nonsense-Mediated Decay (NMD) Prediction](#5-nonsense-mediated-decay-nmd-prediction)
  - [6. Conservation & Missense Pathogenicity Scoring](#6-conservation--missense-pathogenicity-scoring)
  - [7. Splice-Site-Strength Scoring](#7-splice-site-strength-scoring)
  - [8. Intrinsic Disorder Scoring](#8-intrinsic-disorder-scoring)
  - [9. Interactive Isoform Viewer](#9-interactive-isoform-viewer)
- [References](#references)

## Installation
 

#### Building Requirements

* Python 3.10 or higher 
* Python packages (numpy, more-itertools, intervaltree, biopython, attrs, tqdm, pandas, pysam)
* Database (sqlalchemy >=1.4)
* Vizualization (matplotlib, brokenaxes, plotly)

Optional extras (base install works fully without any of these -- only needed for the
corresponding analysis in [Usage](#usage) below):

* `pip install biosurfer[conservation]` -- phyloP/phastCons conservation scoring (`pyBigWig`)
* `pip install biosurfer[splicing]` -- SpliceAI splice-site-strength scoring (`spliceai`, `tensorflow`)
* `pip install biosurfer[disorder]` -- metapredict intrinsic disorder scoring (`metapredict`)

#### Local building (without installation)


Clone the project repository (using SSH if need be) and create a [new conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) if needed. 

```
# Clone the repository
git clone https://github.com/sheynkman-lab/biosurfer
    
# Move to the folder
cd biosurfer
    
# Run setup 
pip install --editable .
``` 

## Usage

### Biosurfer command line options:
```
Usage: biosurfer [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...

Options:
  --help  Show this message and exit.

Commands:
  hybrid_alignment  This script runs hybrid alignment on the provided...
  load_db           Loads transcript and protein isoform information from...
  plot              Plot isoforms from a single gene, specified by...
```
* Download the [toy gencode data](https://zenodo.org/record/7297008/files/biosurfer_gencode_toy_data.zip?download=1) from Zenodo into the project directory.

### 1. Load database:

```
Usage: biosurfer load_db [OPTIONS]

  Loads transcript and protein isoform information from provided files into a
  Biosurfer database. A new database is created if the target database does
  not exist.

Options:
  -v, --verbose              Will print verbose messages
  -d, --db_name TEXT         Database name  [required]
  --source [GENCODE|PacBio]  Source of input data  [required]
  --gtf PATH                 Path to gtf file  [required]
  --tx_fasta PATH            Path to transcript sequence fasta file
                             [required]
  --tl_fasta PATH            Path to protein sequence fasta file  [required]
  --sqanti PATH              Path to SQANTI classification tsv file (only for
                             PacBio isoforms)
  --help                     Show this message and exit.
```

 #### Load database using GENCODE reference (toy version) 
 
```bash
biosurfer load_db --source=GENCODE --gtf biosurfer_gencode_toy_data/gencode.v38.toy.gtf --tx_fasta biosurfer_gencode_toy_data/gencode.v38.toy.transcripts.fa --tl_fasta biosurfer_gencode_toy_data/gencode.v38.toy.translations.fa --db_name gencode_toy
``` 
Running GENCODE files without ```--ref``` will 
 #### Load database using PacBio data without reference (WTC11 data)
 
```bash
biosurfer load_db --source=PacBio --gtf biosurfer_wtc11_data/wtc11_with_cds.gtf --tx_fasta biosurfer_wtc11_data/wtc11_corrected.fasta  --tl_fasta biosurfer_wtc11_data/wtc11_orf_refined.fasta --sqanti biosurfer_wtc11_data/wtc11_classification.txt --db_name wtc11_db
``` 

### 2. Hybrid alignment
* Run hybdrid alignment script on the created database. Create a directory to store the output files.

```shell
biosurfer hybrid_alignment -d gencode_toy -o output/gencode_toy -- gencode
```

```
Usage: biosurfer hybrid_alignment [OPTIONS]

  This script runs hybrid alignment on the provided database.

Options:
  -v, --verbose           Print verbose messages
  -d, --db_name TEXT      Database name  [required]
  -o, --output DIRECTORY  Directory for output files
  --gencode               Also compare all GENCODE isoforms of a gene against
                          its anchor isoform
  --anchors FILE          TSV file with gene names in column 1 and anchor
                          isoform IDs in column 2
  --help                  Show this message and exit.
```
> Please note that in the code, the terms *`anchor`* and *`other`* correspond to the *`reference`* and *`alternative`* isoforms mentioned in the manuscript.

### 3. Visualize protein isoforms
* To visualization isoforms of *CRYBG2* gene, run the following snippet.

```shell
biosurfer plot -d gencode_toy --gene CRYBG2
```

```
Usage: biosurfer plot [OPTIONS] [TRANSCRIPT_IDS]...

  Plot isoforms from a single gene, specified by TRANSCRIPT_IDS.

Options:
  -v, --verbose           Print verbose messages
  -o, --output DIRECTORY  Directory in which to save plots
  -d, --db_name TEXT      Database name  [required]
  --gene TEXT             Name of gene for which to plot all isoforms;
                          overrides TRANSCRIPT_IDS
  --help                  Show this message and exit.
```


### 4. Genetics & Risk Analysis
Analyze specific genes to find GWAS hits located in unique N-terminal regions (e.g., PPARG).

This module requires a VCF file (bgzipped & tabix indexed) and a GWAS summary statistics file (TSV format).

```shell
biosurfer analyze_nterm \
    --db_name gencode_v42 \
    --gene PPARG \
    --vcf path/to/genotypes.vcf.gz \
    --gwas path/to/gwas_summary.tsv \
    --output results_folder \
    --verbose
```

```
Usage: biosurfer analyze_nterm [OPTIONS]

  Analyzes N-terminal differences for a specific gene to identify 
  GWAS hits located in unique N-terminal regions.

Options:
  -v, --verbose            Print verbose messages
  -d, --db_name TEXT       Database name  [required]
  --gene TEXT              Target gene to analyze (e.g., PPARG) [required]
  --vcf TEXT               Path to VCF file (bgzipped & tabix indexed) [required]
  --gwas TEXT              Path to GWAS summary statistics (TSV) [required]
  -o, --output DIRECTORY   Directory for output tables
  --help                   Show this message and exit.
```

### 5. Nonsense-Mediated Decay (NMD) Prediction
Predicts NMD susceptibility (the standard 50-nt rule: a stop codon ≥50nt upstream of the
last exon-exon junction) per transcript, per isoform-vs-anchor delta, or per genetic
variant.

```shell
biosurfer analyze_nmd -d gencode_toy --gene CRYBG2 --mode isoform -o results_folder
```

```
Usage: biosurfer analyze_nmd [OPTIONS]

  Predicts nonsense-mediated decay (NMD) susceptibility for a gene's
  transcripts, using the standard 50-nt rule (already implemented as
  ORF.nmd in biosurfer.core.models.biomolecules).

Options:
  -v, --verbose           Print verbose messages
  -d, --db_name TEXT       Database name  [required]
  --gene TEXT              Target gene to analyze (e.g. PPARG)  [required]
  --mode [transcript|isoform|variant]
                           'transcript': per-transcript NMD status for every
                           transcript of the gene. 'isoform': NMD status of
                           every alternative isoform vs. the gene's anchor
                           isoform. 'variant': for each SNP in --vcf
                           overlapping the gene, whether it flips a
                           transcript's NMD status.  [required]
  --vcf PATH               Path to VCF file (bgzipped & tabix indexed).
                           Required for --mode=variant.
  -o, --output DIRECTORY   Directory for output tables  [required]
  --help                   Show this message and exit.
```

### 6. Conservation & Missense Pathogenicity Scoring
Scores nucleotide conservation (phyloP/phastCons) and/or AlphaMissense missense
pathogenicity over the genomic span of each pblock/cblock in an isoform-vs-anchor
comparison. Requires `pip install biosurfer[conservation]` for the conservation half;
AlphaMissense scoring uses `pysam` (already a base dependency).

```shell
biosurfer analyze_conservation -d gencode_toy --gene CRYBG2 --bigwig -o results_folder
```

```
Usage: biosurfer analyze_conservation [OPTIONS]

  Scores a gene's alternative-vs-anchor isoform differences
  (pblocks/cblocks) for nucleotide conservation (phyloP/phastCons) and/or
  AlphaMissense missense pathogenicity, over the genomic span of each block.

Options:
  -v, --verbose            Print verbose messages
  -d, --db_name TEXT       Database name  [required]
  --gene TEXT              Target gene to analyze (e.g. PPARG)  [required]
  --bigwig TEXT            Path or URL to a phyloP/phastCons bigWig
                           conservation track. Pass with no value to use
                           the UCSC hg38 100-way phyloP track; omit
                           entirely to skip conservation scoring.
  --no-conservation        Skip conservation scoring even if --bigwig is not
                           given (default behavior).
  --alphamissense PATH     Path to a bgzipped + tabix-indexed AlphaMissense
                           TSV.
  -o, --output DIRECTORY   Directory for output tables  [required]
  --help                   Show this message and exit.
```

### 7. Splice-Site-Strength Scoring
Scores SpliceAI donor/acceptor probability at every splice junction unique to an
alternative isoform or unique to the gene's anchor isoform, to compare predicted
splice-site strength between them. Requires `pip install biosurfer[splicing]` and an
indexed reference genome FASTA matching the database's assembly.

> SpliceAI's trained models are distributed under CC BY-NC 4.0 (academic/non-commercial
> use only); see [github.com/Illumina/SpliceAI](https://github.com/Illumina/SpliceAI).

```shell
biosurfer analyze_splicing -d gencode_toy --gene CRYBG2 --genome-fasta hg38.fa -o results_folder
```

```
Usage: biosurfer analyze_splicing [OPTIONS]

  Scores SpliceAI donor/acceptor probability at every splice junction
  unique to an alternative isoform or unique to the gene's anchor isoform,
  to compare predicted splice-site strength between them.

Options:
  -v, --verbose            Print verbose messages
  -d, --db_name TEXT       Database name  [required]
  --gene TEXT              Target gene to analyze (e.g. PPARG)  [required]
  --genome-fasta PATH      Indexed reference genome FASTA (.fai companion
                           required)  [required]
  --context INTEGER        SpliceAI receptive-field context size (bp)
                           [default: 10000]
  -o, --output DIRECTORY   Directory for output tables  [required]
  --help                   Show this message and exit.
```

### 8. Intrinsic Disorder Scoring
Scores intrinsic disorder (metapredict V3) for every alternative isoform vs. the gene's
anchor isoform, over each protein/codon alignment block. Complementary to (not a
replacement for) the MobiDB-lite IDR annotations already loaded via
`load_feature_mappings`. Requires `pip install biosurfer[disorder]`.

```shell
biosurfer analyze_disorder -d gencode_toy --gene CRYBG2 -o results_folder
```

```
Usage: biosurfer analyze_disorder [OPTIONS]

  Scores intrinsic disorder (metapredict V3) for every alternative isoform
  vs. the gene's anchor isoform, over each protein/codon alignment block,
  to see whether alternative splicing disrupts disordered regions.

Options:
  -v, --verbose            Print verbose messages
  -d, --db_name TEXT       Database name  [required]
  --gene TEXT              Target gene to analyze (e.g. PPARG)  [required]
  -o, --output DIRECTORY   Directory for output tables  [required]
  --help                   Show this message and exit.
```

### 9. Interactive Isoform Viewer
Renders isoform exon/CDS/UTR structure (with start/stop-codon markers) as a
hoverable/zoomable/pannable HTML figure -- a companion to the static `plot` command
above, not a replacement for it. Uses `plotly`, already a base dependency.

```shell
biosurfer plot_interactive -d gencode_toy --gene CRYBG2 -o results_folder
```

```
Usage: biosurfer plot_interactive [OPTIONS] [TRANSCRIPT_IDS]...

  Plot isoforms from a single gene (or a list of TRANSCRIPT_IDS) as an
  interactive, hoverable/zoomable HTML figure -- a companion to the
  existing static `plot` command, not a replacement for it.

Options:
  -v, --verbose           Print verbose messages
  -o, --output DIRECTORY  Directory in which to save the interactive plot
  -d, --db_name TEXT      Database name  [required]
  --gene TEXT             Name of gene for which to plot all isoforms;
                          overrides TRANSCRIPT_IDS
  --help                  Show this message and exit.
```