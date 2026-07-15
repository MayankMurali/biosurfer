# Changelog

All notable changes to this fork are logged here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/). See `CLAUDE.md` in the parent
`biosurfer_analysis` folder for the full working log behind these entries.

## [v0.4.0] - 2026-07-15

### Verification
- `pytest tests -q`: 188 passed, 0 failed (164 baseline + 24 new across the four
  modules below).
- `golden_master/compare_to_baseline.sh`: CLEAN, only the two known-noise diffs.
- `biosurfer analyze_conservation --help` / `analyze_splicing --help` /
  `analyze_disorder --help` / `plot_interactive --help`: all four render correctly.
- Two real bugs caught and fixed during Claude Code CLI verification, not present in
  the original PyPI/conda-free Cowork draft's syntax-check-only pass:
  - `tests/test_splicing.py` and `tests/test_disorder.py` each had a test-only bug:
    their `pytest.raises(match=...)` regex chained
    `.replace('[', r'\[').replace(']', r'\]')` onto an already-escaped raw string,
    double-escaping the brackets so the regex no longer matched. The underlying
    `ImportError` messages were correct throughout. Fixed both to use the plain
    escaped raw string directly.
  - `analyze_conservation`'s `--bigwig` option had `is_flag=False,
    flag_value=DEFAULT_PHYLOP_HG38_URL, default=None`. Click 8.3.1 checks `self.default
    is UNSET` (not `is None`) to decide whether an option needs its "pass with no
    value" behavior -- explicitly passing `default=None` defeated that check, so
    `--bigwig` silently swallowed the next CLI token (e.g. `-o`) as its own value
    instead of falling back to the default phyloP URL. Reproduced even against Click's
    own documented example for this pattern, confirming it's a real Click-version
    interaction, not user error. Fixed by dropping the redundant `default=None` (Click
    already defaults non-required options to `None`). Verified `--bigwig` with no
    value now works correctly both mid-argument-list and at the end.

### Added — Stage B: second new analysis module
- `biosurfer/analysis/conservation/` — nucleotide conservation (phyloP/phastCons) and
  AlphaMissense missense-pathogenicity scoring over pblock/cblock genomic spans, zero
  edits to `core/`:
  - `predict.py`: `score_genomic_range`/`score_protein_range_conservation` (bigWig,
    via the new *optional* `pyBigWig` dependency -- `pip install
    biosurfer[conservation]`, base install unaffected) and
    `query_alphamissense`/`query_alphamissense_for_protein_range` (bgzipped +
    tabix-indexed AlphaMissense TSV, via `pysam` -- already a hard dependency, so this
    half needs no new dependency at all). Defaults to UCSC's remote hg38 phyloP100way
    bigWig URL, which `pyBigWig` queries by HTTP range request rather than downloading
    the multi-GB genome-wide file.
  - `report.py`: gene-level TSV report, mirroring `genetics_analyzer.py`'s
    alignment-block iteration and DataFrame-to-TSV pattern.
  - New CLI subcommand `biosurfer analyze_conservation -d <db> --gene <gene>
    [--bigwig <path/url>] [--alphamissense <tabix path>] -o <output_dir>`, added via
    `cli.add_command` -- no existing command changed.
  - `tests/test_conservation.py` -- conservation tests use a duck-typed fake bigWig
    handle (no real `pyBigWig` needed to run them); AlphaMissense tests build a real
    tiny bgzipped+tabix-indexed TSV via `pysam` rather than mocking it. Coordinate
    math (protein range -> transcript range -> per-exon genomic ranges) was
    independently re-derived with a standalone script before writing the multi-exon
    aggregation test, catching an arithmetic slip in an early draft.
  - `setup.py` gained its first `extras_require` entry (`conservation`) for this.
    See the Verification section above for full pytest/golden-master/CLI results.

### Added — Stage B: third and fourth new analysis modules
- `biosurfer/analysis/splicing/` — SpliceAI donor/acceptor probability at every splice
  junction unique to an alternative isoform vs. unique to the anchor, answering "is this
  alternative splice site predicted stronger/weaker than the sites the anchor isoform
  actually uses?" Zero edits to `core/`.
  - `predict.py`: `score_position`/`score_junction`/`compare_isoform_splice_sites`, via
    the new *optional* `spliceai` + `tensorflow` dependencies (`pip install
    biosurfer[splicing]`). Fetches a real (not synthetic-N-padded) genomic window
    around each position via `pysam.FastaFile` -- already a hard dependency -- and
    reverse-complements it for minus-strand positions before scoring.
  - New CLI subcommand `biosurfer analyze_splicing -d <db> --gene <gene> --genome-fasta
    <path> [--context <bp>] -o <output_dir>`.
  - **Licensing note**: SpliceAI's trained models are CC BY-NC 4.0 (academic/
    non-commercial use only); its source is under the PolyForm Strict License. Fine for
    this project's academic use, flagged since it's a real constraint on downstream use.
  - `tests/test_splicing.py` -- window-fetching/reverse-complement tests against a real
    pysam-indexed FASTA; model-ensemble averaging and junction set-difference labeling
    tested via monkeypatched stand-ins (no real Keras/TensorFlow needed to run them).
- `biosurfer/analysis/disorder/` — intrinsic disorder (metapredict V3) scoring per
  protein/codon alignment block, comparing anchor vs. alternative isoform disorder over
  the same block. Zero edits to `core/`.
  - **Scope note**: `core.data_loaders.load_feature_mappings` already loads MobiDB-lite
    IDR annotations as `ProteinFeature` rows (same pipeline as Pfam domains), and
    `core.alignments` already projects them onto alternative isoforms via
    `ProjectedFeature`/`altered_residues` -- confirmed this is live, working
    infrastructure (not dormant like `ORF.nmd` was) before writing anything. This module
    is a deliberately independent, complementary data source (a different, newer
    ML-based per-residue score) rather than a duplicate of that machinery, and doesn't
    touch the database at all.
  - `predict.py`: `score_disorder`/`compare_isoform_disorder`, via the new *optional*
    `metapredict` dependency (`pip install biosurfer[disorder]`). Scores each protein's
    *full* sequence once (real flanking context everywhere) rather than isolated
    per-block fragments, then slices the cached array per block.
  - New CLI subcommand `biosurfer analyze_disorder -d <db> --gene <gene> -o <output_dir>`.
  - `tests/test_disorder.py` -- exercises the real alignment-block iteration against the
    toy GENCODE `session` fixture's actual CRYBG2 isoforms, with `score_disorder`
    monkeypatched to a deterministic stand-in (no real `metapredict`/PyTorch needed).
- `setup.py` extras_require gained `splicing` (`spliceai`, `tensorflow`) and `disorder`
  (`metapredict`) groups, alongside the existing `conservation` group.

### Added — interactive visualization (companion to the existing static plot)
- `biosurfer/plots/interactive.py` + `biosurfer/plots/interactive_cli.py` — a new
  `biosurfer plot_interactive -d <db> [--gene <gene> | TRANSCRIPT_IDS...] -o
  <output_dir>` subcommand that renders isoform exon/CDS/UTR structure (with
  start/stop-codon markers) as a hoverable/zoomable/pannable HTML figure via Plotly,
  alongside -- not replacing -- the existing static `plot` command. `plotly` is already
  a base `install_requires` dependency, so this needed no new dependency at all.
  `plots/plotting.py` and `plots/canvas.py` are completely untouched; the existing
  `plot` command's behavior is unchanged.
  - v1 scope is exon/CDS/UTR structure only -- pblock/cblock/domain overlays aren't
    included yet (flagged as a natural follow-up, not silently dropped).
  - `tests/test_interactive.py` -- hand-built Transcript/Exon/ORF/Protein fixtures with
    *real* codon sequences (ATG/GCT-repeat/TAA) so `Protein.residues`' amino-acid-to-
    nucleotide linking actually succeeds and start/stop-codon coordinates resolve;
    `plotly` is a hard dependency so, unlike the other three modules above, this test
    file needs no monkeypatching to run.

## [v0.3.0] - 2026-07-15

### Added — Stage B: first new analysis module
- `biosurfer/analysis/nmd/` — nonsense-mediated decay (NMD) prediction, surfacing the
  pre-existing `ORF.nmd` 50-nt-rule property (was implemented but never called from
  anywhere in the codebase). Three additive capabilities, zero edits to `core/`:
  - `predict.py`: `get_nmd_status` (per-transcript), `compare_isoform_nmd`
    (anchor-vs-alternative NMD delta), `predict_variant_nmd_effect` (SNP-triggered NMD
    status flip; indels explicitly rejected, not silently mishandled).
  - `report.py`: gene-level TSV reports for all three modes, mirroring
    `genetics_analyzer.py`'s structure. The variant mode reads a VCF directly via
    `pysam` with no database writes.
  - New CLI subcommand `biosurfer analyze_nmd -d <db> --gene <gene> --mode
    {transcript,isoform,variant} [--vcf <path>] -o <output_dir>`, added via
    `cli.add_command` exactly like `illustrate` -- no existing command changed.
  - `tests/test_nmd.py` -- 9 tests against hand-built Transcript/ORF/Exon fixtures.
  - **Verification status**: `pytest tests/test_nmd.py -q` (9 passed);
    `pytest tests -q` (164 passed, 0 failed); `golden_master/compare_to_baseline.sh`
    (CLEAN -- only the two known-noise diffs, confirming this module touched nothing
    existing); `biosurfer analyze_nmd --help` and both transcript/isoform modes run
    against the toy GENCODE dataset (`--gene CRYBG2`), producing real output TSVs.

## [v0.2.0] - 2026-07-15

### Changed — Stage A modernization (no observable behavior change)
- Removed `core/alignments_old.py`: a superseded implementation with zero references
  anywhere in either repo, confirmed via repo-wide grep before deletion.
- Removed dead code and unused imports across `core/alignments.py`, `core/database.py`,
  `core/splice_events.py`, `core/models/{features,genetics,nonpersistent}.py`,
  `analysis/*.py`, `plots/plotting.py`, and `biosurfer.py`.
- Restructure Phase 1: split `core/helpers.py`'s grab-bag of utilities into focused
  modules (`core/collections_utils.py`, `core/algorithms.py`, `core/io_helpers.py`);
  `core/helpers.py` is now a thin re-export shim so existing
  `from biosurfer.core.helpers import X` call sites keep working unchanged.
- Restructure Phase 2: decomposed the `Database` god-class into a thin lifecycle class
  plus dedicated loader functions/modules, with every public method's name and
  signature preserved.
- Restructure Phase 3: split `TranscriptAlignment.from_transcripts` and
  `CodonAlignment.from_proteins` (`core/alignments.py`) into named helper
  methods/functions (`_map_events_to_intervals`, `_build_compound_event_maps`,
  `_build_alignment_blocks`, `_tx_blocks_to_orf_coords`, `_classify_codon_blocks`,
  `_adjust_boundary_edges`, `_merge_codon_blocks`,
  `_link_codon_blocks_to_transcript_blocks`, `_compare_ranges`/`_split_paired_ranges`),
  with both classmethods' public signatures unchanged.
- Restructure Phase 4: split `plots/plotting.py` into `plots/canvas.py` (`PlotCanvas`:
  generic drawing primitives) and a slimmed `plots/plotting.py`
  (`IsoformPlot(PlotCanvas)`: isoform/alignment-specific drawing). `generate_subtracks`
  moved to `core/algorithms.py` (Phase 1) and is re-exported from `plots/plotting.py`
  for existing callers.
- Restructure Phase 5: extracted `analysis/genome_wide_alignment_analysis.py`'s
  `process_gene` closure into an explicit module-level `_process_gene` function with
  named parameters instead of closure capture.
- `helpers.py`/`io_helpers.py`'s `get_ids_from_gencode_fasta`: added mouse Ensembl ID
  support (`ENSMUST`/`ENSMUSP`, alongside human `ENST`/`ENSP`) and a `None` default
  instead of raising an unhandled `StopIteration` on unmatched headers.
- Removed hardcoded `PPARG`/`T2D` defaults from the GWAS N-terminal-risk analysis code
  and parameterized the trait label (`--trait` CLI option).
- Bumped the Python floor from 3.9 to 3.10 to match what's actually tested/supported by
  the `graph-tool` dependency; updated `setup.py` and the README accordingly.
- Pinned `biopython < 1.87` to avoid a FASTA-comment `ValueError` regression surfaced
  only when testing against a truly fresh dependency resolution.
- Added the `illustrate` CLI subcommand (`biosurfer illustrate`), backed by
  `analysis/illustrate_figs.py`, plus supporting `analysis/frameshift_positions.py` and
  `analysis/st_positions.py` modules.
- Removed a byte-identical duplicate `analyze_nterm` CLI command definition left over
  from a merge-conflict resolution (Click was silently using only the second copy;
  functionally inert, but dead code).

### Fixed
- `.github/workflows/python-package-conda.yml` never installed `graph-tool`,
  `sqlalchemy`, `biopython`, or `biosurfer` itself before running `pytest` — every test
  would have failed at collection (`import biosurfer...`). Added a step that runs
  `conda install --channel conda-forge graph-tool` + `pip install --editable .` first.

### Added
- `golden_master/` — a regression harness that runs the toy-GENCODE pipeline
  (`load_db` → `hybrid_alignment` → `plot`) and diffs the result against a captured
  baseline. This is the permanent Stage A/B regression gate; see
  `golden_master/README.md`.
- Sourced the pytest fixtures `biosurfer/data/gencode/{pfamA.tsv,prosite.dat,
  grch38-protein-features.tsv}` that `tests/conftest.py` needs — the full suite
  (`pytest tests -q`, 155 tests) now runs and passes instead of erroring at collection.

### Unchanged by design
- CLI surface (`load_db`, `hybrid_alignment`, `plot`, `illustrate`, `analyze_nterm` —
  every flag, default, and help string), database schema, output TSV formats
  (pblocks/tblocks/cblocks), and the `conda create -c conda-forge ... && conda install
  -c conda-forge graph-tool && pip install --editable biosurfer` install workflow are
  all unchanged. Verified via a golden-master regression harness (schema/row-count
  snapshot, every TSV, and a rendered plot PNG, byte/visually compared before and after
  every change) and the full 155-test pytest suite passing throughout. This release is
  internal cleanup only — no new user-facing behavior.
