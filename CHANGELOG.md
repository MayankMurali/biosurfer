# Changelog

All notable changes to this fork are logged here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/). See `CLAUDE.md` in the parent
`biosurfer_analysis` folder for the full working log behind these entries.

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
