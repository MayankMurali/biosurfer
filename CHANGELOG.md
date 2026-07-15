# Changelog

All notable changes to this fork are logged here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/). See `CLAUDE.md` in the parent
`biosurfer_analysis` folder for the full working log behind these entries.

## [Unreleased]

### Changed — Stage A modernization (no observable behavior change)
- Removed `core/alignments_old.py`: a superseded implementation with zero references
  anywhere in either repo, confirmed via repo-wide grep before deletion.
- Split `TranscriptAlignment.from_transcripts` (`core/alignments.py`) into
  `_map_events_to_intervals`, `_build_compound_event_maps`, `_build_alignment_blocks`.
- Split `CodonAlignment.from_proteins` (`core/alignments.py`) into
  `_tx_blocks_to_orf_coords`, `_classify_codon_blocks`, `_adjust_boundary_edges`,
  `_merge_codon_blocks`, `_link_codon_blocks_to_transcript_blocks`, plus module-level
  `_compare_ranges`/`_split_paired_ranges`.
- Split `plots/plotting.py` into `plots/canvas.py` (`PlotCanvas`: generic drawing
  primitives) and a slimmed `plots/plotting.py` (`IsoformPlot(PlotCanvas)`: isoform/
  alignment-specific drawing). Moved `generate_subtracks` out of `plots/plotting.py`.
- <!-- CLAUDE-CODE-CLI: confirm final destination module for generate_subtracks and
  update this line + the "Fixed" CI note below before tagging a release. -->

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
- CLI surface (`load_db`, `hybrid_alignment`, `plot`, `illustrate` — every flag,
  default, and help string), database schema, output TSV formats (pblocks/tblocks/
  cblocks), and the `conda create -c conda-forge ... && conda install -c conda-forge
  graph-tool && pip install --editable biosurfer` install workflow are all unchanged.
  This release is internal cleanup only.

<!--
CLAUDE-CODE-CLI: before tagging, diff this against `git log v0.1.0..HEAD --oneline` and
correct anything above that doesn't match what actually landed. Then replace this
comment and the "[Unreleased]" header with:

## [v0.2.0] - YYYY-MM-DD
-->
