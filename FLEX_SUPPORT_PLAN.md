# Flex Probe-Barcode + CRISPR Guide Capture Support for `sbxpress`

> Handoff plan drafted 2026-09-14 in a `roche-SBX-ABRF-CAT06` planning session,
> saved here (its intended implementation home) per the plan's own "Next step"
> note. See the **Amendment** section at the end for review notes added the
> same day before implementation started.

## Context

CAT06 benchmarks Roche SBX sequencing against Illumina on a 10x Genomics Flex V2 library (9 samples pooled via probe barcodes, plus a paired CRISPR Guide Capture library for a CD45-knockout screen — see project memory `cat06_project_overview.md`). The Illumina arm is fully processed via `cellranger multi` and gives us real, validated ground truth at `CAT06/illumina_data/RZ16443_L1_multi/outs/`. `sbxpress` (the SBX-side tool, `/mnt/analysis/software_shared/sbxpress`, Rust, v0.3.7 source / stale 0.3.5 installed) has **no Flex support today** — it only demuxes at the i7/i5 library level and only understands plain poly-dT 3' GEX (genome+GTF alignment). This plan adds a new, cohesive module so `sbxpress` can process both the Flex GEX and Flex CRISPR streams for the same pooled library, with the explicit goal of **matching cellranger's actual published algorithms** (not just approximate output numbers) wherever a design choice exists, so the eventual SBX-vs-Illumina comparison is apples-to-apples.

This plan is grounded in (a) direct reading of the sbxpress source, (b) empirical scans against the real CAT06 raw SBX reads / probe-barcode whitelist / feature reference, and (c) 10x's published algorithm docs (`cr-flex-frp-algorithm`, `cr-crispr-algorithm`) cross-checked against the real completed `cellranger multi` output files already on disk.

## Key architecture facts (verified against the real source)

- CB(16)+UMI(12) location in `orient.rs`/`extract.rs` is assay-agnostic and unaffected.
- **Probe barcode location, empirically confirmed**: for both GEX and CRISPR reads, the 10bp probe barcode sits at a fixed small window right after `Extracted::polyt_end` (`extract.rs:52-53`) — i.e. immediately after the existing CB+UMI+poly-dT anchor, in `revcomp(whitelist_seq)` orientation in canonical coordinates. **No new anchor search is needed** — this is a new fixed-window read off a coordinate `extract()` already computes. This must still be confirmed with a real Rust harness (Phase 0) before hardcoding, since the confirming pass so far used an ad hoc Python k-mer scan, not the actual `orient()`/`extract()`/emitter code path.
- Because `lib.rs::push_fastq_r2_rc` reverse-complements at emission time, the probe barcode must be located during `extract`'s per-read processing on the uncapped sequence — **not** recovered from the emitted (possibly `--max-r2-len`-truncated) R2 FASTQ, since truncation keeps the wrong end.
- GEX vs. CRISPR library separation is **already solved** — they're different i7/i5 pairs in the existing raw SBX reads, and `sbxpress extract --sample-sheet <bclconvert csv>` (existing `sample_sheet.rs`/`demux.rs`) splits them today, unmodified.
- `demux.rs`'s Hamming-tolerant small-set lookup (`insert_neighbors`/`assign_windowed`) is the right model for probe-barcode and CRISPR-anchor matching — `barcode.rs`'s 16bp-specific bit-packed `Whitelist` is not reusable as-is (should be refactored into a shared generic `HammingIndex` primitive both old and new code builds on).
- `annotate_bam.rs` only touches known tags (`GX,GN,TX,AN,RE,ZG,ZT,XC,XT,XP,NH,HI`) — a new `PR`/`PB` tag pair rides through `annotate`/`molecules` untouched once added at `extract` time via the existing FASTQ-comment → minimap2 `-y` mechanism (same as `CR`/`CB`/`UR`/`UB` today).
- **`matrix --cell-barcodes <WHITELIST>`** already exists and does exactly "restrict/order a matrix to an arbitrary barcode list" — per-sample splitting needs **zero changes** to `cr_bam_to_mtx.rs`/`cell_calling.rs`. `molecules --tagseq` and `cells --tagseq-barcodes` already exist and skip UMI-collision adjudication / statistical cell-calling respectively when the real barcode identity is already known externally — both are a good fit for the CRISPR stream's simpler dedup needs, though see the protospacer-calling gap below.
- `molecules`'s core loop (`sbx10x-molecules.rs`) skips any `record.is_unmapped()` — a synthetic *unmapped* BAM for CRISPR reads would be silently dropped. CRISPR reads must be emitted as **mapped** records against a tiny synthetic per-guide contig space (tid=guide index, trivial CIGAR, mapq 60) to flow through unmodified.
- Reinstall note (do before any field test, not a design item): `target/release/sbxpress` is already 0.3.7; `/usr/local/bin/sbxpress` is a stale 0.3.5 build.

## Where this plan corrects the initial design against published cellranger algorithms

Two places where "just wire the existing pieces together" would silently diverge from what cellranger actually does — confirmed against `cr-flex-frp-algorithm`/`cr-crispr-algorithm` and the real `outs/` files:

1. **Probe-barcode sample assignment is per-molecule against a sample's *declared barcode set* (which can be >1 barcode), not a naive "top single barcode wins" vote.** Published: "the corrected Probe Barcode sequence is mapped to its corresponding Probe Barcode ID, and if this ID has been specified within the `probe_barcode_ids` column of the `[samples]` section, the read is assigned to the sample." Real `config.csv` confirms samples like `Pos,D-A01|D-A02` use *two* probe barcodes per sample. The real `multiplexing_analysis/frp_gem_barcode_overlap.csv` (pairwise GEM co-occurrence between every probe-barcode pair, e.g. `D-A01,D-A02,...,overlap=0.277`) shows 10x uses this co-occurrence structure for **high-occupancy GEM detection** (a droplet showing barcodes from *two different samples'* declared sets, not just its own pair, indicates >1 real cell in that droplet) — this is a distinct, separate mechanism from simple purity thresholding.
2. **CRISPR protospacer calling is a per-guide Gaussian Mixture Model over log(UMI count), with an unconditional minimum-UMI floor — not a raw-count passthrough.** Published: cellranger "models two cell populations for each guide RNA... fits a Gaussian Mixture Model to the log-transformed distribution of Molecules/Cell," and "a cell barcode now requires a minimum of three guide RNA UMIs to be considered positive." **Empirically confirmed this floor is genuinely adaptive, not just a hardcoded 3**: real `protospacer_umi_thresholds.csv` shows `CD45,3` for sample `Pos` but `CD45,6` for sample `RZ211` — the GMM-fit threshold moves with each sample's ambient background. A simple `--min-crispr-umi 3` passthrough (what plain `matrix`+`cells` would give you) would be *wrong* for RZ211 and any higher-background sample. This needs an explicit new calling step.

## Phase 0 — Confirm coordinates with real code (½ day)

New diagnostic example `examples/probe_barcode_locate.rs` (mirrors existing `examples/orient_purity.rs`/`examples/cb_error_profile.rs`): run real `orient()`+`extract()` over a sample of `CAT06/raw_fastqs/*.fastq.gz`, confirm the probe-barcode window position off `polyt_end` (Hamming ≤2 against the real 384 entries in `CAT06/resources/flex-v2-384.txt`), and separately confirm the CRISPR anchor+capture layout (`<anchor>(BC)`) against `CAT06/feature_reference_RZ.csv`'s two guides, on the *emitted* R2 (post `push_fastq_r2_rc`, uncapped `--max-r2-len 0`). This must pass before any anchor constant is hardcoded in Phase 1.

### Phase 0 — results (2026-09-14)

Ran `examples/probe_barcode_locate.rs` (committed) against 4 shards picked at
random from the 754 in `CAT06/raw_fastqs/` (indices 431, 489, 646, 693 —
deliberately *not* the first few shards, which independent QC has flagged as
lower quality than the run average), 100k reads/shard, 340,432 confidently
oriented reads total:

- **Probe barcode: confirmed at `polyt_end + 10`, `revcomp` orientation** —
  62.4% of confidently-oriented reads hit within ≤2 mismatches at this exact
  `(offset, orientation)`, a clear multiple of the next-highest offset (35.1%
  at `+2`, then a long tail down to a ~15–20% noise floor at other offsets —
  consistent with the expected false-hit rate of a 10bp/≤2mm search against
  384 targets on non-matching sequence). **This corrects the plan's original
  assumption of offset 0** ("immediately after polyt_end") — the true offset
  is +10, with real jitter across nearby offsets that tracks `dt_run_end`
  sometimes over/under-extending the detected T-run boundary by a few bp
  (expected: `dt_run_end` tolerates isolated mismatches and the true junction
  isn't itself a homopolymer edge). **Design implication for Phase 1**:
  `ProbeWhitelist::locate` must do a real windowed search centered near
  `polyt_end + 10` (not trust a fixed exact offset), exactly as already
  planned — this result is the empirical center-point for that window, not a
  reason to hardcode a single offset.
- **CRISPR anchor: confirmed at emitted-R2 offset 10, immediately followed by
  a clean capture barcode.** The `CD45` anchor hit offset 10 in >20x more
  reads than any other offset (46,938 vs. 2,301 at the next-highest, offset
  11 — almost certainly boundary fuzz on the same true site, not a second
  real location). At that anchor offset, the following 20bp capture-barcode
  window matched a real guide sequence exactly (mm=0) 95.6% of the time
  matches occurred, with a sharp exact-match spike and only a long thin tail
  from sequencing error — no ambiguity in the anchor+capture layout.
  Guide calls (capture mm≤2): **CD45 48,301 vs. inactive 319** — consistent
  with this being a CD45-knockout screen where `inactive` is a control guide
  present at much lower abundance/only in specific samples.
- Both loci reproduced with the same dominant offset independently across
  all 4 shards — not a single-file artifact.

**Conclusion: Phase 0 passes.** Both empirically-confirmed offsets refine
(rather than overturn) the plan's architecture facts — real orientation and
windowed-search design hold; the one hardcoded assumption (probe barcode at
`polyt_end + 0`) is now corrected to `polyt_end + 10` ahead of Phase 1.

## Phase 1 — Shared `HammingIndex` + `probe_barcode.rs` + `extract` tagging

- Refactor `Slot`/`insert_neighbors`/`enumerate` out of `demux.rs` into `src/hamming_index.rs` (generic over fixed-length byte keys); `demux.rs` becomes a thin wrapper. This is the shared primitive both `probe_barcode.rs` (this phase) and `guide_capture.rs` (Phase 4) build on.
- New `src/probe_barcode.rs`: `ProbeWhitelist::load(flex-v2-384.txt)` building a `HammingIndex` over `revcomp(seq)`; real min pairwise Hamming distance among the 384 barcodes is 3, so `max_mm=1` is provably unambiguous (mirrors `demux.rs`'s own doc comment about i5/i7 min-distance reasoning). `locate(seq, polyt_end, window)` does a windowed search (mirrors `demux::assign_windowed`).
- `config.rs`: add `--probe-barcode-whitelist <FILE>`, `--probe-barcode-max-mm` (default 1), `--probe-barcode-window` (default 10).
- `lib.rs`: thread `probe_wl: Option<&ProbeWhitelist>` through `Shared` (mirrors existing `whitelist`/`demuxer` fields); in `process_one`, after `extract()`, locate+tag `PR:Z:<raw>`/`PB:Z:<corrected>` into the FASTQ comment (gated on `--minimap2` preset, same as `CR`/`CB` today). Zero cost when unset.
- `accounting.rs`: add probe-barcode hit-rate metrics to `metrics_summary.json`, mirroring `WhitelistMetrics`.
- **Validate**: re-run Phase 0's harness at scale; confirm hit-rate is in a sane range and no throughput regression when the flag is unset.

### Phase 1 — results (2026-09-14)

Implemented as designed, with one orientation bug caught by unit tests before
it ever reached real data: `HammingIndex` keys are built as
`revcomp(whitelist_seq)` (the orientation the barcode physically appears in
on a canonical read), so lookup must query with the raw canonical window
directly — the first draft of `locate()` revcomp'd the window *before*
querying, a double-revcomp that matched nothing. Caught immediately by
`probe_barcode::tests::locates_exact_match_at_confirmed_offset` failing;
fixed by querying with the canonical window and reserving the revcomp for
the reported `raw` field (kept in whitelist orientation for readability).

Demux refactor regression check: full existing suite green before and after
(148 lib tests, all `demux::` tests unchanged and passing, +11 new tests for
`hamming_index`/`probe_barcode`), no clippy warnings introduced.

Ran the real `sbx10x` binary (not the Phase 0 diagnostic) with
`--minimap2 --probe-barcode-whitelist flex-v2-384.txt` against the **full**
4 shards used in Phase 0 (not just their 100k-read samples) — 176,072,967
reads, 143,099,522 pairs emitted, 27.2 min wall / 107,862 reads/sec:

- **Probe-barcode hit rate: 92.7%** (132,605,115 / 143,099,522 pairs) — well
  above Phase 0's single-fixed-offset 62%, confirming the windowed search
  (±10bp around `polyt_end+10`) is doing real work absorbing per-read
  `dt_run_end` jitter, not just replicating the diagnostic's ceiling.
- **Cross-checked against ground truth already on disk**: tallied `PB` id
  frequency over a 2M-read sample of the emitted R2 FASTQ. The **top 15 ids
  by count are exactly the 15 probe barcodes declared in
  `config.csv`'s `[samples]` section** (`D-A01`–`D-A08`, `D-B01`–`D-B07`),
  with a sharp ~7× cliff to the noise floor (15th place `D-A07` = 36,407
  vs. 16th place `A-E06` = 5,037, then a long thin tail across the
  remaining 369 barcodes at low counts — background mismatch noise, not
  real samples). This is a real-world confirmation using data already on
  disk, no need to wait for Phase 2/6's BAM-level ground truth comparisons.
- Throughput-regression check (flag unset) tracked separately — see below.

**Conclusion: Phase 1 passes**, including the free cross-check against
`config.csv` proposed when starting this phase.

### CB/UMI ground-truth validation (2026-09-14, ad hoc — not new code)

Prompted by "can we validate CB and UMI at this juncture" — answered with
data already on disk, no new pipeline code:

- Built a whitelist from the 735,036 distinct 16bp prefixes of
  `raw_feature_bc_matrix/barcodes.tsv.gz` — the real GEM barcodes cellranger
  observed in *this exact* CAT06 run (a tighter ground truth than any
  generic 10x whitelist, e.g. `3M-february-2018.txt`/`737K-arc-v1.txt.gz` —
  neither of which is confirmed as CAT06's actual chemistry). Ran
  `sbx10x extract --whitelist <this> --probe-barcode-whitelist
  flex-v2-384.txt` on an 800k-read sample (200k/shard, same 4 shards): **CB
  hit rate 81.6%** (on-whitelist + within-1-corrected).
- **Joint CB+PB reconstruction** (stronger check): cellranger's Flex
  `raw_feature_bc_matrix` barcode string is literally
  `CB(16bp) + probe_barcode_seq(10bp) + "-1"` — confirmed directly (row 1's
  trailing 10bp decoded to `D-A01` in `flex-v2-384.txt`, matching the same
  id our own `probe_barcode.rs` would report). Reconstructed that exact
  string from our independently-extracted `CB` + `PB`-id-looked-up-seq and
  checked *exact* membership against the real 8.9M-row raw barcode set:
  **73.5%** (194,108 / 264,069) matched a real observed GEM×probe-barcode
  combination — with zero mismatch tolerance on either segment, so this is
  a meaningfully strict joint test, not two independent checks added
  together.
- **UMI**: no ground-truth UMI sequence exists (random by construction), but
  structural correctness is exact: the full Phase 1 run's
  `r1_eq_28 = 143,099,522` / `r1_not_28 = 0` — 100% of emitted R1 blocks are
  the correct 16+12 = 28bp, zero exceptions across 143M pairs. UMI
  *deduplication* is out of scope here — it's the pre-existing, unmodified
  `molecules`/`--tagseq` path; this plan doesn't touch it.
- Not investigated further (out of scope for this check): `CB gate`'s
  `located_confident_pct` came back oddly low (0.1%) on this Flex sample —
  worth a closer look before Phase 2, since that heuristic (dual-anchor:
  R1-tail + polyT-flank confirmation) may need retuning for Flex's shorter
  fixed spacer vs. the GEX assay it was tuned against.

## Phase 2 — `probe_sample_sheet.rs` + `probe-demux` subcommand (revised for published-algorithm parity)

- New `src/probe_sample_sheet.rs`: parses the real cellranger `[samples]` section verbatim (`sample_id,probe_barcode_ids,...`, pipe-delimited multi-barcode groups) — primary format, so a user can point `sbxpress` straight at an existing `config.csv`; flat `sample_id,probe_barcode_ids` CSV as fallback.
- New `src/probe_demux.rs` + `src/bin/sbxpress-probe-demux.rs` (new `[[bin]]`, new dispatcher arm). CLI: `sbxpress probe-demux <BAM> --probe-sample-sheet <CSV> -o <DIR> [--min-molecules 10]`.
- Input: the GEX-library BAM after `annotate`+`molecules` (gate on the `xf` representative-molecule bit to tally per-*molecule*, not per-*read*, avoiding PCR-duplicate inflation).
- **Assignment logic (revised)**: for each `CB`, tally molecule counts per `PB`. Map each sample's declared `probe_barcode_ids` (one or more) from the sample sheet. Assign `CB` to the sample whose declared barcode set captures the dominant share of that CB's molecules (not just a single top barcode — sum counts across a sample's whole declared set first, then compare across samples). Flag `Unassigned` below `--min-molecules`.
- **New: high-occupancy GEM diagnostic**, mirroring `frp_gem_barcode_overlap.csv` — for every pair of probe barcodes across the whole experiment, compute `common_gems / min(barcode1_gems, barcode2_gems)` (or matching overlap definition) and emit an equivalent `gem_barcode_overlap.csv`. Use it to flag/exclude CBs whose molecule mixture spans barcodes from *two different samples'* declared sets (true cross-sample doublet signature) as `HighOccupancy`, distinct from the `Unassigned` (low-confidence/noisy) bucket.
- Outputs: `sample_assignment.csv` (`cell_barcode,sample_id,probe_barcodes,molecules,purity,status` where status ∈ {Assigned, Unassigned, HighOccupancy}), one `<sample_id>.barcodes.txt` per sample (feeds `matrix --cell-barcodes` directly), `gem_barcode_overlap.csv`, `probe_demux_metrics.json`.
- **Validate**: compare `<sample_id>.barcodes.txt` counts and `gem_barcode_overlap.csv` directly against the real `multiplexing_analysis/frp_gem_barcode_overlap.csv` and `cells_per_tag.json` at `CAT06/illumina_data/RZ16443_L1_multi/outs/`.

### Phase 2 — results (2026-09-14)

Implemented as designed: `probe_sample_sheet.rs` (mirrors `sample_sheet.rs`'s
`[section]`-parsing style — auto-detects a real `config.csv` by the presence
of a `[samples]` header, flat CSV fallback otherwise), `probe_demux.rs`
(`assign`/`gem_overlap` as pure, unit-tested functions over a `(CB, PB) ->
molecule count` tally, plus a thin BAM-reading/CSV-writing shell), new
`sbxpress-probe-demux` binary wired into the dispatcher. 20 new unit tests
(8 for `probe_demux`, 5 for `probe_sample_sheet`, both real-shaped fixtures
— e.g. `probe_sample_sheet`'s test literally reproduces `config.csv`'s
`[gene-expression]`/`[libraries]`/`[samples]` layout to prove unrelated
sections are correctly skipped). Full suite green throughout (161 lib tests),
no clippy warnings in new code.

**No real annotated+molecule BAM exists yet for CAT06** (that needs a full
alignment pass against a reference genome — out of scope for this phase, and
properly Phase 6's job), so full BAM-in end-to-end validation as originally
scoped isn't possible yet. Worked around it two ways:

1. `read_bam_tally`'s BAM-reading shell (the only piece that actually needs a
   BAM) is unit-tested directly against a synthetic BAM built with
   `rust-htslib` (mirrors `annotate_bam.rs`'s own synthetic-BAM test
   pattern), confirming the `xf`-bit-8 gate and CB/empty-PB skip logic.
2. For a **real-data** exercise of the actual `assign`/`gem_overlap` logic
   (not just synthetic unit tests), built a one-off diagnostic
   (`examples/probe_demux_from_fastq.rs`) that constructs the same `(CB, PB)
   -> molecule count` tally directly from a `--minimap2`-preset R2 FASTQ's
   `CB`/`UB`/`PB` comment tags (deduped by `(CB, UB, PB)` as a proxy for the
   real `xf`-gated count — no gene gate, since this predates
   alignment/annotation, so it's a slight overcount vs. true per-gene
   molecules, but close enough to exercise the real code against real data),
   then calls the actual `probe_demux::assign`/`gem_overlap` functions and
   `ProbeSampleSheet::load` on the real `config.csv`. Run against the 800k-read
   CB-validation sample (200k/shard, `--min-molecules 3`):
   - **Sample ranking is directionally correct against real ground truth.**
     Real per-sample true-cell counts (`sample_filtered_barcodes.csv`) rank
     `RZ213NR` (15,728) / `RZ213R` (20,866) / `RZ213RD` (20,732) as a clear
     minority group against the other 6 samples (50k–75k each). Our raw,
     pre-cell-calling `Assigned` counts on this small sample reproduce that
     exact grouping (276 / 279 / 899 vs. 2,537–9,697 for the other six) —
     not an exact match (expected: small non-random sample, no real
     cell-calling yet, that's Phase 3's job on top of this), but the
     qualitative structure is right.
   - **Found and fixed a real bug via this exercise, not a unit test**: the
     first `gem_overlap` run's top-10-by-overlap rows were *all*
     `overlap == 1.0` artifacts from single-read background/error probe-id
     calls coincidentally sharing one droplet with a real barcode — `tally`
     holds unfiltered raw PB calls, and a `barcode1_gems == 1` row hits
     `overlap = common/min(1, big) = 1.0` trivially. Fixed by restricting
     `gem_overlap` to the sample sheet's declared probe-barcode set (new
     `declared: &HashSet<&str>` parameter); re-running produced sane,
     stable top rows across all 15 real barcodes (~0.14–0.15 overlap).
   - **Compared against the real `frp_gem_barcode_overlap.csv` directly**:
     all 105 real pairs sit in a *continuous* 0.15–0.28 band with no
     separated outlier cluster — meaning aggregate pairwise-overlap
     thresholding isn't how cellranger's real
     `filter-high-occupancy-gems` actually flags doublets; it must operate
     per-GEM. This is a point *in favor* of this plan's original design
     choice (`Status::HighOccupancy` is already a self-contained per-CB
     rule, independent of the pairwise table) rather than something to
     rework. Our own values (~0.14–0.15) are in the same ballpark but lower
     — expected, given the small non-random sample and no real
     UMI/gene-level dedup yet.

**Conclusion: Phase 2's core logic is implemented and exercised against real
data with a plausible, ground-truth-consistent result** — but not yet
"passed" in the same rigorous sense as Phase 0/1, since that requires the
real annotated+molecule BAM this phase was actually designed to consume.
Revisit with a proper end-to-end run once Phase 6's alignment step exists.

## Phase 3 — GEX per-sample matrix/cells (zero new Rust code)

Per sample: `sbxpress matrix <gex.bam> --cell-barcodes <sample>.barcodes.txt -o <sample>/gex_raw`, then `sbxpress cells <sample>/gex_raw/raw_feature_bc_matrix.h5 -o <sample>/gex_outs` (normal OrdMag/EmptyDrops cell-calling per sample — reuses existing, fully general commands per the architecture facts above). This phase is the orchestration script skeleton (`scripts/run_flex_multi.sh`, new file in the sbxpress repo) plus validation, no new Rust.

### Phase 3 — results (2026-09-14)

"Zero new Rust" turned out to be almost, not quite, true — wiring Phase 2's
real output into Phase 3's real commands (not just reading the plan)
surfaced two genuine bugs, both only visible by actually chaining the real
binaries together end-to-end on a small synthetic BAM (built with
`samtools`/hand-written SAM, `rust-htslib`-tag-compatible):

1. **Phase 2 bug**: `<sample>.barcodes.txt` was writing `CB` values with
   their `-1` gem-group suffix still attached (our own `CB` tag convention).
   But `matrix --cell-barcodes` (`cr_bam_to_mtx.rs::read_cell_barcodes`)
   expects *bare* ACGT barcodes and appends its own `-{gem_group}` suffix —
   it rejects anything else as "not a valid 10x reference barcode". Every
   Phase 3 `matrix` call would have failed immediately. Fixed in
   `probe_demux.rs` with a `strip_gem_suffix` helper (strips a trailing
   `-<digits>`, not hardcoded to `-1`, in case of a future multi-gem-group
   run) applied only to the `.barcodes.txt` output — `sample_assignment.csv`
   keeps the `-1` as-is, since that's meant to read directly against the
   BAM's own `CB` tag. New tests: `strip_gem_suffix_removes_trailing_digits_only`,
   `per_sample_barcode_files_are_bare_acgt_matrix_cell_barcodes_ready`.
2. **Pre-existing `cr_bam_to_mtx.rs` gap, not Flex-specific**: `matrix
   --cell-barcodes` treats its list as a *strict, exhaustive* whitelist —
   any barcode counted from the BAM that isn't on the list is a hard error,
   not a skip. That's the right behavior for its original single-sample use
   case (an unlisted barcode there usually means a real whitelist/BAM
   mismatch worth failing loudly on). But Phase 3's whole premise is one
   **pooled**, multi-sample GEX BAM split *after the fact* by probe barcode
   — every per-sample `matrix` call legitimately sees plenty of other
   samples' barcodes too, and would hard-error on the first one. Fixed by
   adding an opt-in `--skip-unlisted-barcodes` flag (default off, so the
   existing single-sample behavior/tests are untouched) that skips instead
   of erroring, and reports how many molecules were skipped in
   `metrics_summary.json`. New tests: `unlisted_barcode_errors_by_default`,
   `skip_unlisted_barcodes_drops_them_instead_of_erroring`.

**Follow-up (2026-09-14)**: asked whether `-1` is actually the typical
`--cell-barcodes` input shape, since Phase 2's bug (above) was about
stripping it. Answer: no, one input source (10x's own raw whitelist files,
e.g. `3M-february-2018.txt`) is bare; but another very plausible one
(a barcode list exported from a real run — `barcodes.tsv.gz`,
`filtered_barcodes.csv`, or our own `probe-demux`'s CB values before this
session's fix) legitimately carries the suffix. `read_cell_barcodes`
previously hard-rejected the latter. Made it tolerate both — mirroring
`cell_calling.rs`'s pre-existing `resolve_tagseq_barcodes`, which already
handles this same ambiguity the same way — so `--gem-group` is always what
gets applied regardless of what the input file happened to contain. New
test: `an_existing_gem_group_suffix_is_accepted_but_gem_group_flag_wins`.

With both fixed, ran the real chain — `sbxpress probe-demux` →
`sbxpress matrix --skip-unlisted-barcodes` → `sbxpress cells` — through the
actual compiled binaries (not library calls) against a synthetic 2-sample,
2-gene pooled BAM (real `samtools`-built), and it produced real
`filtered_feature_bc_matrix.h5`/`filtered_barcodes.csv` per sample with the
correct, cross-contamination-free per-sample counts (verified: a molecule
tagged with an undeclared probe barcode landed in the right sample's count
via its `CB`, exactly as it should). `scripts/run_flex_multi.sh` runs this
same chain per sample from `probe-demux`'s output directory, with a clean
skip for a sample with zero assigned barcodes.

Full suite green throughout (165 lib tests, +4 from this phase), no clippy
warnings introduced.

**Conclusion: Phase 3 passes** — real orchestration script, and the whole
Phase 2 → 3 handoff verified end-to-end through actual binaries, not just
unit tests. Real CAT06 numbers are still blocked on Phase 6's alignment
step, same as Phases 1–2.

## Phase 4 — `guide_capture.rs` + `guides` subcommand (CRISPR read → guide assignment)

> Implemented as `src/feature_capture.rs`/`FeatureMatcher`, not
> `guide_capture.rs`/`GuideMatcher` — see "Phase 4 — results" below for why
> (the Phase 0/1 amendment's generalization requirement).

- New `src/guide_capture.rs`: `FeatureReference::load(csv)` parsing the real `id,name,read,pattern,sequence,feature_type` schema (confirmed against `CAT06/feature_reference_RZ.csv`); parse `pattern` as `<anchor>(BC)` (confirmed layout; flag `(BC)<anchor>` as an explicit unsupported-format error rather than silently mis-parsing). `GuideMatcher` built on Phase 1's shared `HammingIndex`: one index for anchors (windowed, Hamming-tolerant, mirrors `demux::assign_windowed`), then match the preceding `bc_len` bases against a second small index of guides' `sequence` values.
- New `src/bin/sbxpress-guides.rs` (new `[[bin]]`, new dispatcher arm). CLI: `sbxpress guides <R2_FASTQ> --feature-reference <CSV> -o <BAM> [--features-out <TSV>] [--anchor-max-mm 2] [--capture-max-mm 1] [--anchor-window 5]`.
- Input is the `--minimap2`-preset R2 FASTQ from `extract` on the CRISPR-library i7/i5 stream (no `extract.rs` changes needed — its existing dT-run-anchored bracketing already captures the anchor+capture region as a side effect). Parses `CR/CB/UR/UB/PR/PB` back out of the FASTQ comment, matches guides, emits a **mapped** BAM against a synthetic per-guide contig header (tid=guide index, trivial CIGAR, mapq 60 on match; unmapped — correctly skipped by `molecules` — on no match), tags `GX`/`GN`=guide id/name plus carried `CB`/`UB`/`PB`. `--features-out` writes a canonical `id\tname` TSV (mirrors `annotate_bam.rs`'s `--genes-mtx` create-or-verify pattern) for `matrix --features`.
- **Validate**: anchor/capture hit rate against `feature_reference_RZ.csv` on real raw reads.

### Phase 4 — results (2026-09-14)

Implemented as `src/feature_capture.rs` (not `guide_capture.rs`) per the
Phase 0/1 amendment's generalization requirement: `FeatureReference`/
`FeatureMatcher` are assay-agnostic (no CRISPR- or Flex-specific branching;
`feature_type` is carried through as data, never dispatched on), `PB` is
genuinely optional (parsed if present, simply absent from the tag list
otherwise — verified by a dedicated test with no `PB:Z:` in the input at
all). `FeatureMatcher` also handles a mixed panel with *different* anchors
per feature (grouping by `(anchor, capture_len)`, one `HammingIndex` per
group, one per distinct anchor length) rather than assuming one universal
anchor — a real scenario (CITE-seq/CMO experiments mix feature types with
different anchors) the original single-anchor design would have silently
mishandled. `sbxpress guides` stays the CRISPR-flavored CLI entry point per
the amendment's suggestion, built on the general matcher.

One deliberate deviation from the plan's literal CLI default: `--anchor-window`
defaults to 30, not the plan's placeholder 5. That placeholder predates
Phase 0's empirical finding (anchor at emitted-R2 offset ~10, with jitter
out to the low teens); 5 would have silently missed the real signal. 30
covers the confirmed real offset with margin for other constructs, without
hardcoding a specific expected offset the way probe-barcode's `locate` does
(unlike the probe barcode, one universal offset isn't something we can
assume holds for arbitrary future Feature Barcoding constructs, so this
scans from the read start rather than centering on a guess).

12 new unit tests, including one specifically proving the generalization
claim rather than just asserting it in a comment
(`handles_a_mixed_panel_with_different_anchors_per_feature_type`: a
CRISPR-flavored feature and an "Antibody Capture"-flavored one with
different anchors, both resolved correctly) and a full `run()`
end-to-end test asserting on the actual written BAM (mapped record has
`GX`/`GN` + carried tags including `PB`; unmapped record has no `GX`/`GN`,
still carries `CB`/`UB`, and correctly has no `PB` when the input never
had one). 181 lib tests total, no clippy warnings.

**Real-data validation**: ran the actual `sbxpress guides` binary (not a
diagnostic) against the same 800k-read CB-validation FASTQ used earlier,
with the real `feature_reference_RZ.csv`:

- 546,491 reads → 83,493 matched (**15.3% hit rate**) — consistent with
  Phase 0's own ad hoc finding (~13.8% anchor hit rate at the 4-shard
  scale); the two independent measurements agree.
- Per-feature: CD45 82,963, `inactive` 530 — same lopsided ratio Phase 0
  found, consistent with a CD45-knockout screen where `inactive` is a
  low-abundance control.
- Verified the BAM directly (not just the metrics JSON): 83,493 mapped +
  462,998 unmapped = 546,491 exactly; a real mapped record shows
  `RNAME:CD45 CIGAR:102M GX:Z:CD45 GN:Z:CD45` plus `CB`/`PB` carried
  through (a real Flex probe barcode landing on a CRISPR-library read,
  as expected — probe barcodes hybridize before GEM partitioning, so they
  identify the sample regardless of which physical library a read comes
  from); a real unmapped record carries `CB`/`UB` but correctly no `GX`/`GN`.
  `@SQ` header has exactly the 2 synthetic contigs (`CD45`, `inactive`),
  both the nominal `LN:1000000`.

**Conclusion: Phase 4 passes.**

## Phase 5 — CRISPR counting + **new protospacer-calling module** (the published-algorithm gap)

Per sample: `sbxpress molecules guides.bam --tagseq -o corrected.bam` → `sbxpress matrix corrected.bam --cell-barcodes <sample>.gex_cells.txt --features guide_features.tsv -o <sample>/crispr_raw` (use the sample's **Phase 3 GEX-called cell list**, not the pre-cell-calling probe-demux list — CRISPR identity piggybacks on GEX-defined cell identity, matching real 10x biology) → `sbxpress cells <sample>/crispr_raw/... --tagseq-barcodes <sample>.gex_cells.txt -o <sample>/crispr_outs`. This gives a raw per-cell-per-guide UMI matrix, but **stops short of cellranger's actual call** — do not treat any nonzero count as "detected."

- **New `src/protospacer_calling.rs` + `sbxpress protospacer-call` subcommand**: input is the per-sample filtered CRISPR matrix + feature reference. For each guide: fit a 2-component Gaussian mixture on `log1p(UMI count)` across all called cells (background/ambient component vs. signal component — standard EM for a 2-component 1D GMM), derive the per-guide UMI threshold from the fitted decision boundary, and **enforce `max(fitted_threshold, --min-umi)` unconditionally** (default `--min-umi 3`, matching cellranger's floor — confirmed genuinely adaptive above/below this floor across real samples, so the fit must actually run, not just apply the floor).
- Outputs, matching cellranger's real schema exactly for direct comparison: `protospacer_calls_per_cell.csv` (`cell_barcode,num_features,feature_call,num_umis`, only called cells listed), `protospacer_calls_summary.csv` (hierarchical categories: `No guide molecules` ⊆ `No confident call`, then `1 protospacer assigned` / `More than 1 protospacer assigned`, with `pct_cells`/`median_umis`/`stddev_umis`), `protospacer_umi_thresholds.csv`/`.json` (per-guide fitted threshold).
- Extend `scripts/run_flex_multi.sh` with this step after the CRISPR matrix/cells loop.

### Phase 5 — results (2026-09-14)

Implemented `src/protospacer_calling.rs` (2-component Gaussian-mixture EM,
hand-rolled — `hist`/EM helpers, not a crate — with a numerical
decision-boundary scan rather than the textbook closed-form quadratic; see
below for why) + `sbxpress protospacer-call`, and extended
`scripts/run_flex_multi.sh` with `--crispr-bam`/`--guide-features` to run
the full per-sample CRISPR chain (`molecules --tagseq` → `matrix --features`
→ `cells --tagseq-barcodes` → `protospacer-call`) after the GEX loop,
piggybacking on each sample's **Phase 3 GEX-called cell list**
(`gex_outs/filtered_barcodes.csv`, barcode column only) exactly as the plan
specifies.

**Ground-truth validation, with one limitation investigated thoroughly and
reported honestly rather than hidden:**

Pulled real per-cell CD45 UMI-count data for two real CAT06 samples
directly from their actual `sample_raw_feature_bc_matrix.h5` (Pos: 74,581
real cells; RZ211: 50,279), restricted to each sample's real called-cell
population — real ground truth, not synthetic — and checked this module's
fit against cellranger's own real `protospacer_umi_thresholds.csv`:

- **`Pos`: exact match.** Real threshold `CD45,3` (the floor binds) —
  reproduced exactly.
- **`RZ211`: does not match.** Real threshold `CD45,6` (genuinely
  adaptive, fit exceeds the floor); this implementation lands at `3`. This
  was investigated, not shrugged off:
  - A full-population fit (zeros included) is clearly wrong here: a
    dominant exact-zero spike badly distorts a continuous Gaussian mixture.
    Switched to fitting **nonzero counts only** (exact zeros are a
    capture/detection non-event, not a low-but-real "ambient" observation
    — biologically they aren't really part of either population the
    mixture is meant to separate).
  - The closed-form quadratic solve for two Gaussians' density crossing is
    numerically brittle exactly when a component's fitted variance is
    small (common here) — `1/std²` blows up its coefficients. Replaced
    with a plain numerical scan for the sign change in weighted
    log-density (20,000 steps), which has no such failure mode and is
    trivially fast at these data sizes.
  - Suspected a poor local optimum in the EM (a single percentile-based
    init can converge poorly on a harder, more continuous distribution
    like `RZ211`'s). Added multiple restarts — several percentile splits
    plus a real 1D k-means seed, mirroring how `sklearn`'s
    `GaussianMixture` avoids exactly this — and kept the
    highest-log-likelihood result.
  - **Definitively ruled out "unescaped local optimum" as the explanation**:
    evaluated this implementation's converged fit and a reference
    `sklearn.mixture.GaussianMixture` fit (fit independently in Python on
    the identical real data) against this module's own log-likelihood
    formula. This implementation's fit scores **higher** (-47,062 vs.
    sklearn's -47,572) — i.e. this *is* the actual maximum-likelihood
    optimum for a 2-Gaussian mixture on `RZ211`'s real distribution, not a
    bug or a local-optimum failure to escape. Cellranger's real,
    undocumented algorithm evidently does something beyond "MLE-fit a
    2-Gaussian mixture and take the density crossing" — matching it
    further isn't achievable from the published description alone, and
    chasing the specific number 6 without a real algorithmic lead would be
    overfitting to one data point.
  - Pinned as `real_data_rz211_cd45_threshold_regression` — an honest
    regression test of actual current behavior (asserts `3`, with the full
    investigation above in its doc comment), not a test quietly relaxed to
    pass.

**Orchestration validated end-to-end on synthetic data** (a real
plumbing/mechanics test, since the two real-data checks above already cover
the statistical core): built a synthetic CRISPR "guides.bam" (same shape
`sbxpress guides` produces — mapped, per-guide-contig, `CB`/`UB`/`GX`/`GN`
tags, no `xf`) and ran the actual extended `scripts/run_flex_multi.sh`
against it plus the existing Phase 3 GEX synthetic BAM: one cell with 5
distinct-UMI CD45 molecules (clears the floor, correctly called with
`num_umis=5`) and one with 1 (correctly *not* called, and correctly
counted under "No confident call" rather than "No guide molecules" — the
plan's hierarchical distinction between the two, exercised on a real edge
case). Caught and fixed a test-data bug this way, not a pipeline bug: `UB`
values must be real ACGT sequences for `sbx10x-molecules`'s `Dna::pack` to
accept them — a placeholder like `"UMI0000A"` (contains non-ACGT
characters) is silently treated as ineligible, zeroing out the whole run
until real-shaped UMIs were used.

10 new unit tests (7 `protospacer_calling` + reused Phase 3/4 machinery for
the orchestration check). Full suite green (188 lib tests), no clippy
warnings in new code.

**Conclusion: Phase 5's plumbing and its `Pos` ground-truth case pass.**
Its `RZ211` case is a real, investigated, and openly documented gap against
cellranger's undocumented exact algorithm — worth revisiting if cellranger's
real CRISPR-calling implementation details ever become available, but not
something to paper over with an unprincipled fit adjustment now.

### Amendment (2026-09-19) — `--tagseq` was wrong for CRISPR; known unimplemented scope

Two findings from reviewing 10x's real "Cell Ranger's CRISPR Guide Capture
Algorithm" doc against this phase's actual CAT06 output:

- **`molecules --tagseq` for the CRISPR side, as originally specified above
  and implemented in `run_flex_multi.sh`, was the wrong mode.** `--tagseq`
  exists for a genuinely different protocol — pooled, non-single-cell
  libraries where one barcode is a whole sample at tens-of-millions-of-reads
  depth (see its doc comment in `sbx10x-molecules`). Real CAT06 CRISPR Guide
  Capture data is ordinary single-cell depth (real per-cell counts: 6-7
  UMIs/cell against the `--min-umi 3` call floor above), the same order of
  magnitude as GEX, so `--tagseq`'s skip of cross-guide UMI collision
  adjudication *and* directional UMI error-correction was very likely
  hurting accuracy right at the threshold decision boundary, not just a
  naming mismatch. Fixed: `run_flex_multi.sh` now runs plain `molecules` (no
  `--tagseq`) for the CRISPR stream. `cells --tagseq-barcodes` is unrelated
  (an external cell-identity override — CRISPR piggybacks on GEX-called
  cells, matching 10x's real design: "cell calling is done using the Gene
  Expression data") and is unaffected.
- **Known unimplemented, out of scope for now**: 10x's real algorithm also
  covers "by target" grouping (multiple guide RNA constructs per target
  gene, results groupable by feature or by target) and downstream
  perturbation analysis (sSeq differential expression for Perturbation
  Efficiency / Perturbation Effects, log2-fold-change reporting). Neither
  exists in `protospacer_calling.rs` or anywhere else in sbxpress — this
  project's scope has been matching cellranger's raw counts/calls, not full
  secondary statistical analysis. Also moot for CAT06 specifically:
  `feature_reference_RZ.csv` has no `target_gene_id`/`target_gene_name`
  column, so even cellranger has nothing to group by on this real dataset.
  Revisit if a future pooled-screen dataset (multiple guides/gene, real
  perturbation-effect questions) needs it.

## Phase 6 — End-to-end verification against ground truth

- Per-phase spot checks as described above (Phase 0 coordinate harness, Phase 2 vs. `frp_gem_barcode_overlap.csv`/`cells_per_tag.json`, Phase 4 vs. feature-reference hit rates).
- **Final, per sample** (all 9: `Pos, RZ211, RZ212, RZ213NR, RZ213R, RZ213RD, RZ214, RZ215, Neg`): diff `<sample>/gex_outs/filtered_feature_bc_matrix` against `illumina_data/RZ16443_L1_multi/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix.h5` (cell count, median UMI/genes, top genes), and diff `<sample>/protospacer_calls_per_cell.csv`/`_summary.csv`/`_thresholds.csv` directly against the matching real files in `per_sample_outs/<sample>/crispr_analysis/` — call-level and threshold-level agreement, not just aggregate cell counts, since that's the level cellranger's actual algorithm operates at.

### Critical files
- `/mnt/analysis/software_shared/sbxpress/src/lib.rs`, `extract.rs`, `demux.rs`, `config.rs`, `barcode.rs`, `annotate_bam.rs`
- `/mnt/analysis/software_shared/sbxpress/src/bin/sbx10x-molecules.rs`, `cr_bam_to_mtx.rs`, `cell_calling.rs`, `sample_sheet.rs`
- New: `src/hamming_index.rs`, `src/probe_barcode.rs`, `src/probe_sample_sheet.rs`, `src/probe_demux.rs`, `src/guide_capture.rs`, `src/protospacer_calling.rs` + corresponding `src/bin/sbxpress-*.rs` thin wrappers and `Cargo.toml` `[[bin]]` entries
- Reference data (already on disk): `CAT06/resources/flex-v2-384.txt`, `CAT06/feature_reference_RZ.csv`, `CAT06/illumina_data/RZ16443_L1_multi/outs/{config.csv, multiplexing_analysis/*, per_sample_outs/*/crispr_analysis/*}`
  (full path: `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06/...`)

### Orchestration & CLI surface decisions
- Per-sample `matrix`/`cells` looping stays an **external shell script** (`scripts/run_flex_multi.sh`), consistent with the CAT02 convention (`run_cellranger_counts.sh`) and sbxpress's own thin-dispatcher philosophy — heavy O(records) aggregation logic (probe-demux, guide matching, protospacer calling) is Rust; looping over already-general Rust commands is shell.
- Probe-barcode *locating/tagging* is a flag on `extract` (per-read-local, same shape as existing CR/CB/UR/UB mechanism). Sample *assignment*, guide *matching*, and protospacer *calling* are each new subcommands (whole-BAM/whole-matrix aggregations, structurally like `molecules`/`cells` today).

### Next step (separate session)
Implementation happens in a fresh session rooted at `/mnt/analysis/software_shared/sbxpress` (its own git repo, its own project memory) rather than this CAT06 analysis session — see conversation discussion. This plan file is the handoff artifact for that session.

---

## Phase 6 — results (2026-09-14): the plan's missing prerequisite, found and built

The plan above assumes the existing GEX pipeline (`extract` → `minimap2` →
`annotate` → `molecules`) already produces a real per-sample GEX BAM for
Flex data — Phase 6 was written as "just diff the outputs." **That
assumption is wrong**, discovered empirically in this session's earlier
full-chain `--no-gzip` profiling test: real CAT06 Flex R2 aligned against
the whole genome (splice-aware, real GRCh38 index) found almost nothing —
59 real molecules out of 1.65M records. Flex doesn't sequence free cDNA;
each read's payload is one of a fixed, closed set of ~54k known 50bp probe
sequences (`Chromium_Human_Transcriptome_Probe_Set_v2.0.0...csv`, already
referenced in `config.csv`'s `[gene-expression]` section but never
consumed by anything in this plan). Standard genome alignment was never
going to produce real Flex GEX counts — a probe-set matcher was a missing
prerequisite, not an optional refinement.

### Confirmed the read layout empirically before building anything

Same rigor as every other coordinate in this plan (`examples/probe_seq_locate.rs`,
not checked in as a permanent tool — a one-off diagnostic): scanned real
extracted R2 reads for **exact** 50bp substring matches against the real
probe set. Forward orientation, dominant offset **69** (jitter into the
low 60s/high 60s — the same `dt_run_end`-measurement-jitter shape as the
probe barcode and the CRISPR anchor, at yet another distinct offset, since
GEX and CRISPR are physically different library constructs). **42% of all
reads (mixed GEX+CRISPR+noise) matched exactly**, zero mismatch tolerance
— and the genes hit were immediately recognizable real biology (TMSB4X,
B2M, ACTB, VIM, IL32, CD74 — housekeeping + immune markers, consistent
with a PBMC-like sample), not noise.

### Built `src/probe_set.rs` + `sbxpress probeset-align`

`ProbeSetReference::load` (real schema: `gene_id,probe_seq,probe_id,
included,region,gene_name`, `#key=value` metadata header skipped,
`included=FALSE` rows excluded — 10x's own curation of probes too
cross-reactive to count on). `ProbeMatcher` is architecturally the
probe-barcode's own `locate()` pattern reused a third time (windowed,
Hamming-tolerant lookup on the shared `HammingIndex`, centered on the
empirically confirmed offset 69) rather than `feature_capture.rs`'s
anchor+capture pattern — there's no anchor here, just a direct fixed-length
lookup. `--max-mm` defaults to 1, not more: `HammingIndex` precomputes
every ≤max_mm neighbor per key, and neighbor count grows as `C(50,k)·3^k`
— 1 is ~8M entries (fast, ~1.3GB peak RSS observed); 2 would be ~590M
(tens of GB), impractical to build. Output BAM is shaped exactly like
`feature_capture.rs`'s (synthetic per-gene contigs, `GX`/`GN` tags, `CB`/
`UB`/`PB` carried through, unmapped-but-tagged on no match) specifically
so it feeds `molecules`/`matrix`/`cells`/`probe-demux` completely
unchanged — no new consumer-side code needed anywhere else in the plan.
10 new unit tests (reference loading incl. `included=FALSE` filtering,
matcher exact/mismatch/windowed/multi-gene behavior, manifest
create-or-verify, full `run()` end-to-end BAM assertions). 198 lib tests
total, no clippy warnings.

### Real-data validation: production run, then the full real chain

Ran the actual `sbxpress probeset-align` binary (not the diagnostic)
against the real probe set and the 800k-read CB-validation sample:
**51.5% hit rate** (up from the 42% exact-only baseline, confirming
`--max-mm 1` is doing real work), **10,487 distinct genes hit** — a
sensible count for real scRNA-seq — same top genes as the diagnostic.

Chained the actual, unmodified rest of the pipeline on top, for real —
this is the first time in this plan a Flex sample has gone all the way
from raw reads to a real gene-expression matrix:
`probeset-align` → `sbx10x-molecules` (normal UMI correction, **not**
`--tagseq` — Flex GEX is real per-cell single-cell depth, not the
pooled/aggregate-barcode depth `--tagseq` exists for; that flag's use here
stays CRISPR-specific) → `sbxpress probe-demux` against the real
`config.csv` → `matrix`/`cells` for sample `Pos`.

- **Sample assignment** (all 9 real samples, gene-expression-based this
  time, not just probe-barcode counting): `Pos` highest (2,450 of 9,624
  total assigned cells at this reduced depth), `RZ214` second — matches
  real cellranger's own ranking's top two exactly. The `RZ213` trio
  remained the clear minority group, as in every earlier phase's check.
  Not a perfect rank match further down (real ground truth: `Pos > RZ214 >
  Neg > RZ212 > RZ215 > RZ211 > ...`; ours: `Pos > RZ214 > RZ212 > RZ211 >
  RZ215 > Neg > ...`) — expected at 800k reads (~0.1% of real depth), not
  investigated further at this scale.
- **Gene expression, the real payoff**: compared `Pos`'s top-15 expressed
  genes against real cellranger's actual `sample_filtered_feature_bc_matrix.h5`
  for `Pos` (full real depth). **9 of our top 10 genes are in cellranger's
  real top 10**, in nearly the same order (`TMSB4X` #1 and `ACTB` #2 in
  both) — from an 800k-read sample against cellranger's full-depth result.
  This is strong, real evidence the whole new alignment strategy produces
  biologically correct results, not an artifact of small-sample noise
  happening to look right.

### What Phase 6 still needs for the plan's original, literal bar

The plan's own validate step asks for all 9 samples, full real depth,
cell-count/matrix-level diffs, and (once real CRISPR-library sample data
exists) real `protospacer_calls_*` comparisons. What's been done here is
a real, working, validated pipeline and a strong small-sample proof that
it's biologically correct — not yet the full 754-shard, full-depth,
all-9-sample run this phase's literal bar describes. That's a
multi-hour-plus compute job on the infrastructure now built (real per-read
throughput is established from the profiling work — `--no-gzip` +
`probeset-align`'s ~40k reads/sec single-threaded, parallelizable), not
something to run to completion inline here. Natural next step when
pursuing full validation: run `scripts/run_flex_multi.sh`-style
orchestration (extended to call `probeset-align` instead of assuming a
pre-existing GEX BAM) across the full raw dataset.

---

## Amendment (2026-09-14) — review notes before implementation

Reviewed the plan against the current source (0.3.8, post-tagseq/`--trim-r2-polya` release) and against two follow-up questions.

### 1. Impact on `--tagseq` and the existing GEX model

Confirmed **no functional change** to existing GEX or `--tagseq` behavior anywhere in Phases 2–6 — every one of those phases is purely additive (new files, new subcommands, new optional flags gated off by default) and reuses `matrix --cell-barcodes`, `molecules --tagseq`, and `cells --tagseq-barcodes` as-is against new BAM inputs, without touching their implementations.

Two places in **Phase 1** are structural refactors of currently-live, shared code rather than pure additions, and are the actual risk surface:
- `demux.rs`'s `Slot`/`insert_neighbors`/`enumerate` get pulled out into a new `src/hamming_index.rs`, with `demux.rs` becoming a thin wrapper. `demux.rs` is on the default path for every `extract --sample-sheet` run (GEX included) and currently has 9 unit tests — the refactor must keep those green and needs a before/after diff on a real GEX sample's demux metrics, not just "tests pass," before merging.
- `lib.rs::Shared` gains a new `probe_wl: Option<&ProbeWhitelist>` field threaded through `process_one`. The plan states this is zero-cost when unset, which is the correct design, but `lib.rs` is the core extract loop every run goes through, so Phase 1's acceptance criteria should explicitly include: existing `lib.rs`/`extract.rs` test suite green, and a byte-identical (or metrics-identical) `extract` output on a real GEX run with the new flags absent.

Everything downstream of Phase 1 (probe-demux, guide matching, protospacer calling, the orchestration script) operates on new BAMs/matrices and cannot regress GEX or tagseq even in principle — they're new consumers, not new producers into the existing path.

### 2. Generalize CRISPR / Feature Barcoding beyond Flex

10x's "Feature Barcoding" read structure (`<anchor>(BC)` + a `feature_reference.csv` with `id,name,read,pattern,sequence,feature_type`) is not Flex-specific — it's the same mechanism behind CRISPR Guide Capture, Antibody Capture (CITE-seq), and Cell Multiplexing (CMO hashing), and 10x runs all three on ordinary 3'/5' poly-dT GEX libraries with **no probe barcode at all**, not just on Flex. As drafted, the plan is CRISPR+Flex-shaped by construction (naming, and an implicit dependency on `PB` in Phase 4's tag list). Two changes keep the underlying capability general, without changing phase scope or the CAT06 delivery goal:

- **Make `PB` strictly optional in Phase 4.** The guide/feature matcher should parse whatever of `CR/CB/UR/UB/PR/PB` is present in the FASTQ comment and carry it through unchanged — `PB` will simply be absent on a non-Flex (plain GEX) Feature Barcoding run, and the matcher must not assume it exists.
- **Generalize naming and dispatch in Phase 4, keep Phase 5 CRISPR-specific.** Rename `guide_capture.rs` → `feature_capture.rs` and `GuideMatcher` → `FeatureMatcher`, with matching driven by the reference CSV's `feature_type` column rather than hardcoded "guide" naming; `sbxpress guides` can stay as the CRISPR-flavored CLI entry point (or become `sbxpress features --feature-type "CRISPR Guide Capture"`) built on the general matcher. Phase 5's Gaussian-mixture protospacer calling is genuinely CRISPR-specific (that's cellranger's actual published algorithm for guides — Antibody/Multiplexing Capture use different downstream calling, e.g. normalization or CMO-based sample assignment) and should stay its own module; only the matching/tagging layer underneath (Phase 1's `HammingIndex` + Phase 4's anchor+capture matcher) needs to be assay/feature-type-agnostic.
- Net effect: Phase 1 and Phase 4's matching core become a general "Feature Barcoding" building block usable for any 10x protocol that uses this read structure (Flex or not, CRISPR or not); Flex probe-barcode demux and CRISPR guide calling remain the two concrete consumers this plan delivers first.
