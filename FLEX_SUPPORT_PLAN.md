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

## Phase 4 — `guide_capture.rs` + `guides` subcommand (CRISPR read → guide assignment)

- New `src/guide_capture.rs`: `FeatureReference::load(csv)` parsing the real `id,name,read,pattern,sequence,feature_type` schema (confirmed against `CAT06/feature_reference_RZ.csv`); parse `pattern` as `<anchor>(BC)` (confirmed layout; flag `(BC)<anchor>` as an explicit unsupported-format error rather than silently mis-parsing). `GuideMatcher` built on Phase 1's shared `HammingIndex`: one index for anchors (windowed, Hamming-tolerant, mirrors `demux::assign_windowed`), then match the preceding `bc_len` bases against a second small index of guides' `sequence` values.
- New `src/bin/sbxpress-guides.rs` (new `[[bin]]`, new dispatcher arm). CLI: `sbxpress guides <R2_FASTQ> --feature-reference <CSV> -o <BAM> [--features-out <TSV>] [--anchor-max-mm 2] [--capture-max-mm 1] [--anchor-window 5]`.
- Input is the `--minimap2`-preset R2 FASTQ from `extract` on the CRISPR-library i7/i5 stream (no `extract.rs` changes needed — its existing dT-run-anchored bracketing already captures the anchor+capture region as a side effect). Parses `CR/CB/UR/UB/PR/PB` back out of the FASTQ comment, matches guides, emits a **mapped** BAM against a synthetic per-guide contig header (tid=guide index, trivial CIGAR, mapq 60 on match; unmapped — correctly skipped by `molecules` — on no match), tags `GX`/`GN`=guide id/name plus carried `CB`/`UB`/`PB`. `--features-out` writes a canonical `id\tname` TSV (mirrors `annotate_bam.rs`'s `--genes-mtx` create-or-verify pattern) for `matrix --features`.
- **Validate**: anchor/capture hit rate against `feature_reference_RZ.csv` on real raw reads.

## Phase 5 — CRISPR counting + **new protospacer-calling module** (the published-algorithm gap)

Per sample: `sbxpress molecules guides.bam --tagseq -o corrected.bam` → `sbxpress matrix corrected.bam --cell-barcodes <sample>.gex_cells.txt --features guide_features.tsv -o <sample>/crispr_raw` (use the sample's **Phase 3 GEX-called cell list**, not the pre-cell-calling probe-demux list — CRISPR identity piggybacks on GEX-defined cell identity, matching real 10x biology) → `sbxpress cells <sample>/crispr_raw/... --tagseq-barcodes <sample>.gex_cells.txt -o <sample>/crispr_outs`. This gives a raw per-cell-per-guide UMI matrix, but **stops short of cellranger's actual call** — do not treat any nonzero count as "detected."

- **New `src/protospacer_calling.rs` + `sbxpress protospacer-call` subcommand**: input is the per-sample filtered CRISPR matrix + feature reference. For each guide: fit a 2-component Gaussian mixture on `log1p(UMI count)` across all called cells (background/ambient component vs. signal component — standard EM for a 2-component 1D GMM), derive the per-guide UMI threshold from the fitted decision boundary, and **enforce `max(fitted_threshold, --min-umi)` unconditionally** (default `--min-umi 3`, matching cellranger's floor — confirmed genuinely adaptive above/below this floor across real samples, so the fit must actually run, not just apply the floor).
- Outputs, matching cellranger's real schema exactly for direct comparison: `protospacer_calls_per_cell.csv` (`cell_barcode,num_features,feature_call,num_umis`, only called cells listed), `protospacer_calls_summary.csv` (hierarchical categories: `No guide molecules` ⊆ `No confident call`, then `1 protospacer assigned` / `More than 1 protospacer assigned`, with `pct_cells`/`median_umis`/`stddev_umis`), `protospacer_umi_thresholds.csv`/`.json` (per-guide fitted threshold).
- Extend `scripts/run_flex_multi.sh` with this step after the CRISPR matrix/cells loop.

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
