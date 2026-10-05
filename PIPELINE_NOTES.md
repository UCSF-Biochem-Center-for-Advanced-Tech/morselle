# Pipeline resource & architecture notes

Real, measured operational findings from running the full CAT06 pipeline
(extract → probeset-align/guides → samtools cat → samtools sort/index →
sbx10x-molecules → probe-demux → matrix → cells → protospacer-call) at true
production scale — 754 raw shards, ~32.6B raw reads, on a 512-core / 1.5TB
RAM box (NFS-mounted storage at `/mnt/analysis`, one shared mount, so its
aggregate I/O bandwidth is shared across every concurrent job on the box,
not per-process).

Purpose: capture what actually happened, stage by stage, so a future
**unified Rust pipeline** (replacing today's shell scripts gluing separate
binaries + `samtools` together) can manage CPU and memory deliberately
instead of rediscovering these constraints the hard way again. See
`/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06/read_depth_comparison_plan.md`
for the specific run this was learned on; this doc is the general,
tool-level takeaway, not CAT06-specific.

## Per-stage resource profile (real measurements)

| Stage | CPU behavior | Memory | I/O | Notes |
|---|---|---|---|---|
| `extract` | Scales genuinely well with `--threads`/`--split` | Modest, bounded per-thread buffers | Write-heavy (gzip output) | 32.6B reads in 5h40m at `--threads 128 --split 128` (~1.6M reads/sec aggregate). **`--split` must round-robin, not tie chunk assignment to thread identity** — a real bug found and fixed this session (`5f4f2c8`): thread-affinity chunking let 2 of 128 output chunks balloon to ~4.4x the median size, because real per-shard yield varies and rayon's work-stealing doesn't guarantee even file counts per thread. |
| `probeset-align` / `guides` | **Single-threaded per process, no internal parallelism** | Modest | Read+write moderate | Must be fanned out externally (`xargs -P<cores>`, one process per `--split` chunk) to use multiple cores at all. Real per-chunk timing: ~977s for a 58M-read GEX chunk, ~172s for a 14.6M-read CRISPR chunk (single-threaded each). |
| `samtools cat` (merge) | Low CPU (no recompression, same-codec concatenation) | Low | **Pure I/O, can be fast**: observed 150–400MB/s real sustained throughput on this box's NFS mount | Only concatenates — does **not** sort. Downstream indexed access needs an explicit sort step (see next row); don't assume `cat` output is index-able. |
| `samtools sort` | `-@` mainly parallelizes *compression*, not the core in-memory sort-comparison step — real observed pattern: bursty CPU (idle while filling a buffer, spikes to several cores while compressing+flushing one chunk), not sustained at the requested thread count | **`-m` (memory per thread) is the real lever**, not `-@`, once you have RAM to spare | Directly memory-gated: bigger `-m` → fewer, larger chunks → far less repeated read/sort/compress/flush overhead | Real before/after on this box: `-m 2G` → ~70MB/s effective read throughput (many small ~37.5GB chunks, one every ~10min); `-m 8G` → ~254MB/s (3.6x) on the identical file. **Always check `free -h` before raising `-m`** — this is a real, direct memory commitment, not just a suggestion. **Follow-up, real numbers from watching the actual `-m 8G` run on `TSeq100_L1_GEX.bam` (703GB) to completion of its first two chunks:** an early single-snapshot read (~26min in) caught `tmp.0000.bam` at ~44.2GB and was wrongly logged here as a *finished* chunk size — it was actually a **mid-flush snapshot**, not the final size; that "~18% bigger than `-m 2G`" conclusion was an artifact of reading too early and is **corrected** below. Real finished-chunk sizes, read back from disk once each `tmp.NNNN.bam` stopped growing: `tmp.0000.bam` = **150.7GB**, `tmp.0001.bam` = **150.6GB** — roughly **4x** `-m 2G`'s ~37.5GB chunks, in line with the 4x memory-per-thread increase after all. Revised chunk-count estimate for this 703GB file: **~5 chunks at `-m 8G`** (vs. ~19 at `-m 2G`), not the ~16 first guessed from the bad early sample. `/proc/<pid>/io` cross-check at 2h07m elapsed: 465GB read of 703GB (66%) at an average 62.7MB/s over the full elapsed window (lower than the ~254MB/s *peak* measured earlier, because the average includes idle fill-phase time, not just active-read bursts) — real end-to-end lesson: **don't trust a single early snapshot for chunk-size or throughput conclusions; only read finished-file sizes and full-window `/proc/<pid>/io` deltas.** Process RSS held steady at **~1.1TB** across both the early and later checks (tracks the naive `-m 8G × ~128 threads ≈ ~1TB` ceiling), with system-wide `free -h` "available" holding ~403-404GB throughout — a real, direct memory commitment worth checking against whatever else is running concurrently before choosing `-m`, not a free lunch, but stable rather than climbing further over time. **Monitoring note:** a single `/proc/<pid>/io`-delta snapshot reading ~0.0MB/s does **not** mean the process stalled — external-merge sort alternates fill (reading, low measured I/O) and flush (writing a completed chunk) phases, so a 0MB/s *read*-throughput sample during a flush window is expected; only sustained zero across multiple snapshots spanning both phases would indicate a real stall. |
| `samtools index` | Low | Low | Read-heavy, fast | **Requires coordinate-sorted input** — fails outright ("NO_COOR reads not in a single block") on `samtools cat` output. Sort first, every time. |
| `sbx10x-molecules` | Default `--evidence-workers 1` is **fully serial** and dramatically slower at billion-read scale despite a real parallel path (`--evidence-workers N`, commits `3b7366c`/`9afc664`) existing in the binary since 2026-09-02 — that path additionally **requires the input BAM to be indexed**, so it only activates after the sort+index steps above. | **By far the dominant memory consumer in the whole pipeline.** Real observed: **~976GB RSS** for one ~9B-read / ~4.02B-raw-UMI-node library (`TSeq100_L1_GEX`), serial mode. Root cause (verified by reading the code, not guessed): `evidence` (`HashMap<RawKey,Evidence>`), `roots` (`HashMap<RawKey,RawKey>`), and `decisions` (`HashMap<RawKey,Decision>`) are all genuinely alive *simultaneously* at the point `decisions` gets built — each ~4B entries, ~24–48 bytes/entry raw plus real hashbrown overhead. No cheap safe reduction found; would need real restructuring (e.g. streaming `decisions` construction while draining `roots`/`evidence`) to meaningfully cut peak memory. | Read-heavy (evidence collection), write-heavy (final BAM rewrite) | **Two real libraries' worth of this can't run concurrently** on a 1.5TB box — confirmed L1 alone used ~976GB, leaving only ~552GB free, not enough headroom for a second comparably-sized library's molecules pass at the same time. This is the single biggest scheduling constraint discovered this session. **Answered, and the answer is bad**: `--evidence-workers 64` on this same `TSeq100_L1_GEX` library was real-world **OOM-killed** by the kernel on 2026-09-16 — confirmed via `sudo dmesg`: `Out of memory: Killed process 503511 (sbx10x-molecule) total-vm:1676222284kB, anon-rss:1561676212kB` — i.e. **~1.49TB anon-RSS**, essentially the whole 1.5TB box, right at the `umi_correct` → `molecule_aggregate` phase transition (union-find fully resolved, `roots=4,078,415,170`, aggregation just starting). That's ~525GB **more** than the ~976GB serial-mode figure above, not less — parallel mode's per-worker local tables (each of the 64 evidence-workers builds its own region-local table before the global remap, same pattern used for genes/samples) add real duplicated overhead on top of the same `evidence`/`roots`/`decisions` structures serial mode already needed, rather than replacing them. **Practical implication: `--evidence-workers N>1` is not a safe default at this library's scale on this box** — it trades a real, large memory-safety margin for `umi_correct`-phase wall-clock speed, and at N=64 that trade goes negative (crash, not just "faster but chunkier"). Serial mode's own true peak was never fully confirmed either (that run was killed intentionally, mid-run, for the correctness bug — not because it hit a wall) — so "~976GB" should be read as "highest observed," not "confirmed peak," until a serial run is actually let run to completion. Whether some smaller worker count (e.g. 8–16) keeps most of the speed benefit while staying under the serial figure is untested. |
| `probe-demux` / `matrix` / `cells` / `protospacer-call` | Single-threaded each | Low | Low | Fast, lightweight relative to everything above — never the bottleneck. |

## Concurrency lessons (real, not theoretical)

- **Check `free -h` before scheduling anything memory-heavy concurrently.**
  This session nearly ran a second ~976GB-class `molecules` job while the
  first was still using ~976GB, on a 1.5TB box — caught by checking real
  numbers first, not assumed safe from thread/CPU headroom alone. CPU
  headroom (this box has 512 cores, most runs use ≤128) does **not** imply
  memory headroom — they're independent budgets and both matter.
- **`extract`/`align`-stage jobs and `molecules`-stage jobs have very
  different resource shapes** — the former is CPU+I/O bound with modest
  memory, the latter is memory-dominant. A resource-aware scheduler should
  track these as genuinely different budget categories, not just "cores
  busy."
- **NFS I/O bandwidth is a shared, single pool** across every concurrent job
  on this box (one mount, `/mnt/analysis`) — real throughput observed
  ranged from ~70MB/s (many-small-flushes pattern) to ~400MB/s (large
  sequential transfer pattern) *for the same underlying storage*, purely
  based on the I/O access pattern of the job doing it. A future pipeline
  should treat "network storage bandwidth" as a real, finite, shared
  resource dimension too, not assume it scales with `--threads`.
- **A flag's name doesn't guarantee what it parallelizes.** `-@`/`--threads`
  on `samtools sort` sounds like it should scale the whole operation
  linearly; in practice it mainly scales compression, and the real lever
  for this specific bottleneck was `-m` (memory). Verify real behavior
  (measure `/proc/<pid>/io` deltas, real elapsed time, `ps` thread/CPU%)
  rather than trust a flag's apparent semantics before committing to a
  large run — this session's whole "why is samtools only using 4 CPUs"
  investigation is the concrete example.

## Design considerations for a future unified Rust pipeline

Once the current shell-script-orchestrated version is validated and
"solid" (per-checkpoint molecule/matrix comparisons in
`read_depth_comparison_plan.md` producing real, trusted results), folding
`extract → align → sort/index → molecules → demux → matrix → cells →
protospacer-call` into one Rust-native pipeline should deliberately design
for the constraints found above, not just port the shell logic:

1. **A real resource scheduler, not just a process-pool.** Stages should
   declare an approximate resource shape (CPU-bound-and-parallel,
   memory-dominant-mostly-serial, I/O-bound, etc.), and the scheduler
   should refuse to run two memory-dominant stages concurrently unless
   there's real, checked headroom (query actual available memory, don't
   assume from a static config) — exactly the collision this session found
   and avoided by hand.
2. **Native sort+index, not a shell-out.** If `molecules`' parallel evidence
   path needs indexed BAM access, the pipeline should own that
   sort/index step directly (with a memory-aware chunk size chosen from
   real available RAM at run time, not a hardcoded default) rather than
   relying on a separately-invoked `samtools sort -m <guessed value>`.
3. **Bound `molecules`' peak memory by design, not by hoping it fits — real
   design confirmed 2026-09-16, after a real OOM kill.** `--evidence-workers
   64` on `TSeq100_L1_GEX` was killed by the kernel: `sudo dmesg` confirmed
   `anon-rss:1561676212kB` (~1.49TB) at the `umi_correct`→
   `molecule_aggregate` boundary — ~525GB **more** than serial mode's ~976GB
   high-water mark, because parallel evidence collection's per-worker local
   tables add real duplicated overhead on top of the same global structures
   serial mode already needed, not instead of them. Reading the actual code
   (`build_decisions()`) to find the *right* partition axis, not just any
   axis:
   - **Partitioning by genome/synthetic-contig position (one gene per
     contig, matching this project's pseudo-reference) does NOT work.**
     `groups` (UMI correction) is keyed `(cell, sample, gene)`, so it's
     gene-scoped and *could* stream per-contig — but
     `gene_collision_index`/`gene_collision_resolve` (deciding which gene
     wins when the same UMI shows up under more than one gene for a cell) is
     keyed `(cell, sample, umi)` **with no gene in the key**, and two
     genuinely different genes for the same cell+UMI can sit on two
     completely different, far-apart contigs. Streaming gene-by-gene would
     silently miss real cross-gene UMI collisions.
   - **The right axis is `(cell, sample)`** — both real grouping keys in the
     algorithm share this as their outer key, so partitioning evidence by a
     hash of `(cell, sample)` into K bounded shards, then running the full
     `groups → correct → molecules → collisions → decisions → write`
     sequence **per shard**, discarding each shard's memory before the next,
     is safe by construction: no correction or collision-resolution step
     ever needs to see two different `(cell, sample)` shards at once. Peak
     memory becomes ~`total/K`, a real tunable knob (matching `extract`'s
     `--split` philosophy), independent of total library size. **This is the
     established, proven pattern** — it's how Cell Ranger (Martian pipeline
     chunks/joins by cell-barcode range), STARsolo (per-cell-barcode
     structures released once a cell's counting finishes), and salmon
     alevin/alevin-fry (`collate`-by-barcode then per-cell `quant`,
     built specifically because plain `alevin` hit this same memory wall at
     scale) all solve this same problem. No BAM re-sort is required — bucket
     assignment can happen during the existing per-region parallel read
     pass.
   - **A smaller, lower-risk quick win, not yet applied**: `evidence`'s
     *values* (the heavy `Evidence` struct) are only read up through the
     point `molecules` is built; the final `decisions`-materialization loop
     only needs `evidence`'s *keys*, which are identical to `roots`' keys
     (a much lighter value type). `evidence` is currently passed by
     reference and kept alive by the caller for the whole function, so it
     never gets freed despite being dead weight for the back half of
     `build_decisions` (`gene_collision_index` through
     `decision_materialize`). Passing `evidence` by value and dropping it
     right after `molecules` is built (switching the final loop to iterate
     `roots.keys()` instead) would free one whole ~4B-entry heavyweight map
     for that back half — real, bounded-risk, and worth doing regardless of
     whether/when the full cell-sharded rewrite happens.
4. **Round-robin output sharding by construction, not thread-affinity.**
   The `--split` bug (fixed as a bolt-on `AtomicUsize` counter this
   session) should just be the native design in a unified pipeline: never
   tie an output shard's identity to which worker happened to process a
   given unit of input.
5. **Real telemetry, not assumed flag semantics.** Emit real peak-RSS,
   real per-stage wall time, and real I/O throughput per run (not just
   record counts) so future capacity planning is based on measurements
   like the ones in this doc, not re-derived by hand each time from `ps`
   and `/proc/<pid>/io` the way this session had to.

## Open items to fill in once observed

- ~~Real peak memory for `sbx10x-molecules --evidence-workers 64` (parallel
  mode) vs. the ~976GB serial-mode baseline~~ — **answered, 2026-09-16**: it
  OOM-killed at ~1.49TB, ~525GB worse than serial, not better. See the
  cell/sample-sharding design note above. Serial mode's own true peak is
  still unconfirmed (that run was killed intentionally, mid-run, for the
  correctness fix, not because it hit a wall).
- Real wall-clock improvement from `-m 8G` sort + `--evidence-workers 64`
  molecules, full library, once L1 completes under the new settings —
  **superseded**: the plan moved away from running full-scale (32.6B-read,
  ~11B-reads-per-library) `molecules` jobs at all — see
  `read_depth_comparison_plan.md`'s CP1–CP3 checkpoints (each ≤~8.1B raw
  reads total across all 4 libraries, safely under a working ~10B ceiling
  for what this comparison actually needs) as the real path forward instead.
- ~~Not yet implemented, real next step when there's time: the
  `(cell, sample)`-hash-sharded `molecules` rewrite described above.~~ —
  **implemented and in real use, 2026-09-17**: `--shards 8
  --shard-concurrency 8` (see `sbx10x-molecules --help`) is exactly this
  design. **Correction, 2026-09-19**: CP1's *L1* run did **not** actually use
  `--shard-concurrency` — it started before that flag existed in
  `read_depth_comparison.sh` (see the dated evidence in the section below) —
  only CP1's *L2* run did. So "CP1's checkpoints actually run this" was only
  true from L2 onward, not for the L1 numbers already on record. It solved
  the memory problem; see "`sbx10x-molecules` real wall-clock vs. Cell
  Ranger, 2026-09-19" below for the current, corrected read on real *time*
  cost.

## `sbx10x-molecules` real wall-clock vs. Cell Ranger, 2026-09-19

In-progress investigation, updated same-day as first written. Written from
CAT06's real CP1 run (`read_depth_comparison_plan.md`, restarted
2026-09-18 with `sbxpress` 0.4.0's composite-cell-barcode fix); see that doc
for full CAT06 context if needed, but this specific question stands on its
own.

**Real numbers, all four measured directly, not estimated:**

| | sbxpress `molecules` (CP1 checkpoint subset) | Cell Ranger `ALIGN_AND_COUNT` (real Illumina full-depth run) |
|---|---|---|
| L1 | **33,427.9s (9.28h)** — `cp1/molecules/L1_gex.metrics.json`'s `elapsed_secs` | **12,569.6s (3.49h), 207.7 core-hours** — `illumina_data/RZ16443_L1_multi_full/_perf` |
| L2 | **10,178.0s (2.83h)** — `cp1/molecules/L2_gex.metrics.json`'s `elapsed_secs` | **1,031.3s (17.2min), 30.0 core-hours** — `illumina_data/RZ16443_L2_multi_full/_perf` |

(Cell Ranger figures pulled by flattening `forks[].fork_stats` in each
pipestance's `_perf` JSON — a flat list of 432 stage entries, not a nested
tree — and matching `fqname` ending in
`..._MATRIX_COMPUTER.ALIGN_AND_COUNT`, the real directly-comparable
alignment+UMI-counting combined stage.)

**The sbxpress L1-vs-L2 3.3x gap is explained, and it isn't an algorithm
difference — corrects the "unexplained" claim this section originally
made.** Real, dated evidence:

- `--shard-concurrency` (the flag that runs shards' evidence/correct/
  aggregate concurrently instead of one at a time) was committed at
  **2026-09-17 08:32:06** (`git log` on `src/bin/sbx10x-molecules.rs`,
  commit `e4dd33e`).
- `read_depth_comparison.sh` was edited to actually pass
  `--shard-concurrency 8` on the command line at **2026-09-17 19:04**
  (file mtime).
- `L1_gex.bam`/`L1_gex.metrics.json` finished at **2026-09-17 09:18** — only
  46 minutes after the flag was *committed*, for a run that took 9.28h, so
  it must have *started* around 2026-09-17 00:00, well before either the
  commit or the script edit existed. **L1 ran with shards processed fully
  sequentially, not concurrently**, confirmed directly from
  `cp1/logs/molecules_L1.log`: every `phase=shard_quant start shard=N/8` is
  immediately followed by that same shard's `complete` before the next
  `start` appears anywhere in the file — zero interleaving.
- `L2_gex.bam`/`L2_gex.metrics.json` finished at **2026-09-18 00:56**, after
  the 19:04 script edit. `cp1/logs/molecules_L2.log` confirms real 8-way
  concurrency: all 8 `phase=evidence_collect mode=serial start` lines fire
  back-to-back at the top of the log (lines 68–75), before any shard's
  `complete`.
- So the 3.3x speedup is (at least mostly) "L2 got real 8-way parallelism
  during evidence-collect/correct/aggregate and L1 didn't" — not a per-library
  skew or algorithm difference. **This means L1's 33,427.9s figure is the
  unparallelized baseline, not a fair "current sbxpress" number** — L2's
  10,178.0s is the one actually representative of what `--shards 8
  --shard-concurrency 8` produces today.

**Two real caveats that make the Cell Ranger comparison less clean than
this section originally implied — flagged, not yet resolved:**

1. **Depth mismatch, both L1 and L2.** Cell Ranger's figures are from
   real **full-depth** Illumina runs (`illumina_data/*/outs/qc_library_metrics.csv`,
   `Number of reads` for the GEX fastq ID: L1 = **2,963,748,016**, L2 =
   **2,187,154,595**). sbxpress's CP1 figures are from a **deliberately
   reduced checkpoint subset** (`records` in the `molecules` metrics.json:
   L1 = 1,922,327,467, L2 = 1,637,555,509 — by design, CP1 targets ≤~8.1B
   raw reads total across all 4 libraries, not full depth). Neither
   comparison above is apples-to-apples on input scale; this was true of
   the original L1-only comparison too, just not flagged until now.
2. **Cell Ranger's own L1 vs. L2 numbers don't internally hang together.**
   L1 has only 1.35x more full-depth reads than L2 (2.96B vs. 2.19B), but
   `ALIGN_AND_COUNT` walltime differs by **12.2x** (12,569.6s vs. 1,031.3s)
   and core-hours by 6.9x (207.7 vs. 30.0) — average concurrency was
   actually *higher* for the smaller L2 job (~104.8 cores vs. ~59.5 for
   L1). A 35% read-count difference cannot explain a 12x runtime
   difference by itself. Real cause unknown: could be a genuine property
   of these two libraries (very different cell count or saturation
   profile), or an artifact of Illumina-cluster scheduling/contention
   unrelated to input size. **Until this is understood, neither Cell
   Ranger figure alone should be treated as a stable per-read baseline**,
   and a naive "sbxpress L2 / Cell Ranger L2" ratio (≈9.9x) would be at
   least as likely to reflect this Cell-Ranger-side anomaly as a real
   sbxpress regression — not being quoted as a headline number for that
   reason.

**Net, honest read as of 2026-09-19:** sbxpress's `molecules` stage is
slower than Cell Ranger's combined align+count stage by somewhere in a wide
band (the L1 comparison says ~2.7x; a naive L2 comparison says ~9.9x), but
neither multiplier should be trusted as precise given the depth mismatch
and the unexplained Cell-Ranger-side L1/L2 inconsistency above. The
*qualitative* gap (sbxpress markedly slower) is real and worth continuing
to investigate at the code level (see the shard-granularity and
serial-evidence-collection findings already in this doc's design-
considerations section); the *precise multiplier* is not yet trustworthy
enough to quote on its own.

**Not yet checked, real starting points for continuing this investigation:**

- Where real wall-clock actually goes across `molecules`' own phases
  (`evidence` collection / `umi_correct` / `molecule_aggregate` / final BAM
  write) — **instrumented 2026-09-20.** Every molecule stderr phase/progress
  line now keeps its existing fields and appends `elapsed_secs=<seconds>`
  from one shared run-relative clock, including interleaved concurrent-shard
  work. A fresh, uncontended L2 run is now needed to capture and compare
  those phase durations; old CP1 logs remain order-only because they predate
  this instrumentation.
- Whether `--shards 8` is actually a good choice for wall-clock (it was
  chosen for a real memory constraint — see the OOM finding above — not
  profiled for speed); a real per-shard timing breakdown (once phase
  timestamps exist) would show whether shards are balanced or whether a
  few stragglers dominate (the same kind of imbalance `--split`'s
  round-robin bug caused elsewhere in this pipeline, real precedent worth
  checking for here too). **Partially answered, 2026-09-19**: tested
  `--shards 16 --shard-concurrency 12 --threads 8` on L2 (2x today's shard
  count, 1.5x today's concurrency) — real result **9,430.0s vs. the
  baseline's 10,178.0s, only ~7.3% faster (1.08x)**
  (`molecules/L2_gex_shardtest_s16c12.bam.metrics.json`), far short of what
  "8 of 512 cores busy" would predict if shard-concurrency alone were the
  bottleneck. This identified `collect_evidence`'s single-threaded-per-shard
  hot loop as the real remaining constraint. **Implemented 2026-09-20:**
  `sbx10x-molecules --shard-split N` now round-robins each logical shard's
  collated records into N physical sub-chunks, scans those sub-chunks in
  parallel, and uses the existing remap/merge path to produce the identical
  logical evidence map. This is deliberately independent of
  `--shards`/`--shard-concurrency`: the maximum evidence-scan width is their
  product, so it must be sized against real CPU and memory headroom.
  `shard_split=1` preserves the former one-file-per-shard layout exactly.
  Unit coverage verifies both balanced routing (including concurrent
  `collate` readers) and identical evidence/decisions for split versus
  unsplit input. What remains is an uncontended production-scale timing run
  to quantify the improvement; the one real-scale run so far was confounded
  by an unrelated concurrent job on the shared NFS mount. Minor caveat on
  the original shard-count comparison: the test run's metrics show a
  `no_sample_records` counter the baseline did not (an unrelated codebase
  evolution), accounting for its small `eligible_records` difference but
  not the timing gap.
- Why Cell Ranger's own L1 vs. L2 `ALIGN_AND_COUNT` times differ 12.2x for
  only 1.35x more input reads (caveat 2 above) — real prerequisite for
  trusting either Cell Ranger figure as a baseline; check `_perf`'s
  per-substage breakdown within `ALIGN_AND_COUNT` for both pipestances
  before assuming it's a library-biology difference rather than a
  cluster-scheduling artifact.
- Whether Cell Ranger's `ALIGN_AND_COUNT` is running at comparable real
  parallelism to sbxpress's `--threads`/`--shard-concurrency` or whether
  one side is spending a lot of wall-clock under-parallelized (e.g.
  serialized phases, lock contention) — real `/usr/bin/time -v` or a
  sampling profiler (see `PROFILING.md`) on a `molecules` run would settle
  this rather than guessing from aggregate numbers alone.

## `sbxpress pipeline` orchestrator: real 5-shard validation, 2026-09-20 — handoff, stopped mid-investigation

Real end-of-session state, so a fresh session can pick this up cold. All
work below sits **uncommitted on the `pipeline-orchestrator` branch**
(on top of commit `6c4357f`) — do not assume it's saved anywhere else.

**Goal**: validate `sbxpress pipeline` (built this session) against the
existing shell-script pipeline shape, on 5 real raw FASTQ shards
(`raw_fastqs/*-00000{1..5}.fastq.gz`), before trusting it for a real CP1-CP4
rerun. Test harness and both paths' full output live at
`/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06/orchestrator_validation_5shard/`
(`shell/` = shell-script-style path, `pipeline/` = orchestrator path;
`compare.log`/`A_*.log`/`B_pipeline.log` have the real timestamped logs).

**One real, systemic bug found and fixed this session (applied, tests
pass, not committed)**: every binary in this tree except `sbx10x-annotate`
let `tracing`'s log output default to **stdout** instead of stderr. This
silently corrupted `sbx10x --print-manifest`'s JSON (a log line landed
ahead of the JSON on the same stream), breaking the orchestrator's manifest
fetch. Fixed in `src/main.rs` and `src/bin/{cr-bam-to-mtx,sbxpress-guides,
sbxpress-probe-demux,sbxpress-probeset-align,sbxpress-protospacer-call,
sbx10x-molecules}.rs` — all now route to stderr via `.with_writer(std::io::stderr)`,
matching `sbx10x-annotate.rs`'s existing (correct) pattern. Also added the
missing `--molecule-shard-split` CLI flag to the orchestrator (it predated
`--shard-split`'s existence and was never wired up) — real, needed for the
`--shard-split` win to actually apply when using `sbxpress pipeline`.

**Path A (shell-script style) succeeded end to end**: extract 774.7s →
align ~347s → molecules 964.1s = **2177s total**, real numbers, real output.
This is the trustworthy baseline for the next attempt.

**Path B (orchestrator), after the tracing fix, got further but surfaced
two more real, unresolved problems — investigation stopped here, pick up
fresh:**

1. **CRISPR/feature molecules invocation crashes on real multi-chunk
   input.** `src/orchestrator.rs`'s `other_bams`/"molecules (feature/CRISPR)"
   block (the one building `args` for `sbx10x-molecules ... --tagseq`, no
   `--shards`/`--shard-concurrency` passed at all — noted as a real gap
   earlier this session but not fixed) defaults to `--shards 1`. That's fine
   for a single CRI file but this run had **24** CRI-routed chunks, and
   `sbx10x-molecules`'s non-sharded path only accepts one input — real
   crash: `multiple input BAMs are only supported with --shards > 1 (got 24
   inputs, --shards 1)`. **Fix is straightforward and already has a working
   template two lines above it**: mirror the GEX block's `--shards`/
   `--shard-concurrency`/`--shard-split` wiring (same `shards`/`concurrency`/
   `cli.molecule_shard_split` variables already in scope). Not yet done.

2. **Much bigger, unexplained, and the real priority for next session:**
   GEX molecules' own metrics show a massive, unexplained gap between the
   two paths for what should be equivalent input:
   | | Path A (shell) | Path B (orchestrator) |
   |---|---|---|
   | records | 85,601,730 | 88,563,919 |
   | eligible_records | **70,073,891 (82%)** | **2,481,511 (2.8%)** |
   | molecules_counted | 69,660,540 | 2,472,428 |

   Record counts are close (B has ~3M more, plausibly explained by B's
   manifest routing "Undetermined" chunks to GEX by default — 12 extra
   chunks A's script never fed to `probeset-align` at all, a real, separate,
   smaller discrepancy worth resolving on its own). **But the eligible-rate
   collapse from 82% to 2.8% is not explained by that** — Undetermined
   reads being low-quality wouldn't drag the *entire* population down this
   far given they're a minority of the input.
   
   One concrete, un-investigated clue at the point this session stopped:
   the *same-named* aligned chunk file is a wildly different size between
   the two paths —
   `pipeline/aligned/TSeq100_L1_GEX_S1_L001_R2_001.bam` is **686MB**,
   `shell/aligned/TSeq100_L1_GEX_S1_L001_R2_001.bam` is **225MB** — a ~3x
   size difference for what should be the same chunk of the same raw input
   run through the same `sbxpress-probeset-align` binary with the same
   `--probe-set`. This was the very next thing to investigate (`samtools
   view` both, compare read counts/CB/GX tag presence directly) when the
   session was stopped — **do this first**. Real candidate causes to check,
   none confirmed yet: (a) the orchestrator's align invocation passing
   different/fewer/extra flags to `sbxpress-probeset-align`/`sbxpress-guides`
   than the manual `xargs` loop in `compare.log`'s Path A used (compare
   `run_extract_align_materialize`'s exact spawned args against
   `A_align.log`'s), (b) the extra "Undetermined" chunks somehow being
   miscounted into the wrong sample bucket, (c) something in the
   `--no-gzip`/materialize path corrupting or duplicating records between
   `extract` and `align` that a size difference this large would be
   consistent with. **Not yet root-caused — do not trust `sbxpress
   pipeline`'s molecule counts until this is understood and fixed.**

**Real, trustworthy result from earlier today, unaffected by the above**:
`--shard-split 8` on top of `--shards 16 --shard-concurrency 12`, run
directly via `sbx10x-molecules` (not through the orchestrator), clean
uncontended box: **9,430.0s → 8,340.5s, ~11.6% faster**, byte-identical
aggregate metrics. This result stands regardless of the orchestrator issues
above.

**Next session should, in order**: (1) `samtools view` both
`TSeq100_L1_GEX_S1_L001_R2_001.bam` copies and diff real read/tag content to
root-cause the eligible-rate collapse — this blocks trusting any
orchestrator-produced molecule counts; (2) fix the CRISPR/feature
`--shards` wiring gap (small, mechanical, template already exists); (3) once
both are fixed, rerun this same 5-shard comparison clean before considering
a real CP1 run; (4) only then revisit the original ask (rerun CP1-CP4 with
all changes).
