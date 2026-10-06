# Pipeline resource & architecture notes

> **Status**: 🟡 **`probeset-align` multi-threading completed in v0.5.3**, **`scripts/sbxpress-pipeline.sh` completed in v0.5.4**, and **master roadmap consolidated in v0.5.6**. `(cell, sample)` sharded RAM bounding for `sbx10x-molecules` is tracked in Priority 1 of [FUTURE_ROADMAP.md](FUTURE_ROADMAP.md).

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
  design and is what `read_depth_comparison.sh`'s CP1–CP4 checkpoints
  actually run. See "`sbx10x-molecules` real wall-clock vs. Cell Ranger,
  2026-09-19" below — it solved the memory problem but is now the
  pipeline's dominant real *time* cost, not yet investigated.

## `sbx10x-molecules` real wall-clock vs. Cell Ranger, 2026-09-19

Separate, not-yet-started investigation — self-contained enough to pick up
in a fresh session with no memory of how it was found. Written from CAT06's
real CP1 run (`read_depth_comparison_plan.md`, restarted 2026-09-18 with
`sbxpress` 0.4.0's composite-cell-barcode fix); see that doc for full CAT06
context if needed, but this specific question stands on its own.

**Real numbers, both sides measured directly, not estimated:**

- `sbx10x-molecules` on CP1's `L1_gex.bam` (`--shards 8
  --shard-concurrency 8 --threads 16`): **33,427.9s (9.28h)** —
  `read_depth_comparison_out/cp1/molecules/L1_gex.metrics.json`'s
  `elapsed_secs`.
- Same stage on CP1's `L2_gex.bam`, same flags, same real box, same day:
  **10,178.0s (2.83h)** — `.../cp1/molecules/L2_gex.metrics.json`. **3.3x
  faster than L1, unexplained.** Both libraries came from the same
  159-shard raw-read subset (`read_depth_comparison_out/cp1_shards.txt`),
  so this isn't an input-size difference by construction — worth checking
  real per-library read/molecule/collision counts (`probeset-align`'s own
  metrics, `molecules`' `collision_groups`/`raw_nodes` figures) for what
  actually differs.
- Cell Ranger's real, directly comparable stage — `ALIGN_AND_COUNT`
  (alignment + UMI counting combined, one stage) on the same L1 library,
  real Illumina full-depth run: **12,569.6s (3.49h), 207.7 core-hours** —
  from the real preserved Martian pipestance's `_perf` JSON:
  `illumina_data/RZ16443_L1_multi_full/_perf` (flatten `forks[].fork_stats`
  per entry; fqname
  `...SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.COUNT_GEM_WELL_PROCESSOR._BASIC_SC_RNA_COUNTER._MATRIX_COMPUTER.ALIGN_AND_COUNT`).
  **sbxpress's `molecules` stage alone, on L1, is 2.7x longer than Cell
  Ranger's entire align+count stage** — real, large, unexplained gap; not
  yet investigated at the code/algorithm level (only the earlier memory
  investigation above has been done, which explains RSS, not wall-clock).

**Not yet checked, real starting points for a fresh investigation:**

- Where real wall-clock actually goes across `molecules`' own phases
  (`evidence` collection / `umi_correct` / `molecule_aggregate` / final BAM
  write) — the binary already logs phase-boundary progress lines (see
  `[sbxpress molecules] phase=...` in `cp1/logs/molecules_L{1,2}.log`);
  real per-phase timing hasn't been extracted and diffed between L1 and L2.
  yet.
- Whether `--shards 8` is actually a good choice for wall-clock (it was
  chosen for a real memory constraint — see the OOM finding above — not
  profiled for speed); a real per-shard timing breakdown from the same logs
  would show whether shards are balanced or whether a few stragglers
  dominate (the same kind of imbalance `--split`'s round-robin bug caused
  elsewhere in this pipeline, real precedent worth checking for here too).
- Whether Cell Ranger's `ALIGN_AND_COUNT`'s 207.7 core-hours at 12,569.6s
  wall-clock (~16.5 average concurrent cores) vs. sbxpress's `--threads 16`
  (also nominally 16-wide) is actually running at comparable real
  parallelism, or whether one side is spending a lot of that wall-clock
  under-parallelized (e.g. serialized phases, lock contention) — real
  `/usr/bin/time -v` or a sampling profiler (see `PROFILING.md`) on a
  `molecules` run would settle this rather than guessing from the
  aggregate numbers alone.
