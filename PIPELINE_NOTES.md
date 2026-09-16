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
| `samtools sort` | `-@` mainly parallelizes *compression*, not the core in-memory sort-comparison step — real observed pattern: bursty CPU (idle while filling a buffer, spikes to several cores while compressing+flushing one chunk), not sustained at the requested thread count | **`-m` (memory per thread) is the real lever**, not `-@`, once you have RAM to spare | Directly memory-gated: bigger `-m` → fewer, larger chunks → far less repeated read/sort/compress/flush overhead | Real before/after on this box: `-m 2G` → ~70MB/s effective read throughput (many small ~37.5GB chunks, one every ~10min); `-m 8G` → ~254MB/s (3.6x) on the identical file. **Always check `free -h` before raising `-m`** — this is a real, direct memory commitment, not just a suggestion. |
| `samtools index` | Low | Low | Read-heavy, fast | **Requires coordinate-sorted input** — fails outright ("NO_COOR reads not in a single block") on `samtools cat` output. Sort first, every time. |
| `sbx10x-molecules` | Default `--evidence-workers 1` is **fully serial** and dramatically slower at billion-read scale despite a real parallel path (`--evidence-workers N`, commits `3b7366c`/`9afc664`) existing in the binary since 2026-09-02 — that path additionally **requires the input BAM to be indexed**, so it only activates after the sort+index steps above. | **By far the dominant memory consumer in the whole pipeline.** Real observed: **~976GB RSS** for one ~9B-read / ~4.02B-raw-UMI-node library (`TSeq100_L1_GEX`), serial mode. Root cause (verified by reading the code, not guessed): `evidence` (`HashMap<RawKey,Evidence>`), `roots` (`HashMap<RawKey,RawKey>`), and `decisions` (`HashMap<RawKey,Decision>`) are all genuinely alive *simultaneously* at the point `decisions` gets built — each ~4B entries, ~24–48 bytes/entry raw plus real hashbrown overhead. No cheap safe reduction found; would need real restructuring (e.g. streaming `decisions` construction while draining `roots`/`evidence`) to meaningfully cut peak memory. | Read-heavy (evidence collection), write-heavy (final BAM rewrite) | **Two real libraries' worth of this can't run concurrently** on a 1.5TB box — confirmed L1 alone used ~976GB, leaving only ~552GB free, not enough headroom for a second comparably-sized library's molecules pass at the same time. This is the single biggest scheduling constraint discovered this session. Whether `--evidence-workers N>1` has a meaningfully *different* peak-memory profile than serial mode (same total live data, just more concurrent access to it) is an open question worth confirming once real data is available — plan to note the answer here once observed. |
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
3. **Bound `molecules`' peak memory by design, not by hoping it fits.**
   Given `evidence`/`roots`/`decisions` need to coexist today, a real
   redesign — e.g., partitioning the genome/synthetic-contig space into
   bounded-memory batches, writing `decisions` incrementally per batch
   instead of materializing one global map — would let memory scale with
   `--evidence-chunk-bases`/batch size instead of total library size. This
   is the single highest-value architectural change available, given
   `molecules` is both the slowest and the most memory-hungry stage by a
   wide margin.
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

- Real peak memory for `sbx10x-molecules --evidence-workers 64` (parallel
  mode) vs. the ~976GB serial-mode baseline, once the corrected pipeline
  actually runs a library through the parallel path end to end.
- Real wall-clock improvement from `-m 8G` sort + `--evidence-workers 64`
  molecules, full library, once L1 completes under the new settings —
  update the table above with the actual number instead of the
  in-progress read-throughput sample.
