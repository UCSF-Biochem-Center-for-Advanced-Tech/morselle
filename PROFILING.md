# Profiling

Three complementary tools, from coarsest to finest. Use them in that order:
find which *stage* is slow first, then which *function within that stage*,
then (once available) which *instruction*.

See `PIPELINE_NOTES.md` for real, measured CPU/memory/I/O resource profiles
of every pipeline stage at true production scale (32.6B-read full-dataset
runs) — complementary to this doc's per-function profiling: that one is
about *whole-pipeline resource planning and scheduling*, this one is about
*finding the hot function once you know which stage to look at*.

## 1. Cross-stage comparison (works today, zero setup)

Every stage (`extract`, `annotate`, `molecules`, `matrix`, `cells`,
`probe-demux`) already writes `elapsed_secs`/`process_secs` and a
`reads_per_sec`/`records_per_sec` figure into its own `metrics_summary.json`
(or `probe_demux_metrics.json`), and logs the same numbers as a `RUN SUMMARY`
line at INFO level. Run the pipeline once end-to-end and just read those back
— no instrumentation to add, this is free.

For memory/syscall-level resource stats per stage (max RSS, major/minor page
faults, context switches, CPU%), wrap any stage in `/usr/bin/time -v`
(**not** the bash builtin `time` — use the full path):

```sh
/usr/bin/time -v sbxpress extract --fastqs ... -o out --minimap2 2>&1 | tail -20
```

## 2. Hot-path microbenchmarks: `cargo bench` (works today, zero setup)

`benches/pipeline_hot_paths.rs` (Criterion, `harness = false`) benchmarks the
hot functions shared across stages in isolation, on synthetic-but-realistic
data (same shape as `orient.rs`'s own unit-test fixtures) — portable, no
external files, no OS profiler dependency, gives proper statistics
(mean/variance across many iterations, not one wall-clock number).

```sh
cargo bench                                    # full run, saves a baseline
cargo bench -- orient_extract                  # filter by name substring
cargo bench -- --quick                         # fast/noisy, for iteration
cargo bench -- --save-baseline before          # ...then after a change:
cargo bench -- --baseline before               # regression-checks against it
```

HTML reports (with distributions, not just point estimates) land in
`target/criterion/<bench-name>/report/index.html` when `gnuplot` is
installed; falls back to the `plotters` backend otherwise (still produces
plots, just simpler ones).

**Add a benchmark per new hot path** as new phases land (Phase 4's anchor
matcher, Phase 5's GMM fit, etc.) — each `c.bench_function` in that file is
independently filterable, so the harness scales with the plan without
needing a new file per phase.

## 3. Real CPU-sampling flamegraphs: `samply` (needs one `sudo` command)

`samply` is already installed (`~/.cargo/bin/samply`, via `cargo install
samply`) and gives real instruction/function-level flamegraphs — the thing
to reach for once (1) or (2) has told you *which* stage/function to look
inside. It's currently blocked on this machine:

```
$ samply record -- ./target/release/sbx10x ...
'/proc/sys/kernel/perf_event_paranoid' is currently set to 4.
In order for samply to work with a non-root user, this level needs
to be set to 1 or lower.
```

**One-time fix** (needs sudo — not something this session can do
non-interactively):

```sh
echo '1' | sudo tee /proc/sys/kernel/perf_event_paranoid
```

That's a runtime-only change (resets on reboot); to make it permanent, add
`kernel.perf_event_paranoid = 1` to `/etc/sysctl.d/` and `sudo sysctl -p`.

Once unblocked, build with debug symbols (a dedicated `profiling` Cargo
profile — same optimizations as `--release`, just keeps DWARF info, so the
normal installed `--release` binaries stay lean and unaffected):

```sh
cargo build --profile profiling
samply record -- ./target/profiling/sbx10x --fastqs ... -o out --minimap2
# opens a local web UI (firefox-profiler) with the flamegraph;
# --save-only -o profile.json.gz to save without opening a browser (e.g. over SSH)
```

Do this once per stage binary (`sbx10x`, `sbx10x-annotate`,
`sbx10x-molecules`, `cr-bam-to-mtx`, `sbxpress-cells`,
`sbxpress-probe-demux`) against a representative-sized real input (a handful
of the CAT06 shards, not the full 754-file run) to find that stage's actual
hot spots.

**Headless viewing note**: `samply load` only serves an interactive browser
UI (Firefox Profiler) and needs `--no-open`/a forwarded port on a headless
box. Symbol names are *not* baked into the saved `.json.gz` — the UI
resolves them lazily. To get a text report without a browser, resolve
addresses yourself with `addr2line` against the exact binary path recorded
in the profile's `libs` table (each `frameTable.address` is already a
file-relative offset, no ASLR-base subtraction needed) and aggregate by leaf
frame for self-time. Recipe (used to produce the finding below, not checked
in — ~50 lines, easy to rewrite): load the profile JSON, walk
`threads[].samples.stack -> stackTable.frame -> frameTable.{address,func} ->
funcTable.resource -> resourceTable.lib -> libs[].path`, batch-call
`addr2line -e <path> -f -C <hex addrs...>`, and `Counter()` the results.

**Watch out for thread-count fragmentation.** This box has 512 logical
CPUs, so `rayon`'s default thread pool spins up 512 OS threads — on a small
input, samples end up spread so thin across threads that there's nothing to
aggregate. Cap threads for a profiling run (`--threads 32`–`64`) and use a
few-million-read slice, not a tiny sample and not a full 2GB shard (the
extra wall-clock isn't worth it for a stable profile).

### Real finding (2026-09-14): `extract` is gzip-bound, not algorithm-bound

First real flamegraph run — `sbx10x extract --minimap2
--probe-barcode-whitelist flex-v2-384.txt` on a 2M-read CAT06 slice,
`--threads 64`, 37,613 samples aggregated across all threads:

| category | share |
|---|---|
| zlib (gzip **output** compression — `deflate_slow`/`longest_match`/`slide_hash`) | **53.9%** |
| core algorithm (`orient`/`extract`/`anchors`/`homopolymer`) | 18.0% |
| libc (syscalls, `write`, malloc) | 9.1% |
| pipeline glue (`lib.rs::process_one`, emit) | 4.5% |
| Rust std/alloc internals (inlined into hot loops) | ~8.9% |
| demux + probe-barcode matching (this session's whole Phase 0/1 addition) | **0.12%** |

**Takeaway**: over half of extract's CPU time is spent compressing the
output FASTQ, not in any biology/matching code — confirming in hard numbers
what `config.rs`'s own `--gzip-level` doc comment already suspected
("output compression is the run's wall-clock floor at scale"). This
session's whole Flex addition (probe-barcode locate/tag) is immaterial to
throughput (0.12%) — there was never a real risk of it being a bottleneck.

### Follow-up (2026-09-14): added `--no-gzip`, measured an 8x win

Acted on the finding above rather than just tuning `--gzip-level`/`--split`:
`extract`'s typical next pipeline step is a mapper reading the FASTQ once,
not long-term storage, so compression there is close to pure waste. Added
`--no-gzip` (writes plain `.fastq`; `io.rs`'s `ShardSink` now wraps either a
`GzEncoder` or a plain `BufWriter`, `level: Option<u32>` — `None` means
uncompressed, not "level 0 gzip"). Measured on the same 2M-read slice,
`--threads 64`:

| | elapsed | reads/sec |
|---|---|---|
| gzip (default, level 4) | 19.1s | 104,888 |
| `--no-gzip` | 2.4s | 837,694 |

**~8x, not the ~2.2x the 54% CPU share alone would predict.** The extra
speedup is lock contention, not just CPU: every worker thread shares one
mutex-protected encoder per output target (`PairShardWriter`), so with gzip
on, only one thread can be deflating at a time while the rest queue up
behind that lock — removing the encoder removes a true serialization point,
not just work.

Verified this is a real, usable capability, not just a raw-throughput
number: piped `extract --no-gzip`'s output through a named pipe
(`mkfifo`) straight into `minimap2` with no intermediate file on either
end — `minimap2` read all 1,649,018 records with zero loss, and the
`CR`/`CB`/`UR`/`UB` tags survived through minimap2's `-y` untouched:

```sh
mkfifo out/sbxpress_S1_L001_R2_001.fastq
sbx10x --fastqs in.fastq.gz -o out --minimap2 --no-gzip &
minimap2 -a -y ref.mmi out/sbxpress_S1_L001_R2_001.fastq > aligned.sam
```

`--gzip-level`/`--split` are still there for when the output genuinely needs
to be stored (e.g. `--cellranger` compatibility mode, or archiving), but for
a mapper-fed run, `--no-gzip` (+ a named pipe, if avoiding the intermediate
file matters) should be the default choice going forward — this is now
standard practice for extract runs in this session, and worth adopting
generally.

### Full chain (2026-09-14): the win survives past mapping too

Ran the real production shape — `extract` → `minimap2` (splice:sr,
real GRCh38 + GENCODE v44 index/junctions reused from the CAT02 project,
same reference the README documents) → `annotate` → `samtools sort` →
`molecules` → `matrix` → `cells` — on the same 2M-read CAT06 slice, gzip
vs. `--no-gzip`, `--threads 64` throughout:

| stage | gzip | `--no-gzip` | speedup |
|---|---|---|---|
| extract | 17.4s | 2.4s | **7.3x** |
| `minimap2 \| annotate \| sort` (one pipe) | 45.3s | 25.6s | **1.8x** |
| **extract + map/annotate/sort total** | **62.7s** | **28.0s** | **2.2x** |

Both variants produced byte-identical downstream results (1,649,018 records
in both `annotated.bam`s, same `annotate` metrics) — `--no-gzip` changes
speed, not output. The mapping-stage win is real, not noise: minimap2's own
reported "Real time" dropped 44.2s → 24.5s in lockstep with the pipe's wall
time, meaning a large chunk of minimap2's own runtime was previously spent
inflating our gzip output on the way in — `--no-gzip` removes a
compress/decompress round trip on *both* ends of that handoff, not just the
encode side.

Molecules/matrix/cells ran clean end-to-end too (`sbx10x-molecules` →
`cr-bam-to-mtx` → `sbxpress-cells`, no errors), but surfaced something
independently worth knowing, unrelated to compression: only 59 of the
1,649,018 records were real countable molecules (55 genes, 59 barcodes).
That's expected, not a bug — CAT06 is 10x Flex/CRISPR data, whose R2 is a
probe-barcode + probe/guide payload (see docs/archive/FLEX_SUPPORT_PLAN.md), not free
poly-dT cDNA, so splice-aware whole-genome alignment of that payload
finding almost nothing real is exactly what the architecture predicts.
This run confirms the mechanics of the full binary chain on real data, but
is **not** the "real annotated+molecule GEX BAM" Phase 2/3 are still
waiting on — that needs Flex's actual probe-set-reference alignment
strategy, not standard genome mapping, and is Phase 6 territory.

### `molecules` / `matrix` per-record cost (2026-10-07)

Investigated the open question in PIPELINE_NOTES.md ("`molecules` alone is
2.7x longer than Cell Ranger's whole ALIGN_AND_COUNT"). Callgrind on a
synthetic 200k-record annotated BAM with the real tag layout (minimap2 +
extract comment + annotate tags, 24 per record) showed the time was
per-record overhead, not the UMI algorithms:

| cost (baseline `sbx10x-molecules`) | share of instructions |
|---|---|
| zlib deflate + crc32 (BGZF output) | ~35% |
| `bam_aux_next`/`bam_aux_get` — every `record.aux()` / `push_aux` / `remove_aux` rescans all 24 tags | ~21% |
| SipHash (std `HashMap`) | ~6% |
| `from_utf8` / `AuxIter` decoding tags never used | ~10% after the first fix |

Fixes (all output-identical): one raw aux walk per record
(`bam_helpers::for_each_aux`, `MoleculeTags`), FxHash, scan-free tag
appends (`molecules/aux_write.rs`), rust-htslib `libdeflate`, and sorted
runs instead of `HashMap<_, Vec<_>>` group indexes. Result on a 2M-read
BAM: `molecules` 16.2 s → 7.6 s, 18.0 G → 7.6 G instructions, peak RSS
339 → 256 MB; `matrix` 4.2 s → 2.7 s.

What is left in `molecules`: ~40% is BGZF compression of the output BAM
(htslib default level). A `--compression-level` option (or level 1–4 for
this intermediate) is the next lever, at some cost in file size. Peak
memory is now at decision materialization (evidence + roots + molecules +
decisions); having the write pass look up `roots` + `molecules` directly
instead of building the `decisions` table would cut it roughly another
quarter, but changes `build_decisions`' public return type.

Reproduce: generate a BAM with the real tag set, then
`valgrind --tool=callgrind sbx10x-molecules in.bam -o out.bam --threads 1`
and `callgrind_annotate`. Instruction counts are deterministic, unlike
wall time on a shared box.
