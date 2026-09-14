# Profiling

Three complementary tools, from coarsest to finest. Use them in that order:
find which *stage* is slow first, then which *function within that stage*,
then (once available) which *instruction*.

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
The two existing levers worth tuning first, before touching any algorithm
code: `--gzip-level` (lower = faster, larger files) and `--split N`
(parallelizes one sample's compression across N writer threads). Neither
has been benchmarked against this profile yet — natural next step.
