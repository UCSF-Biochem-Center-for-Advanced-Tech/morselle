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
