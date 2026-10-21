#!/usr/bin/env bash
# Everything a change must pass before it is committed or merged.
# CI (.github/workflows/ci.yml) runs exactly this script.
#
#   scripts/check.sh          # full check
#   scripts/check.sh --quick  # skip clippy (fast local loop)
#
# What each gate guards against (each one was a real 0.5.2-0.5.8 failure):
#   build --all-targets  benches/examples are not built by `cargo test`
#                        (0.5.1: a bench silently stopped compiling)
#   cargo test           unit, stage-contract (pipeline_chain_e2e), golden
#                        output (golden_outputs), docs flags (docs_commands),
#                        pipeline scripts (pipeline_scripts)
#   test inventory       tests deleted during a "refactor" (0.5.7)
#   500-line ceiling     project rule (keep code files modular)
#   rustfmt, clippy      style drift; clippy's deny-by-default lints
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
QUICK=0
[[ ${1:-} == --quick ]] && QUICK=1
step() { printf '\n==> %s\n' "$*"; }

step "rustfmt"
cargo fmt --check

step "build every target (lib, bins, tests, benches, examples)"
cargo build --release --all-targets

step "tests"
cargo test --release

step "test inventory (no silently deleted tests)"
scripts/check_test_inventory.sh

step "500-line ceiling for Rust sources"
over=$(find src tests benches examples -name '*.rs' -print0 | xargs -0 wc -l | awk '$2 != "total" && $1 > 500')
if [[ -n $over ]]; then
    echo "files over 500 lines (split them):"; echo "$over"; exit 1
fi

if (( QUICK == 0 )); then
    step "clippy (deny-by-default lints fail the build)"
    cargo clippy --release --all-targets
fi

printf '\nAll checks passed.\n'
