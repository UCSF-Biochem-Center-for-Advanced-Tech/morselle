# Project Guidelines & Agent Instructions for `sbxpress`

## Code File Structure & Modularization Guidelines
- **Target File Size**: Code files should ideally be kept between **150 and 300 lines**.
- **Hard Ceiling**: No single code file may exceed **500 lines**.
- **Modularization**:
  - Break monolithic files into clean, logical submodules (e.g., `module/mod.rs`, `module/types.rs`, `module/algo.rs`, `module/tests.rs`).
  - Keep module entry points concise (`mod.rs` re-exporting submodules).
  - Ensure zero breaking changes to public APIs, function signatures, or CLI flags.

## Testing & Change Discipline (read before committing)

These rules exist because 0.5.7/0.5.8 shipped with `sbxpress molecules`
crashing on every real BAM while all tests passed (see CHANGELOG 0.5.9).

- **Run `scripts/check.sh` before every commit** (`--quick` skips clippy).
  `git config core.hooksPath scripts/hooks` makes `git push` run it.
  It is exactly what CI runs. Do not commit, bump versions, or update
  CHANGELOG/docs as "done" while it fails.
- **A refactor must not change behaviour.** `tests/golden_outputs.rs` pins
  the exact output of `molecules` + `matrix` in every mode. If it fails
  after a refactor or optimisation, the change is wrong, not the golden
  file. Only regenerate (`SBXPRESS_UPDATE_GOLDEN=1 cargo test --release
  --test golden_outputs`) when the behaviour change is the point of the
  commit, in a separate commit whose message explains why.
- **Never delete, weaken, or rewrite a test to make it pass.** Removing or
  renaming a test fails `scripts/check_test_inventory.sh` until
  `tests/test_inventory.txt` is updated; do that only with a stated reason.
  Moving tests between files is fine (run `--update` and keep every name).
- **Keep comments when moving code.** Comments here record real-data
  findings and design reasons; they are not noise to strip.
- **Test stages through their real neighbours.** A change to what one
  stage writes or reads (BAM tags, matrix format) needs a case in
  `tests/pipeline_chain_e2e.rs`, using `tests/common` fixtures that look
  like real data (existing `UB`/`ZG` tags, multi-gene `GX`, etc.) -- not a
  hand-built record shaped like what the consumer expects.
- **Docs show only real commands.** `tests/docs_commands.rs` checks every
  `sbxpress ...` and pipeline-script flag in README/MANUAL/INSTALL against
  `--help`. Shell scripts are parsed, shellchecked and smoke-run by
  `tests/pipeline_scripts.rs`; add a smoke case for new script modes.
- **One version bump per release**, after `scripts/check.sh` passes, with
  a git tag (`vX.Y.Z`).
