# Agent Instructions & Code Quality Rules for `sbxpress`

## Code File Structure & Modularization Guidelines
- **Target File Size**: Code files should ideally be kept between **150 and 300 lines**.
- **Hard Ceiling**: No single code file may exceed **500 lines**.
- **Modularization**:
  - Break monolithic files into clean, logical submodules (e.g., `module/mod.rs`, `module/types.rs`, `module/algo.rs`, `module/tests.rs`).
  - Keep module entry points concise (`mod.rs` re-exporting submodules).
  - Ensure zero breaking changes to public APIs, function signatures, or CLI flags.
