# Gene-collision resolution for high-depth pooled samples — proposed statistical fix

Status: **proposed, not implemented.** Pinned for later; no code changes made against this doc.

## Background

`sbx10x-molecules`'s gene-collision resolution (`build_decisions`, `gene_collision_index` /
`gene_collision_resolve` in `src/bin/sbx10x-molecules.rs`) groups molecule roots strictly by
`(cell, UMI)`, with gene identity dropped. When more than one gene claims the same `(cell, UMI)`,
it keeps the single highest-read-support gene (`Won`) and discards the rest (`Lost`), or discards
all of them on an exact tie (`Tie`).

This is a reasonable cleanup step at the depth Cell Ranger-style single-cell pipelines are designed
for (~20-50K reads/cell): aggregate UMI-space occupancy per cell is low relative to the 12bp UMI
alphabet (4^12 = 16,777,216 values), so a genuine *random* collision between two independent
molecules of different genes is vanishingly rare. When one is observed, it's almost always a real
artifact (index hopping, an uncorrected sequencing error, ambient contamination), and discarding
the loser costs a negligible fraction of the data.

DZEC03-style runs are a different regime: each "cell barcode" is a whole pooled sample at
2-25M reads, not a single cell at 20-50K reads. Measured on DZEC03's real samples: average raw
`(cell, gene, UMI)` demand is ~35.8M per sample against the 16.8M-value UMI space — ~2.1x
oversubscribed in aggregate. At that occupancy, cross-gene UMI collision is the *expected*
structural outcome of the experiment design, not rare noise. Empirically: 47,087,038 of
57,874,260 observed `(cell, UMI)` pairs (81%) were claimed by more than one gene, and the current
policy discards 368,816,941 "lost" molecules plus every member of 11,689,736 tied groups
(run-wide) as a result. The top colliding genes (MT-CO3/mitochondrial, ALB/chr4, EEF1A1/chr6,
HBG1-2/chr11, ...) are genomically unrelated, so these are not read-mapping ambiguity — they're
independent real molecules losing to each other by chance.

Three directions were discussed:

1. **Drop cross-gene adjudication entirely** — count distinct UMIs per `(cell, gene)`
   independently, closer to standard bulk UMI-tagged RNA-seq quantification. Simple, but accepts
   unresolved double-counting inflation wherever two genes genuinely do share ambiguous evidence.
2. **Restrict adjudication to genuinely plausible ambiguity** (genomic overlap/proximity, or
   reads already flagged multi-mapping at the annotation stage) — chosen as the preferred
   direction, see caveats below (discussed in chat, not this doc).
3. **Statistical / fractional resolution** — this document.

## The proposal (option 3)

Instead of a deterministic winner-take-all discard, treat a `(cell, UMI)` collision as a
statistical estimation problem: apportion the molecule's expected true origin across the
competing genes rather than assigning it entirely to one and zeroing out the rest.

### 3a. Fractional / proportional assignment (simplest version)

For a `(cell, UMI)` collision among genes `g_1..g_k` with read-support counts `r_1..r_k`, assign
gene `g_i` a fractional molecule count of `r_i / sum(r)` instead of `1` to the winner and `0` to
everyone else. Sums to exactly 1 molecule per collision group, same as today, but the credit is
distributed proportionally to evidence rather than all-or-nothing. Ties resolve naturally (equal
shares) instead of being discarded outright.

- Downstream consumers (the `matrix`/`cells` stage, MEX/H5 output) currently assume integer
  counts. Fractional counts either require those stages to accept floats, or a rounding/stochastic
  rounding step at materialization (round each fractional count with probability equal to its
  fractional part, preserving the expected value across genes/cells — standard practice for
  fractional-to-integer count reconciliation).
- Loses the "how much was ambiguous" signal that raw `Won`/`Lost`/`Tie` currently gives in
  `metrics_summary` unless a parallel diagnostic is kept.

### 3b. Occupancy-based statistical correction (closer to the real mechanism)

Rather than resolving each collision individually, correct at the gene level using the same
birthday-paradox/Poisson reasoning used to explain the 81% empirical collision rate:

- For a given cell, let `D` = total raw `(cell, gene, UMI)` demand summed across all its genes,
  and `M` = 4^12 = 16,777,216. Per-slot expected occupancy is `λ = D / M`.
- For a single gene with raw UMI count `n_g` in that cell, the expected number of its UMI draws
  that *also* landed on a UMI independently claimed by at least one other gene is estimable from
  `λ` and `n_g` (a two-part Poisson/birthday model: probability a given UMI value is "shared" by
  something else in the cell is `1 - e^{-λ'}` where `λ'` excludes gene `g`'s own contribution).
- Instead of deciding per-collision winners, inflate/deflate each gene's *final* distinct-UMI
  count by the model-predicted collision rate, similar to how per-gene UMI-saturation correction
  (e.g. the classic `-ln(1 - observed/possible)` estimator) already corrects for *within-gene*
  UMI collision. This would need to be derived and validated as a *cross-gene* extension of that
  same idea — not a small undertaking, but statistically the most principled of the three options
  since it corrects the estimate rather than making a case-by-case guess.
- Requires deciding what to do with the underlying reads/records for output-BAM tagging
  (`xf`/`XW` representative tags) since a statistical correction operates on counts, not on
  which specific read gets flagged `representative` in the output BAM — likely still needs a
  per-collision representative pick for BAM-tagging purposes even if the *count* used downstream
  is corrected.

### Open questions before implementing either variant

- Validate the Poisson/birthday-paradox occupancy model against this dataset's *actual* observed
  collision rate (81%) at the gene level, not just in aggregate, before trusting a correction
  derived from it.
- Decide whether correction happens in `sbx10x-molecules` (changing `Decision`/metrics semantics
  for every consumer of this tool) or as a post-hoc adjustment layer downstream of the current
  integer-count output (safer for other `software_shared` users depending on today's semantics).
- Whatever is chosen needs the same validation rigor as the `correct_group` performance fix: unit
  tests with a hand-constructed collision scenario with a known expected fractional/corrected
  output, plus a real-data before/after comparison — harder here than for `correct_group` because
  there is no single objectively "correct" answer to check against, only a more or less
  statistically defensible one.

## Non-goals of this document

Not proposing anything for option 1 or option 2 here — this file is scoped to the "put a pin in
option 3" request. Option 2's implementation caveats (genomic-proximity as an imperfect proxy for
true ambiguity, NUMT-style distant homology, engineering lift for a GTF-derived overlap/adjacency
lookup) were covered in chat, not written up separately yet.
