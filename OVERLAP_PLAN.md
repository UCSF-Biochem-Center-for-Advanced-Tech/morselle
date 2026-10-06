# Illumina paired-end overlap + YC tag — implementation plan

> **Status**: 🟢 **Engine Completed in v0.3.7**. `sbxpress illumina --overlap` and `YC` duplex-consensus tagging are merged on `main`. Exporting insert-size distribution histograms is tracked in Priority 3 of [FUTURE_ROADMAP.md](FUTURE_ROADMAP.md).

## Purpose

Add an opt-in paired-end overlap/merge stage to the `sbxpress illumina`
converter, and encode the merge structure losslessly in a Roche-compatible `YC`
tag on the emitted read.

Two motivations, in priority order:

1. **Fairer arm-to-arm comparison.** The SBX arm sees whole fragments. The
   Illumina arm currently sees only R2 — [`convert_pair`](src/bin/sbxpress-illumina.rs:107)
   takes `r1_seq[..28]` for CB+UMI and **discards the rest of R1**. On a
   long-R1 run that is real, paid-for cDNA sequence being thrown away. Merging
   recovers a longer insert per read, which means more splice evidence and a
   less pessimistic Illumina baseline.
2. **Adapter removal by geometry rather than by sequence matching.** The R1
   primer sits *before* the CB, so on a short fragment `rc(R2)`'s 5'
   read-through is primer and R1's 3' read-through is adapter — and both fall
   outside the overlap, so the merge drops them without ever consulting an
   adapter sequence. Each mate is independent evidence of where the fragment
   ends. A Hamming-seeded trimmer can only approximate that and carries a
   false-trim rate; this carries none, and it still works when the adapter
   itself was misread. Today the converter cannot detect R2 read-through at all.
3. **Insert-size distribution as a deliverable.** The merge naturally yields an
   insert-length histogram. Illumina insert size next to SBX fragment length is
   a slide that does not exist yet.

## Scope decisions (locked)

| Decision | Value | Rationale |
|---|---|---|
| Read geometry | **(a) long R1** — R1 sequenced past the 28 bp CB+UMI block into cDNA | Confirmed. R1 insert and `rc(R2)` genuinely overlap on short fragments. |
| Default | `--overlap` **off** | Existing CAT02 evaluation must stay reproducible bit-for-bit. On only for the long-read experiment. |
| YC tag | Emitted whenever `--overlap` is on | Presence/absence otherwise conflates "no overlap found" with "tool not run". |
| Output format | FASTQ header tags, as today | No BAM writing needed — see "Why this needs no BAM code". |

## Background: the geometry

**Merge the raw mates, before any trimming.** The RT primer's dT block *is* the
polyA complement — there is no separate templated polyT on top of it — so the
molecule R1 reads and the molecule `rc(R2)` reads are the same sequence:

```
molecule = CB(16) | UMI(12) | T(30) | rc(cDNA)
R1       = molecule                        | adapter        (truncated)
R2       = rc(molecule)                    | rc(R1 primer)  (truncated)
rc(R2)   = R1 primer | CB | UMI | T(30) | rc(cDNA)
```

This is why the merge must run **before** CB/UMI/dT are stripped: strip them from
R1 first and the mates no longer describe the same span, so a gapless alignment
cannot exist. Merging first also buys duplex confirmation of the barcode itself
on short fragments, where R2 reads back through the UMI and CB — the CB+UMI is
also 28 bp of complex, per-read-unique seed material.

R1 can still be *interior* to `rc(R2)`: whenever R1 is sequenced shorter than R2
and R2 spans the whole molecule, `rc(R2)` extends past R1 on both sides. Hence
index-the-longer / probe-the-shorter.

Whether they *reach* each other depends on insert size. With ~150 bp R1 the
usable insert portion is roughly `150 - 28 - dT_run ≈ 90 bp`; against a 150 bp
R2 that merges any fragment under ~240 bp. Typical 10x inserts run longer than
that, so **expect a modest merge rate, not a majority**. That is fine — the
histogram is informative either way, and it is an empirical question the metrics
answer directly. Do not treat a low merge fraction as a bug.

## Reference implementations

- **`hts_Overlapper`** — `~/projects/src/github.com-s4hts/HTStream/hts_Overlapper/src/hts_Overlapper.h`
  plus `common/src/utils.cpp` (`readOneMap`, `threshold_mismatches`,
  `setDefaultParamsOverlapping`). Apache-2.0; a port needs attribution, not a
  license change.
- **XOOS `yc_decode`** — `Roche-AXELIOS/XOOS`, `common/yc_decode/src/xoos/yc-decode/`.
  **XOOS Public License**, which restricts the work and its derivatives to data
  originating from Roche AXELIOS instruments. We are encoding YC for *Illumina*
  data, so: **implement from the published specification, do not translate
  `yc-decoder.cpp`.** The spec (GitBook docs) is documentation, not the licensed
  Work. Attribute the tables to the public docs in a comment. Worth an OSPO
  glance if sbxpress is ever distributed beyond UCSF.

---

## Phase 1 — `src/overlap.rs` (merge engine)

Pure, no I/O, no allocation on the reject path. Public surface:

```rust
pub struct OverlapCfg {
    pub kmer: usize,                    // 8
    pub kmer_offset: usize,             // 1
    pub check_lengths: usize,           // 20
    pub min_overlap: usize,             // 20  (NOT HTStream's 8 — see below)
    pub max_mismatch: usize,            // 100
    pub max_mismatch_density: f64,      // 0.25
    pub qual_cap: u8,                   // 40
}

pub enum InsertClass { Short, Medium, Long }

pub struct Merged {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub class: InsertClass,
    pub overlap_len: usize,
    pub left_overhang: usize,        // from R1; always 0 in 10x geometry
    pub right_overhang: usize,       // from rc(R2)
    pub discordant: usize,
    pub yc: Option<YcTag>,           // built here; the duplex slices are in hand
    pub yc_stats: EncodeStats,
}

pub fn overlap(r1: &[u8], r1q: &[u8], r2: &[u8], r2q: &[u8], cfg: &OverlapCfg)
    -> Option<Merged>;
```

Algorithm — port of `hts_Overlapper`, with the deviations below:

1. **Do not swap the merge frame.** Choose which mate to *index* independently
   instead: index the longer, probe the shorter. **Corrected during
   implementation** — the original plan said "no swap at all", which is wrong.
   In 10x geometry R1 loses 28 bp of CB+UMI plus the whole dT run, so its insert
   is routinely *contained* within `rc(R2)`, with overhangs on both sides;
   probing `rc(R2)`'s ends then finds nothing and every medium fragment silently
   fails to merge. Only the shorter read is guaranteed to have an end inside the
   overlap. The merge itself stays in R1's frame, which is what `YC` needs.
2. Index the chosen mate's k-mers at stride `kmer_offset`. **Deviation:** 2-bit encode into a
   `u64` (k ≤ 32) the way [`barcode::encode16`](src/barcode.rs) already does, and
   store in a sorted `Vec<(u64, u32)>` with binary search, or a
   `rustc_hash::FxHashMap<u64, SmallVec<[u32; 2]>>`. No `String` allocation per
   k-mer. Skip any k-mer containing `N`.
3. Probe the *other* mate's k-mers in its first and last `check_lengths` windows
   only, mapping each hit back into `(R1, rc(R2))` coordinates before merging.
   **Deviation:** guard the second loop — HTStream computes
   `seq2.length() - (checkLengths + kmer)` unguarded and underflows on short
   mates.
4. Verify each seed gaplessly. **Deviation:** use existing
   [`anchors::hamming_leq`](src/anchors.rs:53) instead of porting
   `threshold_mismatches`, which decrements before testing and underflows
   `size_t` when the mismatch budget is 0.
5. Merge. Agreeing bases: `min(q1 + q2, qual_cap)`. Disagreeing bases: take the
   higher-quality base with `max(1, |q1 - q2|)`. **Deviation:** HTStream writes
   `std::min(q1 - q2 + offset, 1 + offset)`, which is almost certainly a `min`/`max`
   slip — as written every discordant base collapses to Q1 regardless of the
   quality gap. Build the `YC` tag here too — the duplex slices are already in
   hand, so nothing has to re-walk the alignment.
6. Classify `Short` / `Medium` / `Long`. **Not** by HTStream's rule (merged
   length vs. longest input read) — the mates are not symmetric here, so that
   misclassifies. Use the geometric signal the classification actually means:
   read-through happened iff *either* mate ran past the end of the molecule —
   `reach1 > reach2` (R1 into adapter) **or** `start2 > 0` (rc(R2) began before
   R1's first base, i.e. R2 read back into the R1 primer). Checking only the
   first misses the asymmetric case: a 100 bp R1 that never reaches adapter
   paired with a 150 bp R2 that reads clean through the CB.

**`min_overlap` = 20, not 8.** A false merge here produces a chimeric read that
maps to the wrong place — strictly worse than not merging. HTStream's default of
8 is tuned for generic amplicon data; we have no such tolerance and we do have
plenty of reads.

**Homopolymer masking is mandatory.** The dT run sits immediately 5' of the R1
cDNA and any polyA/polyT in R2 will seed spuriously against it. Two guards:
advance past the run with [`extract::dt_run_end`](src/extract.rs:52) before
indexing (Phase 3), and drop k-mers that are a single repeated base from the
seed index. [`homopolymer::longest_run`](src/homopolymer.rs) is available if a
looser mask is wanted.

### Tests (Phase 1)

Perfect overlap; one mismatch inside budget; mismatch over budget rejects;
overlap shorter than `min_overlap` rejects; short insert with adapter overhangs;
`N`-containing seeds; R2 longer than R1 (swap path); empty and 1 bp inputs; the
underflow case that breaks HTStream (mate shorter than `check_lengths + kmer`).

---

## Phase 2 — `src/yc.rs` (YC tag)

Implement from the spec, cited to
`https://roche-axelios.gitbook.io/xoos/latest/tutorials/interpreting-sbx-reads/interpreting-sbx-d-reads`.

Format: `{left_tail}{delimiter}{duplex}{delimiter}{right_tail}`, always exactly
two delimiters. `+` = that tail came from R1, `-` = from R2; `+` is the default
when a tail is empty. Tails are bare integers (simplex base counts). The duplex
segment alternates integers (runs of concordant bases, never zero, no leading
zeros) with single-letter discordance codes; no integer appears between two
adjacent codes.

Code table, indexed `code - 'A'`:

```
index:      ABCDEFGHIJKLMNOPQRSTUVWXYZ
R1 base:    ACCGTTGTAGG-A--C-ACT-GATC-
R2 base:    AACACGGA--TAC-G-CGGT-CT-TT
complement: TKGYRMCWXPBZF-QJOEVA-SHIDL
```

- Uppercase = consensus has a base; op consumes **1** consensus position.
- Lowercase = consensus has a *gap*; op consumes **0** consensus positions, and
  the reconstructed base takes quality `!` (Q0).
- `A`/`C`/`G`/`T` are identity slots, never emitted. `N`/`U`/`n`/`u` unsupported.

Surface:

```rust
pub struct YcTag { /* left_overhang, duplex_ops, right_overhang, *_is_r1 */ }

impl YcTag {
    pub fn parse(s: &[u8]) -> Result<Self, YcError>;
    pub fn to_string(&self) -> String;
    pub fn reverse_complement(&mut self);
    pub fn trim(&mut self, left: usize, right: usize);
    pub fn sequence_length(&self) -> usize;
    pub fn is_duplex(&self) -> bool;     // duplex_ops non-empty
    pub fn is_simplex(&self) -> bool;    // empty tag counts as simplex
    pub fn count_discordant(&self) -> usize;
    pub fn base_types(&self) -> Vec<BaseType>;  // Concordant | Simplex | Discordant
}

pub fn encode(merged: &Merged) -> YcTag;
```

**Scope split:** the *encoder* only ever needs the 12 substitution codes, because
the overlap is gapless — no `I/P/J/X`, no `L/Q/O/Z`, no lowercase. The *parser*
must handle the full set anyway, both so round-trip tests are meaningful and so
we can ingest real SBX-D tags later.

The discordance code is a pure function of `(r1_base, r2_base)` and is
independent of which base the consensus chose. That is what makes the tag
lossless regardless of our quality-based tie-breaking.

### `YC` is lossless with respect to the merged construct, not the raw mates.

Bases the merge deliberately discards — adapter read-through, the polyT, barcode
read-through — are gone and the tag does not describe them. The invariants that
do hold: R1 contributes a prefix of its insert, `rc(R2)` contributes a suffix.

### Orientation — the one thing that must not go wrong

We build the tag in the merged (R1-insert) frame, then reverse-complement it
together with the sequence when we flip to R2 orientation for output. A YC tag
cannot be reversed by reversing the string: tails swap ends *and* every code must
go through the `complement` row.

Separately, note for any future consumer: alignment does not rewrite the tag, so
a read mapped to the minus strand has a revcomp'd `SEQ` and an un-revcomp'd `YC`.
Anything reading YC off a BAM must call `reverse_complement()` itself when
`FLAG & 16`. Document this loudly; it is the classic way to get a silently
wrong answer.

### Tests (Phase 2)

The four published worked examples as fixtures — all verified to reconcile
against the code table above:

| Tag | R1 | R2 |
|---|---|---|
| `4+13-3` | `CATGACGTACGGTCATG---` | `----ACGTACGGTCATGGCG` |
| `+1WK1I1L+` | `CAGAAG-` | `CTTA-GA` |
| `+1KM3L+2` | `CGACGC-GT` | `CTCCGCA--` |
| `9++` | `CGACGCACG` | `---------` |

Plus: `parse(to_string(t)) == t`; `reverse_complement` twice is identity;
`sequence_length()` equals the merged sequence length for every tag we emit
(cheap invariant, assert it in debug builds at the emit site); `9++`-style
all-simplex round trip; malformed inputs rejected (one delimiter, three
delimiters, leading zeros, zero-length runs, adjacent codes with an integer
between them, `N`/`U` codes).

**Property test:** generate a random gapless R1/R2 alignment, encode, decode,
assert both original strands come back byte-identical. This is the test that
actually proves the module.

---

## Phase 3 — wire into `sbxpress-illumina`

In [`convert_pair`](src/bin/sbxpress-illumina.rs:107), **before** CB/UMI
extraction — the merge runs on the raw mates:

1. `overlap::overlap_in(&mut scratch, r1_seq, r1_qual, r2_seq, r2_qual, &cfg)`.
2. On `Some(merged)`, the merged read *is* the molecule, so everything the
   converter already does now runs against it with duplex-backed qualities:
   CB at `merged.seq[..16]`, UMI at `[16..28]`, and
   `cdna_start = extract::dt_run_end(&merged.seq, 28)`.
3. Trim the tag to the emitted window — `tag.trim(cdna_start, adapter_trimmed)` —
   then `reverse_complement()` it alongside the sequence. This is what
   `YcTag::trim` exists for.
4. On `Some(merged)`, also:
   - **Do not** run [`extract::find_adapter`](src/extract.rs:94) over a merged
     read. The geometry has already removed read-through at both ends using the
     mates as independent evidence; layering a Hamming-seeded matcher on top can
     only contribute false trims. `find_adapter` stays on the unmerged (`Long`)
     path, where there is only one observation and matching is all there is.
   - `merged.yc` → `trim(cdna_start, 0)` → `reverse_complement()`;
   - emit `rc(merged.seq)` with reversed quality, so output stays in **R2
     orientation** and every downstream stage is untouched. [`lib.rs:682`](src/lib.rs:682)
     documents that `annotate_bam`'s strand logic depends on this convention.
5. On `None`: emit R2 exactly as today, plus a fully-simplex YC tag
   (`{len}-+`, left tail from R2).

Tag emission joins the existing `CR/CB/UR/UB/CY/UY` list on the FASTQ header
line:

```
@{id}\tCR:Z:..\tCB:Z:..\tUR:Z:..\tUB:Z:..\tCY:Z:..\tUY:Z:..\tYC:Z:{tag}
```

### Why this needs no BAM code

The converter writes tagged FASTQ consumed by `minimap2 -y`, which lifts the
whitespace-delimited comment into SAM aux fields verbatim. Appending
`\tYC:Z:...` is all that is required for the tag to land in the BAM. No
rust-htslib writing, no changes to `annotate_bam` — YC is pure passthrough.

### CLI

```
--overlap                          (default false)
--overlap-min-overlap    <N>       (default 20)
--overlap-kmer           <N>       (default 8)
--overlap-kmer-offset    <N>       (default 1)
--overlap-check-lengths  <N>       (default 20)
--overlap-max-mismatch   <N>       (default 100)
--overlap-max-mismatch-density <F> (default 0.25)
--no-emit-yc                       (YC on by default when --overlap is on)
```

Names deliberately mirror HTStream's so anyone who knows `hts_Overlapper` can
read them.

---

## Phase 4 — metrics and docs

Add to `Metrics` / `PairMetrics` in the converter:

```
overlap_enabled: bool
overlap_attempted: u64
overlap_short: u64          // read-through present on either mate
overlap_medium: u64         // insert >= read length
overlap_long: u64           // no overlap found; emitted as R2 alone
r1_insert_too_short: u64
merged_bases_total: u64
discordant_bases_total: u64
insert_length_histogram: BTreeMap<usize, u64>   // cDNA only: merged len - cdna_start
```

The histogram mirrors HTStream's `overlap_histogram` and is the presentation
artifact. `discordant_bases_total / merged_bases_total` is a free empirical
per-base disagreement rate for the Illumina arm — directly comparable to the
SBX-D discordance rate from real `YC` tags.

Docs: a section in `TECHNICAL_DESCRIPTION.md` covering the geometry, the merge
rules, and the YC encoding (including the minus-strand caveat); a note in
`MANUAL.md` for the new flags; `CHANGELOG.md` entry.

---

## Phase 5 — validation

1. **Cross-validate against `hts_Overlapper` directly.** HTStream is already
   checked out at `~/projects/src/github.com-s4hts/HTStream`. Build it, feed the
   same R1-insert/R2 pairs to both implementations, and diff the merged
   sequences. Expect near-identical output modulo the three deliberate deviations
   (discordant quality, the two underflow guards) — investigate any sequence-level
   difference. This is much stronger evidence than unit tests alone.
2. **YC round-trip on real data.** If any AXELIOS BAM is reachable, parse its YC
   tags and assert `base_types()` agrees position-for-position with the QUAL
   string under the documented intra-molecular encoding (`H` = Q39 concordant,
   `&` = Q5 discordant, `7` = Q22 simplex). Free correctness check against
   Roche's own encoder.
3. **Regression:** re-run `scripts/run_mapping_evaluation.sh` **without**
   `--overlap` and confirm byte-identical output against the committed CAT02
   results. The flag defaulting to off has to actually mean something.
4. **With `--overlap`:** compare mapping rate, mean mapped length, splice-junction
   evidence, and final cell/UMI counts against the no-overlap run.

---

## Risks

| Risk | Mitigation |
|---|---|
| Low merge rate makes the work look pointless | Expected on standard insert sizes; the histogram is the deliverable regardless. Set expectations before running. |
| False merges create chimeric reads | `min_overlap` 20, homopolymer-masked seeds, density cap. Report `overlap_short` separately — an implausible spike there means false merging. |
| dT run bleeds into the insert and seeds spuriously | `dt_run_end` skip plus homopolymer k-mer masking. |
| Orientation flip done wrong | Emit `rc(merged)`; revcomp the YC tag with it; assert `yc.sequence_length() == merged.len()` at the emit site. |
| XOOS license contamination | Implement from spec, never translate `yc-decoder.cpp`. Attribution comment citing the public docs. |
| YC on non-Roche data is semantically novel | Honest: it *is* two strands of one duplex fragment, same as SBX-D. Record the provenance in an `@PG` line so a merged BAM is never ambiguous about which arm produced a tag. |

## Status

* **Phase 1 — done.** `src/overlap.rs`, 29 tests.
* **Phase 2 — done.** `src/yc.rs`, 26 tests.
* Phases 3–5 outstanding. Full suite: 120 passing, clippy and `cargo fmt` clean.

## Estimate

Phase 1 ~250 lines + tests · Phase 2 ~300 lines + tests · Phase 3 ~120 lines ·
Phases 4–5 mostly docs and scripting. Roughly two focused days including
cross-validation.

## Open question

Should the merged read also carry SBX-D's flat intra-molecular quality encoding
(Q39/Q5/Q22) instead of real computed qualities, for strict comparability with
the SBX arm? Recommendation: **no by default** — real merged qualities are more
informative and are what mappers and variant callers want. Revisit as an
optional flag only if a specific comparison demands it.
