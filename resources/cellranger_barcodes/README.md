# Cell Ranger barcode whitelists

10x Genomics barcode whitelists / translation tables used by `sbx10x` for the
cell-barcode validate-gate (plan §4 step 5, M4). **These files are not committed**
— they are external reference data shipped with Cell Ranger. Populate this folder
locally before running with `--whitelist`.

## Where to get them

Bundled in the Cell Ranger distribution under:

```
cellranger-<version>/lib/python/cellranger/barcodes/
```

(plus the `translation/` subfolder for the 2-column `_NXT` files). Copy the files
below into this directory.

## Expected files

| File | Chemistry / role | Format | Rows |
|---|---|---|---|
| `3M-february-2018_NXT.txt.gz` | **3' GEX v3 / v3.1** — primary whitelist for this project | 2-col TSV: `sequenced<TAB>translated` | 3,686,400 |
| `3M-3pgex-may-2023_NXT.txt.gz` | 3' GEX 2023 chemistry | 2-col TSV | 3,686,400 |
| `9K-LT-march-2021.txt.gz` | 3' GEX low-throughput (LT) | 2-col TSV | 9,216 |
| `flex-96-b1.txt` | Flex 96-plex probe barcodes (batch 1) | 3-col: `bc1<TAB>bc2<TAB>id` | 96 |
| `flex-v2-384.txt` | Flex v2 384-plex probe barcodes | 3-col | 384 |
| `overhang.txt` | Overhang barcodes | 3-col | 4 |
| `probe-barcodes-fixed-rna-profiling-rna.txt` | Fixed-RNA profiling, RNA probe BCs | 3-col | 128 |
| `probe-barcodes-fixed-rna-profiling-rna-r1.txt` | …RNA, R1 variant | 3-col | 64 |
| `probe-barcodes-fixed-rna-profiling-ab.txt` | …antibody-capture probe BCs | 3-col | 128 |
| `probe-barcodes-fixed-rna-profiling-ab-r1.txt` | …antibody, R1 variant | 3-col | 64 |
| `probe-barcodes-fixed-rna-profiling-crispr.txt` | …CRISPR probe BCs | 3-col | 16 |

## ⚠️ Needed for M4: the standard GEX whitelist (not present)

The `_NXT` files here are **the wrong list for the GEX validate-gate**. Verified
2026-06-24: neither column contains the 3' v3 cell barcodes — even Roche's own
known-good R1 CBs (from `M2_GEX/`) hit these files at <1% (chance level). They are
NextSeq/translation tables (3,686,400 rows), not the 6,794,880-barcode GEX list.

**Action:** add the single-column `3M-february-2018.txt.gz` (6,794,880 barcodes,
one 16 bp barcode per line) from
`cellranger-<version>/lib/python/cellranger/barcodes/3M-february-2018.txt.gz`
(NOT the `translation/` subfolder). M4 loads that.

(Extraction was independently validated without a whitelist: our M1+M2 CBs match
Roche's CB on 87.6% of shared reads by read id.)

## Notes

- The `_NXT` whitelists are translation tables (two 16 bp columns; `col1 != col2`
  on every row): column 1 is the barcode as sequenced, column 2 the reference
  barcode. They are kept for reference / possible translation use, but are not the
  GEX gate set.
- The relevant file for the current SBX runs (10x 3' GEX v3/v3.1) is
  **`3M-february-2018_NXT.txt.gz`**.
