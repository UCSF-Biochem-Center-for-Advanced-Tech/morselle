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

## The gate whitelist: `737K-arc-v1.txt.gz` (Multiome ARC)

**This data is 10x Single Cell Multiome (ATAC + GEX), GEX side**, so the cell
barcodes come from the **ARC** gel-bead list `737K-arc-v1.txt.gz` (736,320
barcodes, single column) — NOT any 3' GEX list. M4 loads this via `--whitelist`.

Verified 2026-06-24 against Roche's own known-good CBs: `737K-arc-v1` = **73.3%
as-is** (0% reverse-complemented), vs <2% for every 3' GEX list (`3M-february-2018`
v3, `3M-3pgex-may-2023`). Get it from the **Cell Ranger ARC** distribution:
`cellranger-arc-<version>/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz`.
Match orientation is **as-is** (the ATAC barcode is the reverse complement).

The 3' GEX lists (`3M-february-2018*`, `3M-3pgex-may-2023*`, `9K-LT-*`) do not
apply to this data and are kept only for reference.

## Notes

- The `_NXT` whitelists are translation tables (two 16 bp columns; `col1 != col2`
  on every row): column 1 is the barcode as sequenced, column 2 the reference
  barcode. They are kept for reference / possible translation use, but are not the
  GEX gate set.
- The relevant file for the current SBX runs (10x 3' GEX v3/v3.1) is
  **`3M-february-2018_NXT.txt.gz`**.
