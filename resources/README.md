# sbxpress resources

Reference-adjacent data distributed or staged for sbxpress lives here.

## `10x_barcodes/`

10x Genomics cell-barcode whitelists and related translation tables. Keep the
upstream filenames unchanged so chemistry provenance remains explicit. The
installer copies these files to `<prefix>/share/sbxpress/barcodes/`.

The current Multiome GEX workflow uses `737K-arc-v1.txt.gz`. Other assays must
select the whitelist appropriate to their chemistry; sbxpress does not infer a
chemistry merely from the presence of a file in this directory.

These files originate with 10x Genomics software/reference distributions and
must retain their applicable upstream licensing and redistribution terms.

## Gene manifests

Gene manifests are reference-build outputs rather than bundled resources. Create
one canonical `genes.mtx` beside each minimap2/GTF reference by passing
`--genes-mtx` to `sbxpress annotate`. Reuse that exact file with
`sbxpress matrix --features` for every sample and experiment built against the
reference.
