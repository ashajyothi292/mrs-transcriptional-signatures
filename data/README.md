# Data

`sample_design.csv` lists the 195 E-MTAB-14038 samples used and their tumour type.

`gene_sets/` holds the curated workbook and the GMT file built from it.

The expression data is not stored in this repository. Put it in `data/raw/`:

- **E-MTAB-14038**: the gene-level counts table, cut down to the 195 samples in
  `sample_design.csv` and saved as `data/raw/counts_matrix_195.csv` (gene symbols in the
  first column, one column per sample). Script 02 then writes `counts_matrix_195_fixed.csv`.
- **GSE85217**: downloaded automatically by script 04 into `data/raw/GSE85217/`.
- **Hallmark gene sets** (script 06 only): `h.all.v2024.1.Hs.symbols.gmt` from MSigDB,
  placed in `data/gene_sets/`.
