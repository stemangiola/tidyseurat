# Cell types of 80 PBMC single cells

A dataset containing the barcodes and cell types of 80 PBMC single
cells.

## Usage

``` r
data(cell_type_df)
```

## Format

A tibble containing 80 rows and 2 columns. Cells are a subsample of the
Peripheral Blood Mononuclear Cells (PBMC) dataset of 2,700 single cell.
Cell types were identified with SingleR.

- cell:

  cell identifier, barcode

- first.labels:

  cell type

## Source

<https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html>

## Value

\`tibble\`
