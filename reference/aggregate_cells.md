# Aggregate cells

Combine cells into groups based on shared variables and aggregate
feature counts.

## Usage

``` r
# S4 method for class 'Seurat'
aggregate_cells(
  .data,
  .sample = NULL,
  slot = "data",
  assays = NULL,
  aggregation_function = Matrix::rowSums,
  ...
)
```

## Arguments

- .data:

  A tidyseurat object

- .sample:

  A vector of variables by which cells are aggregated

- slot:

  The slot to which the function is applied

- assays:

  The assay to which the function is applied

- aggregation_function:

  The method of cell-feature value aggregation

- ...:

  Used for future extendibility

## Value

A tibble object

## Examples

``` r
data(pbmc_small)
pbmc_small_pseudo_bulk <- pbmc_small |>
  aggregate_cells(c(groups, letter.idents), assays="RNA")
#> Joining with `by = join_by(letter.idents, groups)`
```
