# Keep distinct/unique rows

Keep only unique/distinct rows from a data frame. This is similar to
[`unique.data.frame()`](https://rdrr.io/r/base/unique.html) but
considerably faster.

## Usage

``` r
# S3 method for class 'Seurat'
distinct(.data, ..., .keep_all = FALSE)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional variables to use when determining uniqueness. If there are
  multiple rows for a given combination of inputs, only the first row
  will be preserved. If omitted, will use all variables in the data
  frame.

- .keep_all:

  If `TRUE`, keep all variables in `.data`. If a combination of `...` is
  not distinct, this keeps the first row of values.

## Value

An object of the same type as `.data`. The output has the following
properties:

- Rows are a subset of the input but appear in the same order.

- Columns are not modified if `...` is empty or `.keep_all` is `TRUE`.
  Otherwise, `distinct()` first calls [`mutate()`](mutate.md) to create
  new columns.

- Groups are not modified.

- Data frame attributes are preserved.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`), plotly
([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## Examples

``` r
data("pbmc_small")
pbmc_small |> distinct(groups)
#> tidyseurat says: A data frame is returned for independent data analysis.
#> # A tibble: 2 × 1
#>   groups
#>   <chr> 
#> 1 g2    
#> 2 g1    
```
