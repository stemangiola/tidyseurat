# Sample n rows from a table

**\[superseded\]** `sample_n()` and `sample_frac()` have been superseded
in favour of
[`slice_sample()`](https://dplyr.tidyverse.org/reference/slice.html).
While they will not be deprecated in the near future, retirement means
that we will only perform critical bug fixes, so we recommend moving to
the newer alternative.

These functions were superseded because we realised it was more
convenient to have two mutually exclusive arguments to one function,
rather than two separate functions. This also made it to clean up a few
other smaller design issues with `sample_n()`/`sample_frac`:

- The connection to [`slice()`](slice.md) was not obvious.

- The name of the first argument, `tbl`, is inconsistent with other
  single table verbs which use `.data`.

- The `size` argument uses tidy evaluation, which is surprising and
  undocumented.

- It was easier to remove the deprecated `.env` argument.

- `...` was in a suboptimal position.

## Usage

``` r
# S3 method for class 'Seurat'
sample_n(tbl, size, replace = FALSE, weight = NULL, .env = NULL, ...)

# S3 method for class 'Seurat'
sample_frac(tbl, size = 1, replace = FALSE, weight = NULL, .env = NULL, ...)
```

## Arguments

- tbl:

  A data.frame.

- size:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  For `sample_n()`, the number of rows to select. For `sample_frac()`,
  the fraction of rows to select. If `tbl` is grouped, `size` applies to
  each group.

- replace:

  Sample with or without replacement?

- weight:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Sampling weights. This must evaluate to a vector of non-negative
  numbers the same length as the input. Weights are automatically
  standardised to sum to 1.

- .env:

  DEPRECATED.

- ...:

  ignored

## Examples

``` r
data(pbmc_small)
pbmc_small |> sample_n(50)
#> # A Seurat-tibble abstraction: 50 × 15
#> # Features=230 | Cells=50 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 AGAT… SeuratPro…        187           61 0               A             g2    
#>  2 TACA… SeuratPro…        108           44 0               A             g2    
#>  3 CATG… SeuratPro…         51           26 0               A             g2    
#>  4 GCAC… SeuratPro…        292           71 1               B             g2    
#>  5 CGTA… SeuratPro…        371           75 1               B             g1    
#>  6 TTAC… SeuratPro…        298           65 1               B             g1    
#>  7 ATAA… SeuratPro…         99           42 1               B             g2    
#>  8 TGGT… SeuratPro…         64           36 0               A             g1    
#>  9 GTTG… SeuratPro…        221           67 0               A             g2    
#> 10 GGCA… SeuratPro…        172           29 0               A             g1    
#> # ℹ 40 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
pbmc_small |> sample_frac(0.1)
#> # A Seurat-tibble abstraction: 8 × 15
#> # Features=230 | Cells=8 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 GATAG… SeuratPro…        328           72 1               B             g1    
#> 2 GGCAT… SeuratPro…        126           53 0               A             g1    
#> 3 ATGCC… SeuratPro…         70           47 0               A             g2    
#> 4 AGATA… SeuratPro…        187           61 0               A             g2    
#> 5 TACAA… SeuratPro…        108           44 0               A             g2    
#> 6 CATGA… SeuratPro…         51           26 0               A             g2    
#> 7 GCACT… SeuratPro…        292           71 1               B             g2    
#> 8 CGTAG… SeuratPro…        371           75 1               B             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
