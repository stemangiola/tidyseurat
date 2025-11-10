# Order rows using column values

`arrange()` orders the rows of a data frame by the values of selected
columns.

Unlike other dplyr verbs, `arrange()` largely ignores grouping; you need
to explicitly mention grouping variables (or use `.by_group = TRUE`) in
order to group by them, and functions of variables are evaluated once
per data frame, not once per group.

## Usage

``` r
# S3 method for class 'Seurat'
arrange(.data, ..., .by_group = FALSE)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Variables, or functions of variables. Use
  [`desc()`](https://dplyr.tidyverse.org/reference/desc.html) to sort a
  variable in descending order.

- .by_group:

  If `TRUE`, will sort first by grouping variable. Applies to grouped
  data frames only.

## Value

An object of the same type as `.data`. The output has the following
properties:

- All rows appear in the output, but (usually) in a different place.

- Columns are not modified.

- Groups are not modified.

- Data frame attributes are preserved.

## Details

### Missing values

Unlike base sorting with [`sort()`](https://rdrr.io/r/base/sort.html),
`NA` are:

- always sorted to the end for local data, even when wrapped with
  [`desc()`](https://dplyr.tidyverse.org/reference/desc.html).

- treated differently for remote data, depending on the backend.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`), plotly
([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## See also

Other single table verbs: [`mutate()`](mutate.md),
[`rename()`](rename.md), [`slice()`](slice.md),
[`summarise()`](summarise.md)

## Examples

``` r
data(pbmc_small)
pbmc_small |>
    arrange(nFeature_RNA)
#> Warning: `arrange()` was deprecated in tidyseurat 0.7.5.
#> ℹ tidyseurat says: arrange() is temporarly deprected as it is not clear that
#>   Seurat allows reordering of cells.
#> # A Seurat-tibble abstraction: 80 × 15
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 ATGC… SeuratPro…         70           47 0               A             g2    
#>  2 CATG… SeuratPro…         85           52 0               A             g1    
#>  3 GAAC… SeuratPro…         87           50 1               B             g2    
#>  4 TGAC… SeuratPro…        127           56 0               A             g2    
#>  5 AGTC… SeuratPro…        173           53 0               A             g2    
#>  6 TCTG… SeuratPro…         70           48 0               A             g1    
#>  7 TGGT… SeuratPro…         64           36 0               A             g1    
#>  8 GCAG… SeuratPro…         72           45 0               A             g1    
#>  9 GATA… SeuratPro…         52           36 0               A             g1    
#> 10 AATG… SeuratPro…        100           41 0               A             g1    
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
