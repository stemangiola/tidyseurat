# Nest rows into a list-column of data frames

Nesting creates a list-column of data frames; unnesting flattens it back
out into regular columns. Nesting is implicitly a summarising operation:
you get one row for each group defined by the non-nested columns. This
is useful in conjunction with other summaries that work with whole
datasets, most notably models.

Learn more in
[`vignette("nest")`](https://tidyr.tidyverse.org/articles/nest.html).

## Usage

``` r
# S3 method for class 'Seurat'
nest(.data, ..., .names_sep = NULL)
```

## Arguments

- .data:

  A data frame.

- ...:

  \<[`tidy-select`](https://tidyr.tidyverse.org/reference/tidyr_tidy_select.html)\>
  Columns to nest; these will appear in the inner data frames.

  Specified using name-variable pairs of the form
  `new_col = c(col1, col2, col3)`. The right hand side can be any valid
  tidyselect expression.

  If not supplied, then `...` is derived as all columns *not* selected
  by `.by`, and will use the column name from `.key`.

  **\[deprecated\]**: previously you could write `df %>% nest(x, y, z)`.
  Convert to `df %>% nest(data = c(x, y, z))`.

- .names_sep:

  If `NULL`, the default, the inner names will come from the former
  outer names. If a string, the new inner names will use the outer names
  with `names_sep` automatically stripped. This makes `names_sep`
  roughly symmetric between nesting and unnesting.

## Value

\`tidyseurat_nested\`

## Details

If neither `...` nor `.by` are supplied, `nest()` will nest all
variables, and will use the column name supplied through `.key`.

## New syntax

tidyr 1.0.0 introduced a new syntax for `nest()` and
[`unnest()`](unnest.md) that's designed to be more similar to other
functions. Converting to the new syntax should be straightforward
(guided by the message you'll receive) but if you just need to run an
old analysis, you can easily revert to the previous behaviour using
[`nest_legacy()`](https://tidyr.tidyverse.org/reference/nest_legacy.html)
and
[`unnest_legacy()`](https://tidyr.tidyverse.org/reference/nest_legacy.html)
as follows:

    library(tidyr)
    nest <- nest_legacy
    unnest <- unnest_legacy

## Grouped data frames

`df %>% nest(data = c(x, y))` specifies the columns to be nested; i.e.
the columns that will appear in the inner data frame.
`df %>% nest(.by = c(x, y))` specifies the columns to nest *by*; i.e.
the columns that will remain in the outer data frame. An alternative way
to achieve the latter is to `nest()` a grouped data frame created by
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
The grouping variables remain in the outer data frame and the others are
nested. The result preserves the grouping of the input.

Variables supplied to `nest()` will override grouping variables so that
`df %>% group_by(x, y) %>% nest(data = !z)` will be equivalent to
`df %>% nest(data = !z)`.

You can't supply `.by` with a grouped data frame, as the groups already
represent what you are nesting by.

## Examples

``` r
data(pbmc_small)
pbmc_small |> 
    nest(data=-groups) |> 
    unnest(data)
#> # A Seurat-tibble abstraction: 80 × 8
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell        orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
#>    <chr>        <chr>           <dbl>        <int> <chr>           <chr>        
#>  1 ATGCCAGAACG… SeuratPro…         70           47 0               A            
#>  2 GAACCTGATGA… SeuratPro…         87           50 1               B            
#>  3 TGACTGGATTC… SeuratPro…        127           56 0               A            
#>  4 AGTCAGACTGC… SeuratPro…        173           53 0               A            
#>  5 AGGTCATGAGT… SeuratPro…         62           31 0               A            
#>  6 GGGTAACTCTA… SeuratPro…        101           41 0               A            
#>  7 CATGAGACACG… SeuratPro…         51           26 0               A            
#>  8 TACGCCACTCC… SeuratPro…         99           45 0               A            
#>  9 GTAAGCACTCA… SeuratPro…         67           33 0               A            
#> 10 TACATCACGCT… SeuratPro…        109           41 0               A            
#> # ℹ 70 more rows
#> # ℹ 2 more variables: RNA_snn_res.1 <chr>, groups <chr>
```
