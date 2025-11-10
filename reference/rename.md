# Rename columns

`rename()` changes the names of individual variables using
`new_name = old_name` syntax;
[`rename_with()`](https://dplyr.tidyverse.org/reference/rename.html)
renames columns using a function.

## Usage

``` r
# S3 method for class 'Seurat'
rename(.data, ...)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  For `rename()`:
  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Use `new_name = old_name` to rename selected variables.

  For
  [`rename_with()`](https://dplyr.tidyverse.org/reference/rename.html):
  additional arguments passed onto `.fn`.

## Value

An object of the same type as `.data`. The output has the following
properties:

- Rows are not affected.

- Column names are changed; column order is preserved.

- Data frame attributes are preserved.

- Groups are updated to reflect new names.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`), plotly
([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## See also

Other single table verbs: [`arrange()`](arrange.md),
[`mutate()`](mutate.md), [`slice()`](slice.md),
[`summarise()`](summarise.md)

## Examples

``` r
data(pbmc_small)
pbmc_small |> rename(s_score=nFeature_RNA)
#> # A Seurat-tibble abstraction: 80 × 15
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell      orig.ident nCount_RNA s_score RNA_snn_res.0.8 letter.idents groups
#>    <chr>      <fct>           <dbl>   <int> <fct>           <fct>         <chr> 
#>  1 ATGCCAGAA… SeuratPro…         70      47 0               A             g2    
#>  2 CATGGCCTG… SeuratPro…         85      52 0               A             g1    
#>  3 GAACCTGAT… SeuratPro…         87      50 1               B             g2    
#>  4 TGACTGGAT… SeuratPro…        127      56 0               A             g2    
#>  5 AGTCAGACT… SeuratPro…        173      53 0               A             g2    
#>  6 TCTGATACA… SeuratPro…         70      48 0               A             g1    
#>  7 TGGTATCTA… SeuratPro…         64      36 0               A             g1    
#>  8 GCAGCTCTG… SeuratPro…         72      45 0               A             g1    
#>  9 GATATAACA… SeuratPro…         52      36 0               A             g1    
#> 10 AATGTTGAC… SeuratPro…        100      41 0               A             g1    
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
