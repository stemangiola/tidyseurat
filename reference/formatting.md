# Printing tibbles

One of the main features of the `tbl_df` class is the printing:

- Tibbles only print as many rows and columns as fit on one screen,
  supplemented by a summary of the remaining rows and columns.

- Tibble reveals the type of each column, which keeps the user informed
  about whether a variable is, e.g., `<chr>` or `<fct>` (character
  versus factor). See `vignette("types")` for an overview of common type
  abbreviations.

Printing can be tweaked for a one-off call by calling `print()`
explicitly and setting arguments like `n` and `width`. More persistent
control is available by setting the options described in
[pillar::pillar_options](https://pillar.r-lib.org/reference/pillar_options.html).
See also `vignette("digits")` for a comparison to base options, and
`vignette("numbers")` that showcases
[`num()`](https://tibble.tidyverse.org/reference/num.html) and
[`char()`](https://tibble.tidyverse.org/reference/char.html) for
creating columns with custom formatting options.

As of tibble 3.1.0, printing is handled entirely by the pillar package.
If you implement a package that extends tibble, the printed output can
be customized in various ways. See
[`vignette("extending", package = "pillar")`](https://pillar.r-lib.org/articles/extending.html)
for details, and
[pillar::pillar_options](https://pillar.r-lib.org/reference/pillar_options.html)
for options that control the display in the console.

## Usage

``` r
# S3 method for class 'Seurat'
print(x, ..., n = NULL, width = NULL, n_extra = NULL)
```

## Arguments

- x:

  Object to format or print.

- ...:

  Passed on to
  [`tbl_format_setup()`](https://pillar.r-lib.org/reference/tbl_format_setup.html).

- n:

  Number of rows to show. If `NULL`, the default, will print all rows if
  less than the `print_max`
  [option](https://pillar.r-lib.org/reference/pillar_options.html).
  Otherwise, will print as many rows as specified by the `print_min`
  [option](https://pillar.r-lib.org/reference/pillar_options.html).

- width:

  Width of text output to generate. This defaults to `NULL`, which means
  use the `width`
  [option](https://pillar.r-lib.org/reference/pillar_options.html).

- n_extra:

  Number of extra columns to print abbreviated information for, if the
  width is too small for the entire tibble. If \`NULL\`, the default,
  will print information about at most \`tibble.max_extra_cols\` extra
  columns.

## Value

Prints a message to the console describing the contents of the
\`tidyseurat\`.

## Examples

``` r
data(pbmc_small)
print(pbmc_small)
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
