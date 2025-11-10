# Extract a character column into multiple columns using regular expression groups

**\[superseded\]**

`extract()` has been superseded in favour of
[`separate_wider_regex()`](https://tidyr.tidyverse.org/reference/separate_wider_delim.html)
because it has a more polished API and better handling of problems.
Superseded functions will not go away, but will only receive critical
bug fixes.

Given a regular expression with capturing groups, `extract()` turns each
group into a new column. If the groups don't match, or the input is NA,
the output will be NA.

## Usage

``` r
# S3 method for class 'Seurat'
extract(
  data,
  col,
  into,
  regex = "([[:alnum:]]+)",
  remove = TRUE,
  convert = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame.

- col:

  \<[`tidy-select`](https://tidyr.tidyverse.org/reference/tidyr_tidy_select.html)\>
  Column to expand.

- into:

  Names of new variables to create as character vector. Use `NA` to omit
  the variable in the output.

- regex:

  A string representing a regular expression used to extract the desired
  values. There should be one group (defined by `()`) for each element
  of `into`.

- remove:

  If `TRUE`, remove input column from output data frame.

- convert:

  If `TRUE`, will run
  [`type.convert()`](https://rdrr.io/r/utils/type.convert.html) with
  `as.is = TRUE` on new columns. This is useful if the component columns
  are integer, numeric or logical.

  NB: this will cause string `"NA"`s to be converted to `NA`s.

- ...:

  Additional arguments passed on to methods.

## Value

\`tidyseurat\`

## See also

[`separate()`](https://tidyr.tidyverse.org/reference/separate.html) to
split up by a separator.

## Examples

``` r
data(pbmc_small)
pbmc_small |>
  extract(groups, 
    into="g", 
    regex="g([0-9])", 
    convert=TRUE)
#> # A Seurat-tibble abstraction: 80 × 15
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents     g
#>    <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <int>
#>  1 ATGCC… SeuratPro…         70           47 0               A                 2
#>  2 CATGG… SeuratPro…         85           52 0               A                 1
#>  3 GAACC… SeuratPro…         87           50 1               B                 2
#>  4 TGACT… SeuratPro…        127           56 0               A                 2
#>  5 AGTCA… SeuratPro…        173           53 0               A                 2
#>  6 TCTGA… SeuratPro…         70           48 0               A                 1
#>  7 TGGTA… SeuratPro…         64           36 0               A                 1
#>  8 GCAGC… SeuratPro…         72           45 0               A                 1
#>  9 GATAT… SeuratPro…         52           36 0               A                 1
#> 10 AATGT… SeuratPro…        100           41 0               A                 1
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
