# Separate a character column into multiple columns with a regular expression or numeric locations

**\[superseded\]**

`separate()` has been superseded in favour of
[`separate_wider_position()`](https://tidyr.tidyverse.org/reference/separate_wider_delim.html)
and
[`separate_wider_delim()`](https://tidyr.tidyverse.org/reference/separate_wider_delim.html)
because the two functions make the two uses more obvious, the API is
more polished, and the handling of problems is better. Superseded
functions will not go away, but will only receive critical bug fixes.

Given either a regular expression or a vector of character positions,
`separate()` turns a single character column into multiple columns.

## Usage

``` r
# S3 method for class 'Seurat'
separate(
  data,
  col,
  into,
  sep = "[^[:alnum:]]+",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn",
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

- sep:

  Separator between columns.

  If character, `sep` is interpreted as a regular expression. The
  default value is a regular expression that matches any sequence of
  non-alphanumeric values.

  If numeric, `sep` is interpreted as character positions to split at.
  Positive values start at 1 at the far-left of the string; negative
  value start at -1 at the far-right of the string. The length of `sep`
  should be one less than `into`.

- remove:

  If `TRUE`, remove input column from output data frame.

- convert:

  If `TRUE`, will run
  [`type.convert()`](https://rdrr.io/r/utils/type.convert.html) with
  `as.is = TRUE` on new columns. This is useful if the component columns
  are integer, numeric or logical.

  NB: this will cause string `"NA"`s to be converted to `NA`s.

- extra:

  If `sep` is a character vector, this controls what happens when there
  are too many pieces. There are three valid options:

  - `"warn"` (the default): emit a warning and drop extra values.

  - `"drop"`: drop any extra values without a warning.

  - `"merge"`: only splits at most `length(into)` times

- fill:

  If `sep` is a character vector, this controls what happens when there
  are not enough pieces. There are three valid options:

  - `"warn"` (the default): emit a warning and fill from the right

  - `"right"`: fill with missing values on the right

  - `"left"`: fill with missing values on the left

- ...:

  Additional arguments passed on to methods.

## Value

\`tidyseurat\`

## See also

[`unite()`](https://tidyr.tidyverse.org/reference/unite.html), the
complement,
[`extract()`](https://tidyr.tidyverse.org/reference/extract.html) which
uses regular expression capturing groups.

## Examples

``` r
data(pbmc_small)
un <- pbmc_small |> unite("new_col", c(orig.ident, groups))
un |> separate(new_col, c("orig.ident", "groups"))
#> # A Seurat-tibble abstraction: 80 × 15
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident groups nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
#>    <chr> <chr>      <chr>       <dbl>        <int> <fct>           <fct>        
#>  1 ATGC… SeuratPro… g2             70           47 0               A            
#>  2 CATG… SeuratPro… g1             85           52 0               A            
#>  3 GAAC… SeuratPro… g2             87           50 1               B            
#>  4 TGAC… SeuratPro… g2            127           56 0               A            
#>  5 AGTC… SeuratPro… g2            173           53 0               A            
#>  6 TCTG… SeuratPro… g1             70           48 0               A            
#>  7 TGGT… SeuratPro… g1             64           36 0               A            
#>  8 GCAG… SeuratPro… g1             72           45 0               A            
#>  9 GATA… SeuratPro… g1             52           36 0               A            
#> 10 AATG… SeuratPro… g1            100           41 0               A            
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
