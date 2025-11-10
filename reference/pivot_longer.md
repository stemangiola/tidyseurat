# Pivot data from wide to long

`pivot_longer()` "lengthens" data, increasing the number of rows and
decreasing the number of columns. The inverse transformation is
[`pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)

Learn more in
[`vignette("pivot")`](https://tidyr.tidyverse.org/articles/pivot.html).

## Usage

``` r
# S3 method for class 'Seurat'
pivot_longer(
  data,
  cols,
  names_to = "name",
  names_prefix = NULL,
  names_sep = NULL,
  names_pattern = NULL,
  names_ptypes = NULL,
  names_transform = NULL,
  names_repair = "check_unique",
  values_to = "value",
  values_drop_na = FALSE,
  values_ptypes = NULL,
  values_transform = NULL,
  ...
)
```

## Arguments

- data:

  A data frame to pivot.

- cols:

  \<[`tidy-select`](https://tidyr.tidyverse.org/reference/tidyr_tidy_select.html)\>
  Columns to pivot into longer format.

- names_to:

  A character vector specifying the new column or columns to create from
  the information stored in the column names of `data` specified by
  `cols`.

  - If length 0, or if `NULL` is supplied, no columns will be created.

  - If length 1, a single column will be created which will contain the
    column names specified by `cols`.

  - If length \>1, multiple columns will be created. In this case, one
    of `names_sep` or `names_pattern` must be supplied to specify how
    the column names should be split. There are also two additional
    character values you can take advantage of:

    - `NA` will discard the corresponding component of the column name.

    - `".value"` indicates that the corresponding component of the
      column name defines the name of the output column containing the
      cell values, overriding `values_to` entirely.

- names_prefix:

  A regular expression used to remove matching text from the start of
  each variable name.

- names_sep, names_pattern:

  If `names_to` contains multiple values, these arguments control how
  the column name is broken up.

  `names_sep` takes the same specification as
  [`separate()`](https://tidyr.tidyverse.org/reference/separate.html),
  and can either be a numeric vector (specifying positions to break on),
  or a single string (specifying a regular expression to split on).

  `names_pattern` takes the same specification as
  [`extract()`](https://tidyr.tidyverse.org/reference/extract.html), a
  regular expression containing matching groups (`()`).

  If these arguments do not give you enough control, use
  [`pivot_longer_spec()`](https://tidyr.tidyverse.org/reference/pivot_longer_spec.html)
  to create a spec object and process manually as needed.

- names_ptypes, values_ptypes:

  Optionally, a list of column name-prototype pairs. Alternatively, a
  single empty prototype can be supplied, which will be applied to all
  columns. A prototype (or ptype for short) is a zero-length vector
  (like [`integer()`](https://rdrr.io/r/base/integer.html) or
  [`numeric()`](https://rdrr.io/r/base/numeric.html)) that defines the
  type, class, and attributes of a vector. Use these arguments if you
  want to confirm that the created columns are the types that you
  expect. Note that if you want to change (instead of confirm) the types
  of specific columns, you should use `names_transform` or
  `values_transform` instead.

- names_transform, values_transform:

  Optionally, a list of column name-function pairs. Alternatively, a
  single function can be supplied, which will be applied to all columns.
  Use these arguments if you need to change the types of specific
  columns. For example, `names_transform = list(week = as.integer)`
  would convert a character variable called `week` to an integer.

  If not specified, the type of the columns generated from `names_to`
  will be character, and the type of the variables generated from
  `values_to` will be the common type of the input columns used to
  generate them.

- names_repair:

  What happens if the output has invalid column names? The default,
  `"check_unique"` is to error if the columns are duplicated. Use
  `"minimal"` to allow duplicates in the output, or `"unique"` to
  de-duplicated by adding numeric suffixes. See
  [`vctrs::vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html)
  for more options.

- values_to:

  A string specifying the name of the column to create from the data
  stored in cell values. If `names_to` is a character containing the
  special `.value` sentinel, this value will be ignored, and the name of
  the value column will be derived from part of the existing column
  names.

- values_drop_na:

  If `TRUE`, will drop rows that contain only `NA`s in the `value_to`
  column. This effectively converts explicit missing values to implicit
  missing values, and should generally be used only when missing values
  in `data` were created by its structure.

- ...:

  Additional arguments passed on to methods.

## Value

\`tidyseurat\`

## Details

`pivot_longer()` is an updated approach to
[`gather()`](https://tidyr.tidyverse.org/reference/gather.html),
designed to be both simpler to use and to handle more use cases. We
recommend you use `pivot_longer()` for new code;
[`gather()`](https://tidyr.tidyverse.org/reference/gather.html) isn't
going away but is no longer under active development.

## Examples

``` r
data(pbmc_small)
pbmc_small |> pivot_longer(
  cols=c(orig.ident, groups),
  names_to="name", values_to="value")
#> tidyseurat says: A data frame is returned for independent data analysis.
#> # A tibble: 160 × 29
#>    .cell     nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents RNA_snn_res.1
#>    <chr>          <dbl>        <int> <fct>           <fct>         <fct>        
#>  1 ATGCCAGA…         70           47 0               A             0            
#>  2 ATGCCAGA…         70           47 0               A             0            
#>  3 CATGGCCT…         85           52 0               A             0            
#>  4 CATGGCCT…         85           52 0               A             0            
#>  5 GAACCTGA…         87           50 1               B             0            
#>  6 GAACCTGA…         87           50 1               B             0            
#>  7 TGACTGGA…        127           56 0               A             0            
#>  8 TGACTGGA…        127           56 0               A             0            
#>  9 AGTCAGAC…        173           53 0               A             0            
#> 10 AGTCAGAC…        173           53 0               A             0            
#> # ℹ 150 more rows
#> # ℹ 23 more variables: PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>, PC_4 <dbl>,
#> #   PC_5 <dbl>, PC_6 <dbl>, PC_7 <dbl>, PC_8 <dbl>, PC_9 <dbl>, PC_10 <dbl>,
#> #   PC_11 <dbl>, PC_12 <dbl>, PC_13 <dbl>, PC_14 <dbl>, PC_15 <dbl>,
#> #   PC_16 <dbl>, PC_17 <dbl>, PC_18 <dbl>, PC_19 <dbl>, tSNE_1 <dbl>,
#> #   tSNE_2 <dbl>, name <chr>, value <chr>
```
