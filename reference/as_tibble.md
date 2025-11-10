# Coerce lists, matrices, and more to data frames

`as_tibble()` turns an existing object, such as a data frame or matrix,
into a so-called tibble, a data frame with class
[`tbl_df`](https://tibble.tidyverse.org/reference/tbl_df-class.html).
This is in contrast with
[`tibble()`](https://tibble.tidyverse.org/reference/tibble.html), which
builds a tibble from individual columns. `as_tibble()` is to
[`tibble()`](https://tibble.tidyverse.org/reference/tibble.html) as
[`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) is
to [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html).

`as_tibble()` is an S3 generic, with methods for:

- [`data.frame`](https://rdrr.io/r/base/data.frame.html): Thin wrapper
  around the `list` method that implements tibble's treatment of
  [rownames](https://tibble.tidyverse.org/reference/rownames.html).

- [`matrix`](https://rdrr.io/r/base/matrix.html),
  [`poly`](https://rdrr.io/r/stats/poly.html),
  [`ts`](https://rdrr.io/r/stats/ts.html),
  [`table`](https://rdrr.io/r/base/table.html)

- Default: Other inputs are first coerced with
  [`base::as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

`as_tibble_row()` converts a vector to a tibble with one row. If the
input is a list, all elements must have size one.

`as_tibble_col()` converts a vector to a tibble with one column.

## Usage

``` r
# S3 method for class 'Seurat'
as_tibble(
  x,
  ...,
  .name_repair = c("check_unique", "unique", "universal", "minimal"),
  rownames = NULL
)
```

## Arguments

- x:

  A data frame, list, matrix, or other object that could reasonably be
  coerced to a tibble.

- ...:

  Unused, for extensibility.

- .name_repair:

  Treatment of problematic column names:

  - `"minimal"`: No name repair or checks, beyond basic existence,

  - `"unique"`: Make sure names are unique and not empty,

  - `"check_unique"`: (default value), no name repair, but check they
    are `unique`,

  - `"universal"`: Make the names `unique` and syntactic

  - `"unique_quiet"`: Same as `"unique"`, but "quiet"

  - `"universal_quiet"`: Same as `"universal"`, but "quiet"

  - a function: apply custom name repair (e.g.,
    `.name_repair = make.names` for names in the style of base R).

  - A purrr-style anonymous function, see
    [`rlang::as_function()`](https://rlang.r-lib.org/reference/as_function.html)

  This argument is passed on as `repair` to
  [`vctrs::vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html).
  See there for more details on these terms and the strategies used to
  enforce them.

- rownames:

  How to treat existing row names of a data frame or matrix:

  - `NULL`: remove row names. This is the default.

  - `NA`: keep row names.

  - A string: the name of a new column. Existing rownames are
    transferred into this column and the `row.names` attribute is
    deleted. No name repair is applied to the new column name, even if
    `x` already contains a column of that name. Use
    `as_tibble(rownames_to_column(...))` to safeguard against this case.

  Read more in
  [rownames](https://tibble.tidyverse.org/reference/rownames.html).

## Value

\`tibble\`

## Row names

The default behavior is to silently remove row names.

New code should explicitly convert row names to a new column using the
`rownames` argument.

For existing code that relies on the retention of row names, call
`pkgconfig::set_config("tibble::rownames" = NA)` in your script or in
your package's [`.onLoad()`](https://rdrr.io/r/base/ns-hooks.html)
function.

## Life cycle

Using `as_tibble()` for vectors is superseded as of version 3.0.0,
prefer the more expressive `as_tibble_row()` and `as_tibble_col()`
variants for new code.

## See also

[`tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
constructs a tibble from individual columns.
[`enframe()`](https://tibble.tidyverse.org/reference/enframe.html)
converts a named vector to a tibble with a column of names and column of
values. Name repair is implemented using
[`vctrs::vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html).

## Examples

``` r
data(pbmc_small)
pbmc_small |> as_tibble()
#> # A tibble: 80 × 29
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
#> # ℹ 22 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, PC_6 <dbl>, PC_7 <dbl>, PC_8 <dbl>, PC_9 <dbl>,
#> #   PC_10 <dbl>, PC_11 <dbl>, PC_12 <dbl>, PC_13 <dbl>, PC_14 <dbl>,
#> #   PC_15 <dbl>, PC_16 <dbl>, PC_17 <dbl>, PC_18 <dbl>, PC_19 <dbl>,
#> #   tSNE_1 <dbl>, tSNE_2 <dbl>
```
