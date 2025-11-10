# Unnest a list-column of data frames into rows and columns

Unnest expands a list-column containing data frames into rows and
columns.

## Usage

``` r
# S3 method for class 'tidyseurat_nested'
unnest(
  data,
  cols,
  ...,
  keep_empty = FALSE,
  ptype = NULL,
  names_sep = NULL,
  names_repair = "check_unique",
  .drop,
  .id,
  .sep,
  .preserve
)

unnest_seurat(
  data,
  cols,
  ...,
  keep_empty = FALSE,
  ptype = NULL,
  names_sep = NULL,
  names_repair = "check_unique",
  .drop,
  .id,
  .sep,
  .preserve
)
```

## Arguments

- data:

  A data frame.

- cols:

  \<[`tidy-select`](https://tidyr.tidyverse.org/reference/tidyr_tidy_select.html)\>
  List-columns to unnest.

  When selecting multiple columns, values from the same row will be
  recycled to their common size.

- ...:

  **\[deprecated\]**: previously you could write
  `df %>% unnest(x, y, z)`. Convert to `df %>% unnest(c(x, y, z))`. If
  you previously created a new variable in `unnest()` you'll now need to
  do it explicitly with [`mutate()`](mutate.md). Convert
  `df %>% unnest(y = fun(x, y, z))` to
  `df %>% mutate(y = fun(x, y, z)) %>% unnest(y)`.

- keep_empty:

  By default, you get one row of output for each element of the list
  that you are unchopping/unnesting. This means that if there's a size-0
  element (like `NULL` or an empty data frame or vector), then that
  entire row will be dropped from the output. If you want to preserve
  all rows, use `keep_empty = TRUE` to replace size-0 elements with a
  single row of missing values.

- ptype:

  Optionally, a named list of column name-prototype pairs to coerce
  `cols` to, overriding the default that will be guessed from combining
  the individual values. Alternatively, a single empty ptype can be
  supplied, which will be applied to all `cols`.

- names_sep:

  If `NULL`, the default, the outer names will come from the inner
  names. If a string, the outer names will be formed by pasting together
  the outer and the inner column names, separated by `names_sep`.

- names_repair:

  Used to check that output data frame has valid names. Must be one of
  the following options:

  - `"minimal`": no name repair or checks, beyond basic existence,

  - `"unique`": make sure names are unique and not empty,

  - `"check_unique`": (the default), no name repair, but check they are
    unique,

  - `"universal`": make the names unique and syntactic

  - a function: apply custom name repair.

  - [tidyr_legacy](https://tidyr.tidyverse.org/reference/tidyr_legacy.html):
    use the name repair from tidyr 0.8.

  - a formula: a purrr-style anonymous function (see
    [`rlang::as_function()`](https://rlang.r-lib.org/reference/as_function.html))

  See
  [`vctrs::vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html)
  for more details on these terms and the strategies used to enforce
  them.

- .drop, .preserve:

  **\[deprecated\]**: all list-columns are now preserved; If there are
  any that you don't want in the output use [`select()`](select.md) to
  remove them prior to unnesting.

- .id:

  **\[deprecated\]**: convert `df %>% unnest(x, .id = "id")` to
  `df %>% mutate(id = names(x)) %>% unnest(x))`.

- .sep:

  **\[deprecated\]**: use `names_sep` instead.

## Value

\`tidyseurat\`

## New syntax

tidyr 1.0.0 introduced a new syntax for [`nest()`](nest.md) and
`unnest()` that's designed to be more similar to other functions.
Converting to the new syntax should be straightforward (guided by the
message you'll receive) but if you just need to run an old analysis, you
can easily revert to the previous behaviour using
[`nest_legacy()`](https://tidyr.tidyverse.org/reference/nest_legacy.html)
and
[`unnest_legacy()`](https://tidyr.tidyverse.org/reference/nest_legacy.html)
as follows:

    library(tidyr)
    nest <- nest_legacy
    unnest <- unnest_legacy

## See also

Other rectangling:
[`hoist()`](https://tidyr.tidyverse.org/reference/hoist.html),
[`unnest_longer()`](https://tidyr.tidyverse.org/reference/unnest_longer.html),
[`unnest_wider()`](https://tidyr.tidyverse.org/reference/unnest_wider.html)

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
