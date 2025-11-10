# Keep rows that match a condition

The `filter()` function is used to subset a data frame, retaining all
rows that satisfy your conditions. To be retained, the row must produce
a value of `TRUE` for all conditions. Note that when a condition
evaluates to `NA` the row will be dropped, unlike base subsetting with
`[`.

## Usage

``` r
# S3 method for class 'Seurat'
filter(.data, ..., .preserve = FALSE)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Expressions that return a logical value, and are defined in terms of
  the variables in `.data`. If multiple expressions are included, they
  are combined with the `&` operator. Only rows for which all conditions
  evaluate to `TRUE` are kept.

- .preserve:

  Relevant when the `.data` input is grouped. If `.preserve = FALSE`
  (the default), the grouping structure is recalculated based on the
  resulting data, otherwise the grouping is kept as is.

## Value

An object of the same type as `.data`. The output has the following
properties:

- Rows are a subset of the input, but appear in the same order.

- Columns are not modified.

- The number of groups may be reduced (if `.preserve` is not `TRUE`).

- Data frame attributes are preserved.

## Details

The `filter()` function is used to subset the rows of `.data`, applying
the expressions in `...` to the column values to determine which rows
should be retained. It can be applied to both grouped and ungrouped data
(see [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html)).
However, dplyr is not yet smart enough to optimise the filtering
operation on grouped datasets that do not need grouped calculations. For
this reason, filtering is often considerably faster on ungrouped data.

## Useful filter functions

There are many functions and operators that are useful when constructing
the expressions used to filter the data:

- [`==`](https://rdrr.io/r/base/Comparison.html), `>`, `>=` etc

- `&`, [`|`](https://rdrr.io/r/base/Logic.html),
  [`!`](https://rdrr.io/r/base/Logic.html),
  [`xor()`](https://rdrr.io/r/base/Logic.html)

- [`is.na()`](https://rdrr.io/r/base/NA.html)

- [`between()`](https://dplyr.tidyverse.org/reference/between.html),
  [`near()`](https://dplyr.tidyverse.org/reference/near.html)

## Grouped tibbles

Because filtering expressions are computed within groups, they may yield
different results on grouped tibbles. This will be the case as soon as
an aggregating, lagging, or ranking function is involved. Compare this
ungrouped filtering:

    starwars %>% filter(mass > mean(mass, na.rm = TRUE))

With the grouped equivalent:

    starwars %>% group_by(gender) %>% filter(mass > mean(mass, na.rm = TRUE))

In the ungrouped version, `filter()` compares the value of `mass` in
each row to the global average (taken over the whole data set), keeping
only the rows with `mass` greater than this global average. In contrast,
the grouped version calculates the average mass separately for each
`gender` group, and keeps rows with `mass` greater than the relevant
within-gender average.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`, `ts`), plotly
([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## See also

Other single table verbs:
[`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html),
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
[`reframe()`](https://dplyr.tidyverse.org/reference/reframe.html),
[`rename()`](https://dplyr.tidyverse.org/reference/rename.html),
[`select()`](https://dplyr.tidyverse.org/reference/select.html),
[`slice()`](https://dplyr.tidyverse.org/reference/slice.html),
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)

## Examples

``` r
data("pbmc_small")
pbmc_small |>  filter(groups == "g1")
#> # A Seurat-tibble abstraction: 44 × 15
#> # Features=230 | Cells=44 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 CATG… SeuratPro…         85           52 0               A             g1    
#>  2 TCTG… SeuratPro…         70           48 0               A             g1    
#>  3 TGGT… SeuratPro…         64           36 0               A             g1    
#>  4 GCAG… SeuratPro…         72           45 0               A             g1    
#>  5 GATA… SeuratPro…         52           36 0               A             g1    
#>  6 AATG… SeuratPro…        100           41 0               A             g1    
#>  7 AGAG… SeuratPro…        191           61 0               A             g1    
#>  8 CTAA… SeuratPro…        168           44 0               A             g1    
#>  9 TTGG… SeuratPro…        135           45 0               A             g1    
#> 10 CATC… SeuratPro…         79           43 0               A             g1    
#> # ℹ 34 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# Learn more in ?dplyr_eval
```
