# Keep or drop rows that match a condition

These functions are used to subset a data frame, applying the
expressions in `...` to determine which rows should be kept (for
`filter()`) or dropped ( for
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)).

Multiple conditions can be supplied separated by a comma. These will be
combined with the `&` operator. To combine comma separated conditions
using `|` instead, wrap them in
[`when_any()`](https://dplyr.tidyverse.org/reference/when-any-all.html).

Both `filter()` and
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
treat `NA` like `FALSE`. This subtle behavior can impact how you write
your conditions when missing values are involved. See the section on
`Missing values` for important details and examples.

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
  Expressions that return a logical vector, defined in terms of the
  variables in `.data`. If multiple expressions are included, they are
  combined with the `&` operator. To combine expressions using `|`
  instead, wrap them in
  [`when_any()`](https://dplyr.tidyverse.org/reference/when-any-all.html).
  Only rows for which all expressions evaluate to `TRUE` are kept (for
  `filter()`) or dropped (for
  [`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)).

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

## Missing values

Both `filter()` and
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
treat `NA` like `FALSE`. This results in the following behavior:

- `filter()` *drops* both `NA` and `FALSE`.

- [`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
  *keeps* both `NA` and `FALSE`.

This means that
`filter(data, <conditions>) + filter_out(data, <conditions>)` captures
every row within `data` exactly once.

The `NA` handling of these functions has been designed to match your
*intent*. When your intent is to keep rows, use `filter()`. When your
intent is to drop rows, use
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html).

For example, if your goal with this `cars` data is to "drop rows where
the `class` is suv", then you might write this in one of two ways:

    cars <- tibble(class = c("suv", NA, "coupe"))
    cars
    #> # A tibble: 3 x 1
    #>   class
    #>   <chr>
    #> 1 suv
    #> 2 <NA>
    #> 3 coupe

    cars |> filter(class != "suv")
    #> # A tibble: 1 x 1
    #>   class
    #>   <chr>
    #> 1 coupe

    cars |> filter_out(class == "suv")
    #> # A tibble: 2 x 1
    #>   class
    #>   <chr>
    #> 1 <NA>
    #> 2 coupe

Note how `filter()` drops the `NA` rows even though our goal was only to
drop `"suv"` rows, but
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
matches our intuition.

To generate the correct result with `filter()`, you'd need to use:

    cars |> filter(class != "suv" | is.na(class))
    #> # A tibble: 2 x 1
    #>   class
    #>   <chr>
    #> 1 <NA>
    #> 2 coupe

This quickly gets unwieldy when multiple conditions are involved.

In general, if you find yourself:

- Using "negative" operators like `!=` or `!`

- Adding in `NA` handling like `| is.na(col)` or `& !is.na(col)`

then you should consider if swapping to the other filtering variant
would make your conditions simpler.

### Comparison to base subsetting

Base subsetting with `[` doesn't treat `NA` like `TRUE` or `FALSE`.
Instead, it generates a fully missing row, which is different from how
both `filter()` and
[`filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
work.

    cars <- tibble(class = c("suv", NA, "coupe"), mpg = c(10, 12, 14))
    cars
    #> # A tibble: 3 x 2
    #>   class   mpg
    #>   <chr> <dbl>
    #> 1 suv      10
    #> 2 <NA>     12
    #> 3 coupe    14

    cars[cars$class == "suv",]
    #> # A tibble: 2 x 2
    #>   class   mpg
    #>   <chr> <dbl>
    #> 1 suv      10
    #> 2 <NA>     NA

    cars |> filter(class == "suv")
    #> # A tibble: 1 x 2
    #>   class   mpg
    #>   <chr> <dbl>
    #> 1 suv      10

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

- [`when_any()`](https://dplyr.tidyverse.org/reference/when-any-all.html),
  [`when_all()`](https://dplyr.tidyverse.org/reference/when-any-all.html)

## Grouped tibbles

Because filtering expressions are computed within groups, they may yield
different results on grouped tibbles. This will be the case as soon as
an aggregating, lagging, or ranking function is involved. Compare this
ungrouped filtering:

    starwars |> filter(mass > mean(mass, na.rm = TRUE))

With the grouped equivalent:

    starwars |> filter(mass > mean(mass, na.rm = TRUE), .by = gender)

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
