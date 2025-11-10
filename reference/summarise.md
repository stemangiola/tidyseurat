# Summarise each group down to one row

`summarise()` creates a new data frame. It returns one row for each
combination of grouping variables; if there are no grouping variables,
the output will have a single row summarising all observations in the
input. It will contain one column for each grouping variable and one
column for each of the summary statistics that you have specified.

`summarise()` and `summarize()` are synonyms.

## Usage

``` r
# S3 method for class 'Seurat'
summarise(.data, ...)

# S3 method for class 'Seurat'
summarize(.data, ...)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Name-value pairs of summary functions. The name will be the name of
  the variable in the result.

  The value can be:

  - A vector of length 1, e.g. `min(x)`,
    [`n()`](https://dplyr.tidyverse.org/reference/context.html), or
    `sum(is.na(y))`.

  - A data frame, to add multiple columns from a single expression.

  **\[deprecated\]** Returning values with size 0 or \>1 was deprecated
  as of 1.1.0. Please use
  [`reframe()`](https://dplyr.tidyverse.org/reference/reframe.html) for
  this instead.

## Value

An object *usually* of the same type as `.data`.

- The rows come from the underlying
  [`group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html).

- The columns are a combination of the grouping keys and the summary
  expressions that you provide.

- The grouping structure is controlled by the `.groups=` argument, the
  output may be another
  [grouped_df](https://dplyr.tidyverse.org/reference/grouped_df.html), a
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) or a
  [rowwise](https://dplyr.tidyverse.org/reference/rowwise.html) data
  frame.

- Data frame attributes are **not** preserved, because `summarise()`
  fundamentally creates a new data frame.

## Useful functions

- Center: [`mean()`](https://rdrr.io/r/base/mean.html),
  [`median()`](https://rdrr.io/r/stats/median.html)

- Spread: [`sd()`](https://rdrr.io/r/stats/sd.html),
  [`IQR()`](https://rdrr.io/r/stats/IQR.html),
  [`mad()`](https://rdrr.io/r/stats/mad.html)

- Range: [`min()`](https://rdrr.io/r/base/Extremes.html),
  [`max()`](https://rdrr.io/r/base/Extremes.html),

- Position: [`first()`](https://dplyr.tidyverse.org/reference/nth.html),
  [`last()`](https://dplyr.tidyverse.org/reference/nth.html),
  [`nth()`](https://dplyr.tidyverse.org/reference/nth.html),

- Count: [`n()`](https://dplyr.tidyverse.org/reference/context.html),
  [`n_distinct()`](https://dplyr.tidyverse.org/reference/n_distinct.html)

- Logical: [`any()`](https://rdrr.io/r/base/any.html),
  [`all()`](https://rdrr.io/r/base/all.html)

## Backend variations

The data frame backend supports creating a variable and using it in the
same summary. This means that previously created summary variables can
be further transformed or combined within the summary, as in
[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
However, it also means that summary variables with the same names as
previous variables overwrite them, making those variables unavailable to
later summary variables.

This behaviour may not be supported in other backends. To avoid
unexpected results, consider using new names for your summary variables,
especially when creating multiple summaries.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`, `grouped_df`, `rowwise_df`), plotly
([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## See also

Other single table verbs: [`arrange()`](arrange.md),
[`mutate()`](mutate.md), [`rename()`](rename.md), [`slice()`](slice.md)

## Examples

``` r
data(pbmc_small)
pbmc_small |> summarise(mean(nCount_RNA))
#> tidyseurat says: A data frame is returned for independent data analysis.
#> # A tibble: 1 × 1
#>   `mean(nCount_RNA)`
#>                <dbl>
#> 1               245.
```
