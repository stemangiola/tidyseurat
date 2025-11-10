# Group by one or more variables

Most data operations are done on groups defined by variables.
`group_by()` takes an existing tbl and converts it into a grouped tbl
where operations are performed "by group".
[`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html)
removes grouping.

## Usage

``` r
# S3 method for class 'Seurat'
group_by(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  In `group_by()`, variables or computations to group by. Computations
  are always done on the ungrouped data frame. To perform computations
  on the grouped data, you need to use a separate
  [`mutate()`](mutate.md) step before the `group_by()`. Computations are
  not allowed in
  [`nest_by()`](https://dplyr.tidyverse.org/reference/nest_by.html). In
  [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html),
  variables to remove from the grouping.

- .add:

  When `FALSE`, the default, `group_by()` will override existing groups.
  To add to the existing groups, use `.add = TRUE`.

  This argument was previously called `add`, but that prevented creating
  a new grouping variable called `add`, and conflicts with our naming
  conventions.

- .drop:

  Drop groups formed by factor levels that don't appear in the data? The
  default is `TRUE` except when `.data` has been previously grouped with
  `.drop = FALSE`. See
  [`group_by_drop_default()`](https://dplyr.tidyverse.org/reference/group_by_drop_default.html)
  for details.

## Value

A grouped data frame with class
[`grouped_df`](https://dplyr.tidyverse.org/reference/grouped_df.html),
unless the combination of `...` and `add` yields a empty set of grouping
columns, in which case a tibble will be returned.

## Methods

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `group_by()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- [`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html):
  dplyr (`data.frame`, `grouped_df`, `rowwise_df`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)) .

## Ordering

Currently, `group_by()` internally orders the groups in ascending order.
This results in ordered output from functions that aggregate groups,
such as
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).

When used as grouping columns, character vectors are ordered in the C
locale for performance and reproducibility across R sessions. If the
resulting ordering of your grouped operation matters and is dependent on
the locale, you should follow up the grouped operation with an explicit
call to
[`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html) and
set the `.locale` argument. For example:

    data %>%
      group_by(chr) %>%
      summarise(avg = mean(x)) %>%
      arrange(chr, .locale = "en")

This is often useful as a preliminary step before generating content
intended for humans, such as an HTML table.

### Legacy behavior

Prior to dplyr 1.1.0, character vector grouping columns were ordered in
the system locale. If you need to temporarily revert to this behavior,
you can set the global option `dplyr.legacy_locale` to `TRUE`, but this
should be used sparingly and you should expect this option to be removed
in a future version of dplyr. It is better to update existing code to
explicitly call `arrange(.locale = )` instead. Note that setting
`dplyr.legacy_locale` will also force calls to
[`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html) to use
the system locale.

## See also

Other grouping functions:
[`group_map()`](https://dplyr.tidyverse.org/reference/group_map.html),
[`group_nest()`](https://dplyr.tidyverse.org/reference/group_nest.html),
[`group_split()`](https://dplyr.tidyverse.org/reference/group_split.html),
[`group_trim()`](https://dplyr.tidyverse.org/reference/group_trim.html)

## Examples

``` r
data("pbmc_small")
pbmc_small |>  group_by(groups)
#> tidyseurat says: A data frame is returned for independent data analysis.
#> # A tibble: 80 × 29
#> # Groups:   groups [2]
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
