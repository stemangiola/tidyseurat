# Group input by rows

`rowwise()` allows you to compute on a data frame a row-at-a-time. This
is most useful when a vectorised function doesn't exist.

Most dplyr verbs preserve row-wise grouping. The exception is
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html),
which return a
[grouped_df](https://dplyr.tidyverse.org/reference/grouped_df.html). You
can explicitly ungroup with
[`ungroup()`](https://dplyr.tidyverse.org/reference/group_by.html) or
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
or convert to a
[grouped_df](https://dplyr.tidyverse.org/reference/grouped_df.html) with
[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

## Usage

``` r
# S3 method for class 'Seurat'
rowwise(data, ...)
```

## Arguments

- data:

  Input data frame.

- ...:

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Variables to be preserved when calling
  [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).
  This is typically a set of variables whose combination uniquely
  identify each row.

  **NB**: unlike [`group_by()`](group_by.md) you can not create new
  variables here but instead you can select multiple variables with
  (e.g.)
  [`everything()`](https://tidyselect.r-lib.org/reference/everything.html).

## Value

A row-wise data frame with class `rowwise_df`. Note that a `rowwise_df`
is implicitly grouped by row, but is not a `grouped_df`.

## List-columns

Because a rowwise has exactly one row per group it offers a small
convenience for working with list-columns. Normally,
[`summarise()`](summarise.md) and [`mutate()`](mutate.md) extract a
groups worth of data with `[`. But when you index a list in this way,
you get back another list. When you're working with a `rowwise` tibble,
then dplyr will use `[[` instead of `[` to make your life a little
easier.

## See also

[`nest_by()`](https://dplyr.tidyverse.org/reference/nest_by.html) for a
convenient way of creating rowwise data frames with nested data.

## Examples

``` r
# TODO
```
