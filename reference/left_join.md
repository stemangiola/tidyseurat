# Mutating joins

Mutating joins add columns from `y` to `x`, matching observations based
on the keys. There are four mutating joins: the inner join, and the
three outer joins.

### Inner join

An [`inner_join()`](inner_join.md) only keeps observations from `x` that
have a matching key in `y`.

The most important property of an inner join is that unmatched rows in
either input are not included in the result. This means that generally
inner joins are not appropriate in most analyses, because it is too easy
to lose observations.

### Outer joins

The three outer joins keep observations that appear in at least one of
the data frames:

- A `left_join()` keeps all observations in `x`.

- A [`right_join()`](right_join.md) keeps all observations in `y`.

- A [`full_join()`](full_join.md) keeps all observations in `x` and `y`.

## Usage

``` r
# S3 method for class 'Seurat'
left_join(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...)
```

## Arguments

- x, y:

  A pair of data frames, data frame extensions (e.g. a tibble), or lazy
  data frames (e.g. from dbplyr or dtplyr). See *Methods*, below, for
  more details.

- by:

  A join specification created with
  [`join_by()`](https://dplyr.tidyverse.org/reference/join_by.html), or
  a character vector of variables to join by.

  If `NULL`, the default, `*_join()` will perform a natural join, using
  all variables in common across `x` and `y`. A message lists the
  variables so that you can check they're correct; suppress the message
  by supplying `by` explicitly.

  To join on different variables between `x` and `y`, use a
  [`join_by()`](https://dplyr.tidyverse.org/reference/join_by.html)
  specification. For example, `join_by(a == b)` will match `x$a` to
  `y$b`.

  To join by multiple variables, use a
  [`join_by()`](https://dplyr.tidyverse.org/reference/join_by.html)
  specification with multiple expressions. For example,
  `join_by(a == b, c == d)` will match `x$a` to `y$b` and `x$c` to
  `y$d`. If the column names are the same between `x` and `y`, you can
  shorten this by listing only the variable names, like `join_by(a, c)`.

  [`join_by()`](https://dplyr.tidyverse.org/reference/join_by.html) can
  also be used to perform inequality, rolling, and overlap joins. See
  the documentation at
  [?join_by](https://dplyr.tidyverse.org/reference/join_by.html) for
  details on these types of joins.

  For simple equality joins, you can alternatively specify a character
  vector of variable names to join by. For example, `by = c("a", "b")`
  joins `x$a` to `y$a` and `x$b` to `y$b`. If variable names differ
  between `x` and `y`, use a named character vector like
  `by = c("x_a" = "y_a", "x_b" = "y_b")`.

  To perform a cross-join, generating all combinations of `x` and `y`,
  see
  [`cross_join()`](https://dplyr.tidyverse.org/reference/cross_join.html).

- copy:

  If `x` and `y` are not from the same data source, and `copy` is
  `TRUE`, then `y` will be copied into the same src as `x`. This allows
  you to join tables across srcs, but it is a potentially expensive
  operation so you must opt into it.

- suffix:

  If there are non-joined duplicate variables in `x` and `y`, these
  suffixes will be added to the output to disambiguate them. Should be a
  character vector of length 2.

- ...:

  Other parameters passed onto methods.

## Value

An object of the same type as `x` (including the same groups). The order
of the rows and columns of `x` is preserved as much as possible. The
output has the following properties:

- The rows are affect by the join type.

  - [`inner_join()`](inner_join.md) returns matched `x` rows.

  - `left_join()` returns all `x` rows.

  - [`right_join()`](right_join.md) returns matched of `x` rows,
    followed by unmatched `y` rows.

  - [`full_join()`](full_join.md) returns all `x` rows, followed by
    unmatched `y` rows.

- Output columns include all columns from `x` and all non-key columns
  from `y`. If `keep = TRUE`, the key columns from `y` are included as
  well.

- If non-key columns in `x` and `y` have the same name, `suffix`es are
  added to disambiguate. If `keep = TRUE` and key columns in `x` and `y`
  have the same name, `suffix`es are added to disambiguate these as
  well.

- If `keep = FALSE`, output columns included in `by` are coerced to
  their common type between `x` and `y`.

## Many-to-many relationships

By default, dplyr guards against many-to-many relationships in equality
joins by throwing a warning. These occur when both of the following are
true:

- A row in `x` matches multiple rows in `y`.

- A row in `y` matches multiple rows in `x`.

This is typically surprising, as most joins involve a relationship of
one-to-one, one-to-many, or many-to-one, and is often the result of an
improperly specified join. Many-to-many relationships are particularly
problematic because they can result in a Cartesian explosion of the
number of rows returned from the join.

If a many-to-many relationship is expected, silence this warning by
explicitly setting `relationship = "many-to-many"`.

In production code, it is best to preemptively set `relationship` to
whatever relationship you expect to exist between the keys of `x` and
`y`, as this forces an error to occur immediately if the data doesn't
align with your expectations.

Inequality joins typically result in many-to-many relationships by
nature, so they don't warn on them by default, but you should still take
extra care when specifying an inequality join, because they also have
the capability to return a large number of rows.

Rolling joins don't warn on many-to-many relationships either, but many
rolling joins follow a many-to-one relationship, so it is often useful
to set `relationship = "many-to-one"` to enforce this.

Note that in SQL, most database providers won't let you specify a
many-to-many relationship between two tables, instead requiring that you
create a third *junction table* that results in two one-to-many
relationships instead.

## Methods

These functions are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- [`inner_join()`](inner_join.md): dplyr
  ([`data.frame`](https://dplyr.tidyverse.org/reference/mutate-joins.html)),
  tidyseurat (`Seurat`) .

- `left_join()`: dplyr
  ([`data.frame`](https://dplyr.tidyverse.org/reference/mutate-joins.html)),
  tidyseurat (`Seurat`) .

- [`right_join()`](right_join.md): dplyr
  ([`data.frame`](https://dplyr.tidyverse.org/reference/mutate-joins.html)),
  tidyseurat (`Seurat`) .

- [`full_join()`](full_join.md): dplyr
  ([`data.frame`](https://dplyr.tidyverse.org/reference/mutate-joins.html)),
  tidyseurat (`Seurat`) .

## See also

Other joins:
[`cross_join()`](https://dplyr.tidyverse.org/reference/cross_join.html),
[`filter-joins`](https://dplyr.tidyverse.org/reference/filter-joins.html),
[`nest_join()`](https://dplyr.tidyverse.org/reference/nest_join.html)

## Examples

``` r
data(pbmc_small)
tt <- pbmc_small
tt |> left_join(tt |>  
  distinct(groups) |> 
  mutate(new_column=1:2))
#> tidyseurat says: A data frame is returned for independent data analysis.
#> Joining with `by = join_by(groups)`
#> # A Seurat-tibble abstraction: 80 × 16
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
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
#> # ℹ 9 more variables: RNA_snn_res.1 <fct>, new_column <int>, PC_1 <dbl>,
#> #   PC_2 <dbl>, PC_3 <dbl>, PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
