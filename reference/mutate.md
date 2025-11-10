# Create, modify, and delete columns

`mutate()` creates new columns that are functions of existing variables.
It can also modify (if the name is the same as an existing column) and
delete columns (by setting their value to `NULL`).

## Usage

``` r
# S3 method for class 'Seurat'
mutate(.data, ...)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Name-value pairs. The name gives the name of the column in the output.

  The value can be:

  - A vector of length 1, which will be recycled to the correct length.

  - A vector the same length as the current group (or the whole data
    frame if ungrouped).

  - `NULL`, to remove the column.

  - A data frame or tibble, to create multiple columns in the output.

## Value

An object of the same type as `.data`. The output has the following
properties:

- Columns from `.data` will be preserved according to the `.keep`
  argument.

- Existing columns that are modified by `...` will always be returned in
  their original location.

- New columns created through `...` will be placed according to the
  `.before` and `.after` arguments.

- The number of rows is not affected.

- Columns given the value `NULL` will be removed.

- Groups will be recomputed if a grouping variable is mutated.

- Data frame attributes are preserved.

## Useful mutate functions

- [`+`](https://rdrr.io/r/base/Arithmetic.html),
  [`-`](https://rdrr.io/r/base/Arithmetic.html),
  [`log()`](https://rdrr.io/r/base/Log.html), etc., for their usual
  mathematical meanings

- [`lead()`](https://dplyr.tidyverse.org/reference/lead-lag.html),
  [`lag()`](https://dplyr.tidyverse.org/reference/lead-lag.html)

- [`dense_rank()`](https://dplyr.tidyverse.org/reference/row_number.html),
  [`min_rank()`](https://dplyr.tidyverse.org/reference/row_number.html),
  [`percent_rank()`](https://dplyr.tidyverse.org/reference/percent_rank.html),
  [`row_number()`](https://dplyr.tidyverse.org/reference/row_number.html),
  [`cume_dist()`](https://dplyr.tidyverse.org/reference/percent_rank.html),
  [`ntile()`](https://dplyr.tidyverse.org/reference/ntile.html)

- [`cumsum()`](https://rdrr.io/r/base/cumsum.html),
  [`cummean()`](https://dplyr.tidyverse.org/reference/cumall.html),
  [`cummin()`](https://rdrr.io/r/base/cumsum.html),
  [`cummax()`](https://rdrr.io/r/base/cumsum.html),
  [`cumany()`](https://dplyr.tidyverse.org/reference/cumall.html),
  [`cumall()`](https://dplyr.tidyverse.org/reference/cumall.html)

- [`na_if()`](https://dplyr.tidyverse.org/reference/na_if.html),
  [`coalesce()`](https://dplyr.tidyverse.org/reference/coalesce.html)

- [`if_else()`](https://dplyr.tidyverse.org/reference/if_else.html),
  [`recode()`](https://dplyr.tidyverse.org/reference/recode.html),
  [`case_when()`](https://dplyr.tidyverse.org/reference/case_when.html)

## Grouped tibbles

Because mutating expressions are computed within groups, they may yield
different results on grouped tibbles. This will be the case as soon as
an aggregating, lagging, or ranking function is involved. Compare this
ungrouped mutate:

    starwars %>%
      select(name, mass, species) %>%
      mutate(mass_norm = mass / mean(mass, na.rm = TRUE))

With the grouped equivalent:

    starwars %>%
      select(name, mass, species) %>%
      group_by(species) %>%
      mutate(mass_norm = mass / mean(mass, na.rm = TRUE))

The former normalises `mass` by the global average whereas the latter
normalises by the averages within species levels.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages: dplyr (`data.frame`),
plotly ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
tidyseurat (`Seurat`) .

## See also

Other single table verbs: [`arrange()`](arrange.md),
[`rename()`](rename.md), [`slice()`](slice.md),
[`summarise()`](summarise.md)

## Examples

``` r
data(pbmc_small)
pbmc_small |> mutate(nFeature_RNA=1)
#> # A Seurat-tibble abstraction: 80 × 15
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <dbl> <fct>           <fct>         <chr> 
#>  1 ATGC… SeuratPro…         70            1 0               A             g2    
#>  2 CATG… SeuratPro…         85            1 0               A             g1    
#>  3 GAAC… SeuratPro…         87            1 1               B             g2    
#>  4 TGAC… SeuratPro…        127            1 0               A             g2    
#>  5 AGTC… SeuratPro…        173            1 0               A             g2    
#>  6 TCTG… SeuratPro…         70            1 0               A             g1    
#>  7 TGGT… SeuratPro…         64            1 0               A             g1    
#>  8 GCAG… SeuratPro…         72            1 0               A             g1    
#>  9 GATA… SeuratPro…         52            1 0               A             g1    
#> 10 AATG… SeuratPro…        100            1 0               A             g1    
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
