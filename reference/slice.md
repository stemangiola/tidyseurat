# Subset rows using their positions

`slice()` lets you index rows by their (integer) locations. It allows
you to select, remove, and duplicate rows. It is accompanied by a number
of helpers for common use cases:

- `slice_head()` and `slice_tail()` select the first or last rows.

- `slice_sample()` randomly selects rows.

- `slice_min()` and `slice_max()` select rows with the smallest or
  largest values of a variable.

If `.data` is a
[grouped_df](https://dplyr.tidyverse.org/reference/grouped_df.html), the
operation will be performed on each group, so that (e.g.)
`slice_head(df, n = 5)` will select the first five rows in each group.

## Usage

``` r
# S3 method for class 'Seurat'
slice(.data, ..., .by = NULL, .preserve = FALSE)

# S3 method for class 'Seurat'
slice_sample(
  .data,
  ...,
  n = NULL,
  prop = NULL,
  by = NULL,
  weight_by = NULL,
  replace = FALSE
)

# S3 method for class 'Seurat'
slice_head(.data, ..., n, prop, by = NULL)

# S3 method for class 'Seurat'
slice_tail(.data, ..., n, prop, by = NULL)

# S3 method for class 'Seurat'
slice_min(
  .data,
  order_by,
  ...,
  n,
  prop,
  by = NULL,
  with_ties = TRUE,
  na_rm = FALSE
)

# S3 method for class 'Seurat'
slice_max(
  .data,
  order_by,
  ...,
  n,
  prop,
  by = NULL,
  with_ties = TRUE,
  na_rm = FALSE
)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- ...:

  For `slice()`:
  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Integer row values.

  Provide either positive values to keep, or negative values to drop.
  The values provided must be either all positive or all negative.
  Indices beyond the number of rows in the input are silently ignored.

  For `slice_*()`, these arguments are passed on to methods.

- .by, by:

  **\[experimental\]**

  \<[`tidy-select`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  Optionally, a selection of columns to group by for just this
  operation, functioning as an alternative to
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
  For details and examples, see
  [?dplyr_by](https://dplyr.tidyverse.org/reference/dplyr_by.html).

- .preserve:

  Relevant when the `.data` input is grouped. If `.preserve = FALSE`
  (the default), the grouping structure is recalculated based on the
  resulting data, otherwise the grouping is kept as is.

- n, prop:

  Provide either `n`, the number of rows, or `prop`, the proportion of
  rows to select. If neither are supplied, `n = 1` will be used. If `n`
  is greater than the number of rows in the group (or `prop > 1`), the
  result will be silently truncated to the group size. `prop` will be
  rounded towards zero to generate an integer number of rows.

  A negative value of `n` or `prop` will be subtracted from the group
  size. For example, `n = -2` with a group of 5 rows will select 5 - 2 =
  3 rows; `prop = -0.25` with 8 rows will select 8 \* (1 - 0.25) = 6
  rows.

- weight_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Sampling weights. This must evaluate to a vector of non-negative
  numbers the same length as the input. Weights are automatically
  standardised to sum to 1.

- replace:

  Should sampling be performed with (`TRUE`) or without (`FALSE`, the
  default) replacement.

- order_by:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Variable or function of variables to order by. To order by multiple
  variables, wrap them in a data frame or tibble.

- with_ties:

  Should ties be kept together? The default, `TRUE`, may return more
  rows than you request. Use `FALSE` to ignore ties, and return the
  first `n` rows.

- na_rm:

  Should missing values in `order_by` be removed from the result? If
  `FALSE`, `NA` values are sorted to the end (like in
  [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)), so
  they will only be included if there are insufficient non-missing
  values to reach `n`/`prop`.

## Value

An object of the same type as `.data`. The output has the following
properties:

- Each row may appear 0, 1, or many times in the output.

- Columns are not modified.

- Groups are not modified.

- Data frame attributes are preserved.

## Details

Slice does not work with relational databases because they have no
intrinsic notion of row order. If you want to perform the equivalent
operation, use
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html) and
[`row_number()`](https://dplyr.tidyverse.org/reference/row_number.html).

## Methods

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

These function are **generic**s, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

Methods available in currently loaded packages:

- `slice()`: dplyr (`data.frame`), plotly
  ([`plotly`](https://rdrr.io/pkg/plotly/man/plotly_data.html)),
  tidyseurat (`Seurat`) .

- `slice_head()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_tail()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_min()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_max()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

- `slice_sample()`: dplyr (`data.frame`), tidyseurat (`Seurat`) .

## See also

Other single table verbs: [`arrange()`](arrange.md),
[`mutate()`](mutate.md), [`rename()`](rename.md),
[`summarise()`](summarise.md)

## Examples

``` r
data(pbmc_small)
pbmc_small |> slice(1)
#> # A Seurat-tibble abstraction: 1 × 15
#> # Features=230 | Cells=1 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ATGCC… SeuratPro…         70           47 0               A             g2    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# Slice group-wise using .by
pbmc_small |> slice(1:2, .by=groups)
#> # A Seurat-tibble abstraction: 4 × 15
#> # Features=230 | Cells=4 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ATGCC… SeuratPro…         70           47 0               A             g2    
#> 2 CATGG… SeuratPro…         85           52 0               A             g1    
#> 3 GAACC… SeuratPro…         87           50 1               B             g2    
#> 4 TCTGA… SeuratPro…         70           48 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>


# slice_sample() allows you to random select with or without replacement
pbmc_small |> slice_sample(n=5)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ATGCC… SeuratPro…         70           47 0               A             g2    
#> 2 AGATA… SeuratPro…        187           61 0               A             g2    
#> 3 GGCAT… SeuratPro…        126           53 0               A             g1    
#> 4 TACAA… SeuratPro…        108           44 0               A             g2    
#> 5 GATAG… SeuratPro…        328           72 1               B             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# if using replacement, and duplicate cells are returned, a tibble will be
# returned because duplicate cells cannot exist in Seurat objects
pbmc_small |> slice_sample(n=1, replace=TRUE) # returns Seurat
#> # A Seurat-tibble abstraction: 1 × 15
#> # Features=230 | Cells=1 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 GATAG… SeuratPro…        328           72 1               B             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
pbmc_small |> slice_sample(n=100, replace=TRUE) # returns tibble
#> tidyseurat says: When sampling with replacement a data frame is returned for independent data analysis.
#> # A tibble: 100 × 29
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 ATGC… SeuratPro…         70           47 0               A             g2    
#>  2 ATGC… SeuratPro…         70           47 0               A             g2    
#>  3 ATGC… SeuratPro…         70           47 0               A             g2    
#>  4 CATG… SeuratPro…         85           52 0               A             g1    
#>  5 TCTG… SeuratPro…         70           48 0               A             g1    
#>  6 TGGT… SeuratPro…         64           36 0               A             g1    
#>  7 AATG… SeuratPro…        100           41 0               A             g1    
#>  8 GGGT… SeuratPro…        101           41 0               A             g2    
#>  9 GGGT… SeuratPro…        101           41 0               A             g2    
#> 10 CATG… SeuratPro…         51           26 0               A             g2    
#> # ℹ 90 more rows
#> # ℹ 22 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, PC_6 <dbl>, PC_7 <dbl>, PC_8 <dbl>, PC_9 <dbl>,
#> #   PC_10 <dbl>, PC_11 <dbl>, PC_12 <dbl>, PC_13 <dbl>, PC_14 <dbl>,
#> #   PC_15 <dbl>, PC_16 <dbl>, PC_17 <dbl>, PC_18 <dbl>, PC_19 <dbl>,
#> #   tSNE_1 <dbl>, tSNE_2 <dbl>

# weight by a variable
pbmc_small |> slice_sample(n=5, weight_by=nCount_RNA)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ACTCG… SeuratPro…        231           49 1               B             g2    
#> 2 GGCAT… SeuratPro…        126           53 0               A             g1    
#> 3 CTGCC… SeuratPro…        146           47 0               A             g1    
#> 4 AAGCG… SeuratPro…        443           77 1               B             g1    
#> 5 ACCAG… SeuratPro…        417           75 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# sample by group
pbmc_small |> slice_sample(n=5, by=groups)
#> # A Seurat-tibble abstraction: 10 × 15
#> # Features=230 | Cells=10 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 ATGC… SeuratPro…         70           47 0               A             g2    
#>  2 AGTC… SeuratPro…        173           53 0               A             g2    
#>  3 GCGC… SeuratPro…        213           48 1               B             g2    
#>  4 CATC… SeuratPro…        353           80 1               B             g1    
#>  5 TACT… SeuratPro…        156           48 0               A             g1    
#>  6 GGCA… SeuratPro…        126           53 0               A             g1    
#>  7 TTGC… SeuratPro…        104           40 0               A             g2    
#>  8 ATAC… SeuratPro…        612           69 1               B             g1    
#>  9 GTCA… SeuratPro…        210           33 0               A             g2    
#> 10 TTAC… SeuratPro…        228           39 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# sample using proportions
pbmc_small |> slice_sample(prop=0.10)
#> # A Seurat-tibble abstraction: 8 × 15
#> # Features=230 | Cells=8 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ATGCC… SeuratPro…         70           47 0               A             g2    
#> 2 CATGA… SeuratPro…         51           26 0               A             g2    
#> 3 AGATA… SeuratPro…        187           61 0               A             g2    
#> 4 GGCAT… SeuratPro…        126           53 0               A             g1    
#> 5 TACAA… SeuratPro…        108           44 0               A             g2    
#> 6 CGTAG… SeuratPro…        371           75 1               B             g1    
#> 7 GCACT… SeuratPro…        292           71 1               B             g2    
#> 8 GATAG… SeuratPro…        328           72 1               B             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>


# First rows based on existing order
pbmc_small |> slice_head(n=5)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 ATGCC… SeuratPro…         70           47 0               A             g2    
#> 2 CATGG… SeuratPro…         85           52 0               A             g1    
#> 3 GAACC… SeuratPro…         87           50 1               B             g2    
#> 4 TGACT… SeuratPro…        127           56 0               A             g2    
#> 5 AGTCA… SeuratPro…        173           53 0               A             g2    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>


# Last rows based on existing order
pbmc_small |> slice_tail(n=5)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 GAGTT… SeuratPro…        527           47 0               A             g1    
#> 2 GACGC… SeuratPro…        202           30 0               A             g2    
#> 3 AGTCT… SeuratPro…        157           29 0               A             g1    
#> 4 GGAAC… SeuratPro…        150           30 0               A             g2    
#> 5 CTTGA… SeuratPro…        233           76 1               B             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>


# Rows with minimum and maximum values of a metadata variable
pbmc_small |> slice_min(nFeature_RNA, n=5)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 CATGA… SeuratPro…         51           26 0               A             g2    
#> 2 GGCAT… SeuratPro…        172           29 0               A             g1    
#> 3 GACGC… SeuratPro…        202           30 0               A             g2    
#> 4 AGTCT… SeuratPro…        157           29 0               A             g1    
#> 5 GGAAC… SeuratPro…        150           30 0               A             g2    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# slice_min() and slice_max() may return more rows than requested
# in the presence of ties.
pbmc_small |>  slice_min(nFeature_RNA, n=2)
#> # A Seurat-tibble abstraction: 3 × 15
#> # Features=230 | Cells=3 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 CATGA… SeuratPro…         51           26 0               A             g2    
#> 2 GGCAT… SeuratPro…        172           29 0               A             g1    
#> 3 AGTCT… SeuratPro…        157           29 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# Use with_ties=FALSE to return exactly n matches
pbmc_small |> slice_min(nFeature_RNA, n=2, with_ties=FALSE)
#> # A Seurat-tibble abstraction: 2 × 15
#> # Features=230 | Cells=2 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 CATGA… SeuratPro…         51           26 0               A             g2    
#> 2 GGCAT… SeuratPro…        172           29 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# Or use additional variables to break the tie:
pbmc_small |> slice_min(tibble::tibble(nFeature_RNA, nCount_RNA), n=2)
#> # A Seurat-tibble abstraction: 2 × 15
#> # Features=230 | Cells=2 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 CATGA… SeuratPro…         51           26 0               A             g2    
#> 2 AGTCT… SeuratPro…        157           29 0               A             g1    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

# Use by for group-wise operations
pbmc_small |> slice_min(nFeature_RNA, n=5, by=groups)
#> # A Seurat-tibble abstraction: 10 × 15
#> # Features=230 | Cells=10 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 TGGT… SeuratPro…         64           36 0               A             g1    
#>  2 GATA… SeuratPro…         52           36 0               A             g1    
#>  3 AGGT… SeuratPro…         62           31 0               A             g2    
#>  4 CATG… SeuratPro…         51           26 0               A             g2    
#>  5 CTTC… SeuratPro…         41           32 0               A             g2    
#>  6 GGCA… SeuratPro…        172           29 0               A             g1    
#>  7 TTAC… SeuratPro…        228           39 0               A             g1    
#>  8 GACG… SeuratPro…        202           30 0               A             g2    
#>  9 AGTC… SeuratPro…        157           29 0               A             g1    
#> 10 GGAA… SeuratPro…        150           30 0               A             g2    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>


# Rows with minimum and maximum values of a metadata variable
pbmc_small |> slice_max(nFeature_RNA, n=5)
#> # A Seurat-tibble abstraction: 5 × 15
#> # Features=230 | Cells=5 | Active assay=RNA | Assays=RNA
#>   .cell  orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>   <chr>  <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#> 1 TTGAG… SeuratPro…        787           88 0               A             g1    
#> 2 TTTAG… SeuratPro…        462           86 1               B             g1    
#> 3 GACAT… SeuratPro…        872           96 1               B             g1    
#> 4 ACGTG… SeuratPro…        709           94 1               B             g2    
#> 5 ATTGT… SeuratPro…        745           84 1               B             g2    
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
