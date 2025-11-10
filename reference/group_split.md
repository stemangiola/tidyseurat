# Split data frame by groups

**\[experimental\]**

[`group_split()`](https://dplyr.tidyverse.org/reference/group_split.html)
works like [`base::split()`](https://rdrr.io/r/base/split.html) but:

- It uses the grouping structure from
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
  and therefore is subject to the data mask

- It does not name the elements of the list based on the grouping as
  this only works well for a single character grouping variable.
  Instead, use
  [`group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html)
  to access a data frame that defines the groups.

`group_split()` is primarily designed to work with grouped data frames.
You can pass `...` to group and split an ungrouped data frame, but this
is generally not very useful as you want have easy access to the group
metadata.

## Usage

``` r
# S3 method for class 'Seurat'
group_split(.tbl, ..., .keep = TRUE)
```

## Arguments

- .tbl:

  A tbl.

- ...:

  If `.tbl` is an ungrouped data frame, a grouping specification,
  forwarded to
  [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- .keep:

  Should the grouping columns be kept?

## Value

A list of tibbles. Each tibble contains the rows of `.tbl` for the
associated group and all the columns, including the grouping variables.
Note that this returns a
[list_of](https://vctrs.r-lib.org/reference/list_of.html) which is
slightly stricter than a simple list but is useful for representing
lists where every element has the same type.

## Lifecycle

`group_split()` is not stable because you can achieve very similar
results by manipulating the nested column returned from
[`tidyr::nest(.by =)`](https://tidyr.tidyverse.org/reference/nest.html).
That also retains the group keys all within a single data structure.
`group_split()` may be deprecated in the future.

## See also

Other grouping functions:
[`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
[`group_map()`](https://dplyr.tidyverse.org/reference/group_map.html),
[`group_nest()`](https://dplyr.tidyverse.org/reference/group_nest.html),
[`group_trim()`](https://dplyr.tidyverse.org/reference/group_trim.html)

## Examples

``` r
data(pbmc_small)
pbmc_small |> group_split(groups)
#> [[1]]
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
#> 
#> [[2]]
#> # A Seurat-tibble abstraction: 36 × 15
#> # Features=230 | Cells=36 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
#>  1 ATGC… SeuratPro…         70           47 0               A             g2    
#>  2 GAAC… SeuratPro…         87           50 1               B             g2    
#>  3 TGAC… SeuratPro…        127           56 0               A             g2    
#>  4 AGTC… SeuratPro…        173           53 0               A             g2    
#>  5 AGGT… SeuratPro…         62           31 0               A             g2    
#>  6 GGGT… SeuratPro…        101           41 0               A             g2    
#>  7 CATG… SeuratPro…         51           26 0               A             g2    
#>  8 TACG… SeuratPro…         99           45 0               A             g2    
#>  9 GTAA… SeuratPro…         67           33 0               A             g2    
#> 10 TACA… SeuratPro…        109           41 0               A             g2    
#> # ℹ 26 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
#> 
```
