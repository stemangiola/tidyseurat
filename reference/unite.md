# Unite multiple columns into one by pasting strings together

Convenience function to paste together multiple columns into one.

## Usage

``` r
# S3 method for class 'Seurat'
unite(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)
```

## Arguments

- data:

  A data frame.

- col:

  The name of the new column, as a string or symbol.

  This argument is passed by expression and supports
  [quasiquotation](https://rlang.r-lib.org/reference/topic-inject.html)
  (you can unquote strings and symbols). The name is captured from the
  expression with
  [`rlang::ensym()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
  (note that this kind of interface where symbols do not represent
  actual objects is now discouraged in the tidyverse; we support it here
  for backward compatibility).

- ...:

  \<[`tidy-select`](https://tidyr.tidyverse.org/reference/tidyr_tidy_select.html)\>
  Columns to unite

- sep:

  Separator to use between values.

- remove:

  If `TRUE`, remove input columns from output data frame.

- na.rm:

  If `TRUE`, missing values will be removed prior to uniting each value.

## Value

\`tidyseurat\`

## See also

[`separate()`](https://tidyr.tidyverse.org/reference/separate.html), the
complement.

## Examples

``` r
data(pbmc_small)
pbmc_small |> unite(
  col="new_col", 
  c("orig.ident", "groups"))
#> # A Seurat-tibble abstraction: 80 × 14
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell          new_col  nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
#>    <chr>          <chr>         <dbl>        <int> <fct>           <fct>        
#>  1 ATGCCAGAACGACT SeuratP…         70           47 0               A            
#>  2 CATGGCCTGTGCAT SeuratP…         85           52 0               A            
#>  3 GAACCTGATGAACC SeuratP…         87           50 1               B            
#>  4 TGACTGGATTCTCA SeuratP…        127           56 0               A            
#>  5 AGTCAGACTGCACA SeuratP…        173           53 0               A            
#>  6 TCTGATACACGTGT SeuratP…         70           48 0               A            
#>  7 TGGTATCTAAACAG SeuratP…         64           36 0               A            
#>  8 GCAGCTCTGTTTCT SeuratP…         72           45 0               A            
#>  9 GATATAACACGCAT SeuratP…         52           36 0               A            
#> 10 AATGTTGACAGTCA SeuratP…        100           41 0               A            
#> # ℹ 70 more rows
#> # ℹ 8 more variables: RNA_snn_res.1 <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>,
#> #   PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
    
```
