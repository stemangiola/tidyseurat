# join_features

join_features() extracts and joins information for specific features

## Usage

``` r
# S4 method for class 'Seurat'
join_features(
  .data,
  features = NULL,
  all = FALSE,
  exclude_zeros = FALSE,
  shape = "wide",
  assay = NULL,
  slot = "data",
  ...
)
```

## Arguments

- .data:

  A tidyseurat object

- features:

  A vector of feature identifiers to join

- all:

  If TRUE return all

- exclude_zeros:

  If TRUE exclude zero values

- shape:

  Format of the returned table "long" or "wide"

- assay:

  assay name to extract feature abundance

- slot:

  slot name to extract feature abundance

- ...:

  Parameters to pass to join wide, i.e. assay name to extract feature
  abundance from and gene prefix, for shape="wide"

## Value

A \`tidyseurat\` object containing information for the specified
features.

## Details

This function extracts information for specified features and returns
the information in either long or wide format.

## Examples

``` r
data(pbmc_small)
pbmc_small %>% join_features(
  features=c("HLA-DRA", "LYZ"))
#> # A Seurat-tibble abstraction: 80 × 17
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
#> # ℹ 10 more variables: RNA_snn_res.1 <fct>, `HLA-DRA` <dbl>, LYZ <dbl>,
#> #   PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>, PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>,
#> #   tSNE_2 <dbl>
```
