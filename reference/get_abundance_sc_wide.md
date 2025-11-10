# get abundance wide

get abundance wide

## Usage

``` r
get_abundance_sc_wide(
  .data,
  features = NULL,
  all = FALSE,
  assay = .data@active.assay,
  slot = "data",
  prefix = ""
)
```

## Arguments

- .data:

  A tidyseurat

- features:

  A character

- all:

  A boolean

- assay:

  assay name to extract feature abundance

- slot:

  slot in the assay, e.g. \`data\` and \`scale.data\`

- prefix:

  prefix for the feature names

## Value

A Seurat object

## Examples

``` r
data(pbmc_small)
pbmc_small %>%
  get_abundance_sc_wide(features=c("HLA-DRA", "LYZ"))
#> # A tibble: 80 × 3
#>    .cell          `HLA-DRA`   LYZ
#>    <chr>              <dbl> <dbl>
#>  1 ATGCCAGAACGACT      0     4.97
#>  2 CATGGCCTGTGCAT      4.78  4.78
#>  3 GAACCTGATGAACC      0     4.75
#>  4 TGACTGGATTCTCA      0     0   
#>  5 AGTCAGACTGCACA      4.07  0   
#>  6 TCTGATACACGTGT      4.97  4.97
#>  7 TGGTATCTAAACAG      0     0   
#>  8 GCAGCTCTGTTTCT      4.94  0   
#>  9 GATATAACACGCAT      0     5.26
#> 10 AATGTTGACAGTCA      0     0   
#> # ℹ 70 more rows
```
