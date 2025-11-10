# get abundance long

get abundance long

## Usage

``` r
get_abundance_sc_long(
  .data,
  features = NULL,
  all = FALSE,
  exclude_zeros = FALSE,
  assay = Assays(.data),
  slot = "data"
)
```

## Arguments

- .data:

  A tidyseurat

- features:

  A character

- all:

  A boolean

- exclude_zeros:

  A boolean

- assay:

  assay name to extract feature abundance

- slot:

  slot in the assay, e.g. \`data\` and \`scale.data\`

## Value

A Seurat object

## Examples

``` r
data(pbmc_small)
pbmc_small %>%
  get_abundance_sc_long(features=c("HLA-DRA", "LYZ"))
#> # A tibble: 160 × 3
#>    .feature .cell          .abundance_RNA
#>    <chr>    <chr>                   <dbl>
#>  1 HLA-DRA  ATGCCAGAACGACT           0   
#>  2 HLA-DRA  CATGGCCTGTGCAT           4.78
#>  3 HLA-DRA  GAACCTGATGAACC           0   
#>  4 HLA-DRA  TGACTGGATTCTCA           0   
#>  5 HLA-DRA  AGTCAGACTGCACA           4.07
#>  6 HLA-DRA  TCTGATACACGTGT           4.97
#>  7 HLA-DRA  TGGTATCTAAACAG           0   
#>  8 HLA-DRA  GCAGCTCTGTTTCT           4.94
#>  9 HLA-DRA  GATATAACACGCAT           0   
#> 10 HLA-DRA  AATGTTGACAGTCA           0   
#> # ℹ 150 more rows
```
