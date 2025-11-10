# Get a glimpse of your data

`glimpse()` is like a transposed version of [`print()`](formatting.md):
columns run down the page, and data runs across. This makes it possible
to see every column in a data frame. It's a little like
[`str()`](https://rdrr.io/r/utils/str.html) applied to a data frame but
it tries to show you as much data as possible. (And it always shows the
underlying data, even when applied to a remote data source.)

See
[`format_glimpse()`](https://pillar.r-lib.org/reference/format_glimpse.html)
for details on the formatting.

## Usage

``` r
# S3 method for class 'tidyseurat'
glimpse(x, width = NULL, ...)
```

## Arguments

- x:

  An object to glimpse at.

- width:

  Width of output: defaults to the setting of the `width`
  [option](https://pillar.r-lib.org/reference/pillar_options.html) (if
  finite) or the width of the console.

- ...:

  Unused, for extensibility.

## Value

x original x is (invisibly) returned, allowing `glimpse()` to be used
within a data pipe line.

## S3 methods

`glimpse` is an S3 generic with a customised method for `tbl`s and
`data.frames`, and a default method that calls
[`str()`](https://rdrr.io/r/utils/str.html).

## Examples

``` r
data(pbmc_small)
pbmc_small |> glimpse()
#> Formal class 'Seurat' [package "SeuratObject"] with 13 slots
#>   ..@ assays      :List of 1
#>   .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
#>   ..@ meta.data   :'data.frame': 80 obs. of  7 variables:
#>   .. ..$ orig.ident     : Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..$ nCount_RNA     : num [1:80] 70 85 87 127 173 70 64 72 52 100 ...
#>   .. ..$ nFeature_RNA   : int [1:80] 47 52 50 56 53 48 36 45 36 41 ...
#>   .. ..$ RNA_snn_res.0.8: Factor w/ 2 levels "0","1": 1 1 2 1 1 1 1 1 1 1 ...
#>   .. ..$ letter.idents  : Factor w/ 2 levels "A","B": 1 1 2 1 1 1 1 1 1 1 ...
#>   .. ..$ groups         : chr [1:80] "g2" "g1" "g2" "g2" ...
#>   .. ..$ RNA_snn_res.1  : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>   ..@ active.assay: chr "RNA"
#>   ..@ active.ident: Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..- attr(*, "names")= chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
#>   ..@ graphs      :List of 1
#>   .. ..$ RNA_snn:Formal class 'Graph' [package "SeuratObject"] with 7 slots
#>   ..@ neighbors   : list()
#>   ..@ reductions  :List of 2
#>   .. ..$ pca :Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
#>   .. ..$ tsne:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
#>   ..@ images      : list()
#>   ..@ project.name: chr "SeuratProject"
#>   ..@ misc        : list()
#>   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
#>   .. ..$ : int [1:3] 4 0 0
#>   ..@ commands    :List of 10
#>   .. ..$ NormalizeData.RNA       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ ScaleData.RNA           :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ RunPCA.RNA              :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ BuildSNN.RNA.pca        :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ FindClusters            :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ RunTSNE.pca             :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ JackStraw.RNA.pca       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ ScoreJackStraw.pca      :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ ProjectDim.RNA.pca      :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   .. ..$ FindVariableFeatures.RNA:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
#>   ..@ tools       : list()
```
