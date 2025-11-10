# Efficiently bind multiple data frames by row and column

This is an efficient implementation of the common pattern of
\`do.call(rbind, dfs)\` or \`do.call(cbind, dfs)\` for binding many data
frames into one.

This is an efficient implementation of the common pattern of
\`do.call(rbind, dfs)\` or \`do.call(cbind, dfs)\` for binding many data
frames into one.

## Usage

``` r
# S3 method for class 'Seurat'
bind_rows(..., .id = NULL, add.cell.ids = NULL)

# S3 method for class 'Seurat'
bind_cols(..., .id = NULL)
```

## Arguments

- ...:

  Data frames to combine.

  Each argument can either be a data frame, a list that could be a data
  frame, or a list of data frames.

  When row-binding, columns are matched by name, and any missing columns
  will be filled with NA.

  When column-binding, rows are matched by position, so all data frames
  must have the same number of rows. To match by value, not position,
  see mutate-joins.

- .id:

  Data frame identifier.

  When \`.id\` is supplied, a new column of identifiers is created to
  link each row to its original data frame. The labels are taken from
  the named arguments to \`bind_rows()\`. When a list of data frames is
  supplied, the labels are taken from the names of the list. If no names
  are found a numeric sequence is used instead.

- add.cell.ids:

  from Seurat 3.0 A character vector of length(x = c(x, y)). Appends the
  corresponding values to the start of each objects' cell names.

## Value

\`bind_rows()\` and \`bind_cols()\` return the same type as the first
input, either a data frame, \`tbl_df\`, or \`grouped_df\`.

\`bind_rows()\` and \`bind_cols()\` return the same type as the first
input, either a data frame, \`tbl_df\`, or \`grouped_df\`.

## Details

The output of \`bind_rows()\` will contain a column if that column
appears in any of the inputs.

The output of \`bind_rows()\` will contain a column if that column
appears in any of the inputs.

## Examples

``` r
data(pbmc_small)
tt <- pbmc_small
ttservice::bind_rows(tt, tt)
#> Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
#> # A Seurat-tibble abstraction: 160 × 8
#> # Features=230 | Cells=160 | Active assay=RNA | Assays=RNA
#>    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
#>    <chr> <chr>           <dbl>        <int> <chr>           <chr>         <chr> 
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
#> # ℹ 150 more rows
#> # ℹ 1 more variable: RNA_snn_res.1 <chr>

tt_bind <- tt |> select(nCount_RNA ,nFeature_RNA)
#> tidyseurat says: Key columns are missing. A data frame is returned for independent data analysis.
tt |> ttservice::bind_cols(tt_bind)
#> New names:
#> • `nCount_RNA` -> `nCount_RNA...2`
#> • `nFeature_RNA` -> `nFeature_RNA...3`
#> • `nCount_RNA` -> `nCount_RNA...8`
#> • `nFeature_RNA` -> `nFeature_RNA...9`
#> # A Seurat-tibble abstraction: 80 × 17
#> # Features=230 | Cells=80 | Active assay=RNA | Assays=RNA
#>    .cell          orig.ident    nCount_RNA...2 nFeature_RNA...3 RNA_snn_res.0.8
#>    <chr>          <fct>                  <dbl>            <int> <fct>          
#>  1 ATGCCAGAACGACT SeuratProject             70               47 0              
#>  2 CATGGCCTGTGCAT SeuratProject             85               52 0              
#>  3 GAACCTGATGAACC SeuratProject             87               50 1              
#>  4 TGACTGGATTCTCA SeuratProject            127               56 0              
#>  5 AGTCAGACTGCACA SeuratProject            173               53 0              
#>  6 TCTGATACACGTGT SeuratProject             70               48 0              
#>  7 TGGTATCTAAACAG SeuratProject             64               36 0              
#>  8 GCAGCTCTGTTTCT SeuratProject             72               45 0              
#>  9 GATATAACACGCAT SeuratProject             52               36 0              
#> 10 AATGTTGACAGTCA SeuratProject            100               41 0              
#> # ℹ 70 more rows
#> # ℹ 12 more variables: letter.idents <fct>, groups <chr>, RNA_snn_res.1 <fct>,
#> #   nCount_RNA...8 <dbl>, nFeature_RNA...9 <int>, PC_1 <dbl>, PC_2 <dbl>,
#> #   PC_3 <dbl>, PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>
```
