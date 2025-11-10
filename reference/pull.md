# Extract a single column

`pull()` is similar to `$`. It's mostly useful because it looks a little
nicer in pipes, it also works with remote data frames, and it can
optionally name the output.

## Usage

``` r
# S3 method for class 'Seurat'
pull(.data, var = -1, name = NULL, ...)
```

## Arguments

- .data:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for more
  details.

- var:

  A variable specified as:

  - a literal variable name

  - a positive integer, giving the position counting from the left

  - a negative integer, giving the position counting from the right.

  The default returns the last column (on the assumption that's the
  column you've created most recently).

  This argument is taken by expression and supports
  [quasiquotation](https://rlang.r-lib.org/reference/topic-inject.html)
  (you can unquote column names and column locations).

- name:

  An optional parameter that specifies the column to be used as names
  for a named vector. Specified in a similar manner as `var`.

- ...:

  For use by methods.

## Value

A vector the same size as `.data`.

## Methods

This function is a **generic**, which means that packages can provide
implementations (methods) for other classes. See the documentation of
individual methods for extra arguments and differences in behaviour.

The following methods are currently available in loaded packages: dplyr
(`data.frame`), tidyseurat (`Seurat`) .

## Examples

``` r
data(pbmc_small)
pbmc_small |> pull(groups)
#> tidyseurat says: A data frame is returned for independent data analysis.
#>  [1] "g2" "g1" "g2" "g2" "g2" "g1" "g1" "g1" "g1" "g1" "g2" "g1" "g2" "g2" "g2"
#> [16] "g1" "g2" "g1" "g1" "g2" "g1" "g1" "g2" "g2" "g1" "g2" "g2" "g2" "g2" "g1"
#> [31] "g1" "g1" "g1" "g2" "g1" "g1" "g2" "g1" "g1" "g2" "g1" "g2" "g2" "g2" "g1"
#> [46] "g2" "g1" "g2" "g1" "g2" "g1" "g2" "g2" "g2" "g1" "g1" "g1" "g1" "g2" "g1"
#> [61] "g1" "g1" "g1" "g1" "g1" "g2" "g2" "g1" "g1" "g1" "g2" "g1" "g2" "g2" "g1"
#> [76] "g1" "g2" "g1" "g2" "g1"
```
