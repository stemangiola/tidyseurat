# (DEPRECATED) Extract and join information for transcripts.

join_transcripts() extracts and joins information for specified
transcripts

## Usage

``` r
join_transcripts(
  .data,
  transcripts = NULL,
  all = FALSE,
  exclude_zeros = FALSE,
  shape = "wide",
  ...
)
```

## Arguments

- .data:

  A tidyseurat object

- transcripts:

  A vector of transcript identifiers to join

- all:

  If TRUE return all

- exclude_zeros:

  If TRUE exclude zero values

- shape:

  Format of the returned table "long" or "wide"

- ...:

  Parameters to pass to join wide, i.e. assay name to extract transcript
  abundance from

## Value

A \`tbl\` containing the information.for the specified transcripts

## Details

DEPRECATED, please use join_features()

## Examples

``` r

print("DEPRECATED")
#> [1] "DEPRECATED"

```
