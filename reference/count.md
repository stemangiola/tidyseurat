# Count observations by group

\`count()\` lets you quickly count the unique values of one or more
variables: \`df \`df \`count()\` is paired with \`tally()\`, a
lower-level helper that is equivalent to \`df switching the summary from
\`n = n()\` to \`n = sum(wt)\`.

\`add_count()\` and \`add_tally()\` are equivalents to \`count()\` and
\`tally()\` but use \`mutate()\` instead of \`summarise()\` so that they
add a new column with group-wise counts.

## Usage

``` r
# S3 method for class 'Seurat'
count(
  x,
  ...,
  wt = NULL,
  sort = FALSE,
  name = NULL,
  .drop = group_by_drop_default(x)
)

add_count(x, ..., wt = NULL, sort = FALSE, name = NULL)

# Default S3 method
add_count(x, ..., wt = NULL, sort = FALSE, name = NULL)

# S3 method for class 'Seurat'
add_count(x, ..., wt = NULL, sort = FALSE, name = NULL)
```

## Arguments

- x:

  A data frame, data frame extension (e.g. a tibble), or a lazy data
  frame (e.g. from dbplyr or dtplyr).

- ...:

  \<\[\`data-masking\`\]\[dplyr_data_masking\]\> Variables to group by.

- wt:

  \<\[\`data-masking\`\]\[dplyr_data_masking\]\> Frequency weights. Can
  be \`NULL\` or a variable:

  \* If \`NULL\` (the default), counts the number of rows in each group.
  \* If a variable, computes \`sum(wt)\` for each group.

- sort:

  If \`TRUE\`, will show the largest groups at the top.

- name:

  The name of the new column in the output.

  If omitted, it will default to \`n\`. If there's already a column
  called \`n\`, it will error, and require you to specify the name.

- .drop:

  For \`count()\`: if \`FALSE\` will include counts for empty groups
  (i.e. for levels of factors that don't exist in the data).

## Value

An object of the same type as \`.data\`. \`count()\` and \`add_count()\`
group transiently, so the output has the same groups as the input.

## Examples

``` r
data(pbmc_small)
pbmc_small |> count(groups)
#> tidyseurat says: A data frame is returned for independent data analysis.
#> # A tibble: 2 × 2
#>   groups     n
#>   <chr>  <int>
#> 1 g1        44
#> 2 g2        36
    
```
