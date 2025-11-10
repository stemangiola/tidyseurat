# Create a new `ggplot` from a `tidyseurat`

`ggplot()` initializes a ggplot object. It can be used to declare the
input data frame for a graphic and to specify the set of plot aesthetics
intended to be common throughout all subsequent layers unless
specifically overridden.

## Usage

``` r
# S3 method for class 'Seurat'
ggplot(data = NULL, mapping = aes(), ..., environment = parent.frame())
```

## Arguments

- data:

  Default dataset to use for plot. If not already a data.frame, will be
  converted to one by
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html).
  If not specified, must be supplied in each layer added to the plot.

- mapping:

  Default list of aesthetic mappings to use for plot. If not specified,
  must be supplied in each layer added to the plot.

- ...:

  Other arguments passed on to methods. Not currently used.

- environment:

  **\[deprecated\]** Used prior to tidy evaluation.

## Value

\`ggplot\`

## Details

`ggplot()` is used to construct the initial plot object, and is almost
always followed by a plus sign (`+`) to add components to the plot.

There are three common patterns used to invoke `ggplot()`:

- `ggplot(data = df, mapping = aes(x, y, other aesthetics))`

- `ggplot(data = df)`

- `ggplot()`

The first pattern is recommended if all layers use the same data and the
same set of aesthetics, although this method can also be used when
adding a layer using data from another data frame.

The second pattern specifies the default data frame to use for the plot,
but no aesthetics are defined up front. This is useful when one data
frame is used predominantly for the plot, but the aesthetics vary from
one layer to another.

The third pattern initializes a skeleton `ggplot` object, which is
fleshed out as layers are added. This is useful when multiple data
frames are used to produce different layers, as is often the case in
complex graphics.

The `data =` and `mapping =` specifications in the arguments are
optional (and are often omitted in practice), so long as the data and
the mapping values are passed into the function in the right order. In
the examples below, however, they are left in place for clarity.

## See also

The [first steps chapter](https://ggplot2-book.org/getting-started) of
the online ggplot2 book.

## Examples

``` r
library(ggplot2)
data(pbmc_small)
pbmc_small |> 
  ggplot(aes(groups, nCount_RNA)) +
  geom_boxplot()

```
