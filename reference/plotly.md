# Initiate a plotly visualization

This function maps R objects to
[plotly.js](https://plotly.com/javascript/), an (MIT licensed) web-based
interactive charting library. It provides abstractions for doing common
things (e.g. mapping data values to fill colors (via `color`) or
creating [animation](https://rdrr.io/pkg/plotly/man/animation.html)s
(via `frame`)) and sets some different defaults to make the interface
feel more 'R-like' (i.e., closer to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`ggplot2::qplot()`](https://ggplot2.tidyverse.org/reference/qplot.html)).

## Usage

``` r
plot_ly(
  data = data.frame(),
  ...,
  type = NULL,
  name = NULL,
  color = NULL,
  colors = NULL,
  alpha = NULL,
  stroke = NULL,
  strokes = NULL,
  alpha_stroke = 1,
  size = NULL,
  sizes = c(10, 100),
  span = NULL,
  spans = c(1, 20),
  symbol = NULL,
  symbols = NULL,
  linetype = NULL,
  linetypes = NULL,
  split = NULL,
  frame = NULL,
  width = NULL,
  height = NULL,
  source = "A"
)

# S3 method for class 'tbl_df'
plot_ly(
  data = data.frame(),
  ...,
  type = NULL,
  name = NULL,
  color = NULL,
  colors = NULL,
  alpha = NULL,
  stroke = NULL,
  strokes = NULL,
  alpha_stroke = 1,
  size = NULL,
  sizes = c(10, 100),
  span = NULL,
  spans = c(1, 20),
  symbol = NULL,
  symbols = NULL,
  linetype = NULL,
  linetypes = NULL,
  split = NULL,
  frame = NULL,
  width = NULL,
  height = NULL,
  source = "A"
)

# S3 method for class 'Seurat'
plot_ly(
  data = data.frame(),
  ...,
  type = NULL,
  name = NULL,
  color = NULL,
  colors = NULL,
  alpha = NULL,
  stroke = NULL,
  strokes = NULL,
  alpha_stroke = 1,
  size = NULL,
  sizes = c(10, 100),
  span = NULL,
  spans = c(1, 20),
  symbol = NULL,
  symbols = NULL,
  linetype = NULL,
  linetypes = NULL,
  split = NULL,
  frame = NULL,
  width = NULL,
  height = NULL,
  source = "A"
)
```

## Arguments

- data:

  A data frame (optional) or
  [crosstalk::SharedData](https://rdrr.io/pkg/crosstalk/man/SharedData.html)
  object.

- ...:

  Arguments (i.e., attributes) passed along to the trace `type`. See
  [`schema()`](https://rdrr.io/pkg/plotly/man/schema.html) for a list of
  acceptable attributes for a given trace `type` (by going to `traces`
  -\> `type` -\> `attributes`). Note that attributes provided at this
  level may override other arguments (e.g.
  `plot_ly(x = 1:10, y = 1:10, color = I("red"), marker = list(color = "blue"))`).

- type:

  A character string specifying the trace type (e.g. `"scatter"`,
  `"bar"`, `"box"`, etc). If specified, it *always* creates a trace,
  otherwise

- name:

  Values mapped to the trace's name attribute. Since a trace can only
  have one name, this argument acts very much like `split` in that it
  creates one trace for every unique value.

- color:

  Values mapped to relevant 'fill-color' attribute(s) (e.g.
  [fillcolor](https://plotly.com/r/reference/#scatter-fillcolor),
  [marker.color](https://plotly.com/r/reference/#scatter-marker-color),
  [textfont.color](https://plotly.com/r/reference/#scatter-textfont-color),
  etc.). The mapping from data values to color codes may be controlled
  using `colors` and `alpha`, or avoided altogether via
  [`I()`](https://rdrr.io/r/base/AsIs.html) (e.g., `color = I("red")`).
  Any color understood by
  [`grDevices::col2rgb()`](https://rdrr.io/r/grDevices/col2rgb.html) may
  be used in this way.

- colors:

  Either a colorbrewer2.org palette name (e.g. "YlOrRd" or "Blues"), or
  a vector of colors to interpolate in hexadecimal "#RRGGBB" format, or
  a color interpolation function like
  [`colorRamp()`](https://rdrr.io/r/grDevices/colorRamp.html).

- alpha:

  A number between 0 and 1 specifying the alpha channel applied to
  `color`. Defaults to 0.5 when mapping to
  [fillcolor](https://plotly.com/r/reference/#scatter-fillcolor) and 1
  otherwise.

- stroke:

  Similar to `color`, but values are mapped to relevant 'stroke-color'
  attribute(s) (e.g.,
  [marker.line.color](https://plotly.com/r/reference/#scatter-marker-line-color)
  and [line.color](https://plotly.com/r/reference/#scatter-line-color)
  for filled polygons). If not specified, `stroke` inherits from
  `color`.

- strokes:

  Similar to `colors`, but controls the `stroke` mapping.

- alpha_stroke:

  Similar to `alpha`, but applied to `stroke`.

- size:

  (Numeric) values mapped to relevant 'fill-size' attribute(s) (e.g.,
  [marker.size](https://plotly.com/r/reference/#scatter-marker-size),
  [textfont.size](https://plotly.com/r/reference/#scatter-textfont-size),
  and
  [error_x.width](https://plotly.com/r/reference/#scatter-error_x-width)).
  The mapping from data values to symbols may be controlled using
  `sizes`, or avoided altogether via
  [`I()`](https://rdrr.io/r/base/AsIs.html) (e.g., `size = I(30)`).

- sizes:

  A numeric vector of length 2 used to scale `size` to pixels.

- span:

  (Numeric) values mapped to relevant 'stroke-size' attribute(s) (e.g.,
  [marker.line.width](https://plotly.com/r/reference/#scatter-marker-line-width),
  [line.width](https://plotly.com/r/reference/#scatter-line-width) for
  filled polygons, and
  [error_x.thickness](https://plotly.com/r/reference/#scatter-error_x-thickness))
  The mapping from data values to symbols may be controlled using
  `spans`, or avoided altogether via
  [`I()`](https://rdrr.io/r/base/AsIs.html) (e.g., `span = I(30)`).

- spans:

  A numeric vector of length 2 used to scale `span` to pixels.

- symbol:

  (Discrete) values mapped to
  [marker.symbol](https://plotly.com/r/reference/#scatter-marker-symbol).
  The mapping from data values to symbols may be controlled using
  `symbols`, or avoided altogether via
  [`I()`](https://rdrr.io/r/base/AsIs.html) (e.g.,
  `symbol = I("pentagon")`). Any
  [pch](https://rdrr.io/r/graphics/points.html) value or [symbol
  name](https://plotly.com/r/reference/#scatter-marker-symbol) may be
  used in this way.

- symbols:

  A character vector of [pch](https://rdrr.io/r/graphics/points.html)
  values or [symbol
  names](https://plotly.com/r/reference/#scatter-marker-symbol).

- linetype:

  (Discrete) values mapped to
  [line.dash](https://plotly.com/r/reference/#scatter-line-dash). The
  mapping from data values to symbols may be controlled using
  `linetypes`, or avoided altogether via
  [`I()`](https://rdrr.io/r/base/AsIs.html) (e.g.,
  `linetype = I("dash")`). Any `lty` (see
  [par](https://rdrr.io/r/graphics/par.html)) value or [dash
  name](https://plotly.com/r/reference/#scatter-line-dash) may be used
  in this way.

- linetypes:

  A character vector of `lty` values or [dash
  names](https://plotly.com/r/reference/#scatter-line-dash)

- split:

  (Discrete) values used to create multiple traces (one trace per
  value).

- frame:

  (Discrete) values used to create animation frames.

- width:

  Width in pixels (optional, defaults to automatic sizing).

- height:

  Height in pixels (optional, defaults to automatic sizing).

- source:

  a character string of length 1. Match the value of this string with
  the source argument in
  [`event_data()`](https://rdrr.io/pkg/plotly/man/event_data.html) to
  retrieve the event data corresponding to a specific plot (shiny apps
  can have multiple plots).

## Value

\`plotly\`

## Details

Unless `type` is specified, this function just initiates a plotly object
with 'global' attributes that are passed onto downstream uses of
[`add_trace()`](https://rdrr.io/pkg/plotly/man/add_trace.html) (or
similar). A [formula](https://rdrr.io/r/stats/formula.html) must always
be used when referencing column name(s) in `data` (e.g.
`plot_ly(mtcars, x = ~wt)`). Formulas are optional when supplying values
directly, but they do help inform default axis/scale titles (e.g.,
`plot_ly(x = mtcars$wt)` vs `plot_ly(x = ~mtcars$wt)`)

## References

<https://plotly.com/r/>

## See also

- For initializing a plotly-geo object:
  [`plot_geo()`](https://rdrr.io/pkg/plotly/man/plot_geo.html)

- For initializing a plotly-mapbox object:
  [`plot_mapbox()`](https://rdrr.io/pkg/plotly/man/plot_mapbox.html)

- For translating a ggplot2 object to a plotly object:
  [`ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)

- For modifying any plotly object:
  [`layout()`](https://rdrr.io/r/graphics/layout.html),
  [`add_trace()`](https://rdrr.io/pkg/plotly/man/add_trace.html),
  [`style()`](https://rdrr.io/pkg/plotly/man/style.html)

- For linked brushing:
  [`highlight()`](https://rdrr.io/pkg/plotly/man/highlight.html)

- For arranging multiple plots:
  [`subplot()`](https://rdrr.io/pkg/plotly/man/subplot.html),
  [`crosstalk::bscols()`](https://rdrr.io/pkg/crosstalk/man/bscols.html)

- For inspecting plotly objects:
  [`plotly_json()`](https://rdrr.io/pkg/plotly/man/plotly_json.html)

- For quick, accurate, and searchable plotly.js reference:
  [`schema()`](https://rdrr.io/pkg/plotly/man/schema.html)

## Author

Carson Sievert

## Examples

``` r
data(pbmc_small)
plot_ly(pbmc_small)

{"x":{"visdat":{"22143657d2f2":["function () ","plotlyVisDat"]},"cur_data":"22143657d2f2","attrs":{"22143657d2f2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20]}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true},"yaxis":{"domain":[0,1],"automargin":true},"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"type":"scatter","mode":"markers","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
