
#' Create a new ggplot from a tidyseurat object
#'
#' `ggplot()` initializes a ggplot object. It can be used to
#' declare the input data frame for a graphic and to specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' `ggplot()` is used to construct the initial plot object,
#' and is almost always followed by `+` to add component to the
#' plot. There are three common ways to invoke `ggplot()`:
#'
#' * `ggplot(df, aes(x, y, other aesthetics))`
#' * `ggplot(df)`
#' * `ggplot()`
#'
#' The first method is recommended if all layers use the same
#' data and the same set of aesthetics, although this method
#' can also be used to add a layer using data from another
#' data frame. See the first example below. The second
#' method specifies the default data frame to use for the plot,
#' but no aesthetics are defined up front. This is useful when
#' one data frame is used predominantly as layers are added,
#' but the aesthetics may vary from one layer to another. The
#' third method initializes a skeleton `ggplot` object which
#' is fleshed out as layers are added. This method is useful when
#' multiple data frames are used to produce different layers, as
#' is often the case in complex graphics.
#'
#' @param data Default dataset to use for plot. If not already a data.frame,
#'   will be converted to one by [fortify()]. If not specified,
#'   must be supplied in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot.
#'   If not specified, must be supplied in each layer added to the plot.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param environment DEPRECATED. Used prior to tidy evaluation.
#' @export
#' @examples
#' # Generate some sample data, then compute mean and standard deviation
#' # in each group
#' df <- data.frame(
#'   gp = factor(rep(letters[1:3], each = 10)),
#'   y = rnorm(30)
#' )
#' ds <- do.call(rbind, lapply(split(df, df$gp), function(d) {
#'   data.frame(mean = mean(d$y), sd = sd(d$y), gp = d$gp)
#' }))
#'
#' # The summary data frame ds is used to plot larger red points on top
#' # of the raw data. Note that we don't need to supply `data` or `mapping`
#' # in each layer because the defaults from ggplot() are used.
#' ggplot(df, aes(gp, y)) +
#'   geom_point() +
#'   geom_point(data = ds, aes(y = mean), colour = 'red', size = 3)
#'
ggplot <- function(.data = NULL, mapping = aes(), ..., environment = parent.frame()) {
  UseMethod("ggplot")
}

#' @export
#' 
ggplot.tbl_df <- function(.data = NULL, mapping = aes(), ..., environment = parent.frame()) {
  
  .data %>%
    
    # This is a trick to not loop the call
    drop_class("tbl_df") %>%
    ggplot2::ggplot( mapping = mapping, ..., environment = environment)
  
  
}

#' @export
ggplot.tidyseurat <- function(.data = NULL, mapping = aes(), ..., environment = parent.frame()) {
  .data %>%
    as_tibble() %>%
    ggplot2::ggplot( mapping = mapping)
}

