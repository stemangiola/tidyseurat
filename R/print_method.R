# This file is a replacement of the unexported functions in the tibble package, in order to specify "tibble abstraction in the header"

# This file is a replacement of the unexported functions in the tibble package, in order to specify "tibble abstraction in the header"

#' Format the header of a tibble
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' For easier customization, the formatting of a tibble is split
#' into three components: header, body, and footer.
#' The `tbl_format_header()` method is responsible for formatting the header
#' of a tibble.
#'
#' Override this method if you need to change the appearance
#' of the entire header.
#' If you only need to change or extend the components shown in the header,
#' override or extend [tbl_sum()] for your class which is called by the
#' default method.
#'
#' @importFrom pillar tbl_format_header
#' @inheritParams tbl_format_body
#' @inherit tbl_format_body return
#'
#' @rdname tbl_format_header-methods
#' @name tbl_format_header
#'
#' @export
#'
NULL

#' @importFrom rlang names2
#' @importFrom pillar align
#' @importFrom pillar get_extent
#' @importFrom pillar style_subtle
#' @export
#' @inheritParams tbl_format_header
tbl_format_header.tidySeurat <- function(x, setup, ...){
  
  number_of_features = x |> attr("number_of_features")
  assay_names = x |> attr("assay_names")
  active_assay = x |> attr("active_assay")
  
  named_header <- setup$tbl_sum
  
  # Change name
  names(named_header) = "A Seurat-tibble abstraction"
  
  if (all(names2(named_header) == "")) {
    header <- named_header
  }
  else {
    header <-
      paste0(
        align(paste0(names2(named_header), ":"), space = NBSP),
        " ",
        named_header
      ) %>%
      
      # Add further info single-cell
      append(sprintf(
        "\033[90m Features=%s | Cells=%s | Active assay=%s | Assays=%s\033[39m",
        number_of_features,
        nrow(x),
        active_assay,
        assay_names %>% paste(collapse=", ")
      ), after = 1)
  }
  style_subtle(pillar___format_comment(header, width = setup$width))
  
}


#' Printing tibbles
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' One of the main features of the `tbl_df` class is the printing:
#'
#' * Tibbles only print as many rows and columns as fit on one screen,
#'   supplemented by a summary of the remaining rows and columns.
#' * Tibble reveals the type of each column, which keeps the user informed about
#'   whether a variable is, e.g., `<chr>` or `<fct>` (character versus factor).
#'
#' Printing can be tweaked for a one-off call by calling `print()` explicitly
#' and setting arguments like `n` and `width`. More persistent control is
#' available by setting the options described below.
#' 
#' Only the first 5 reduced dimensions are displayed, while all of them are queriable (e.g. ggplot). All dimensions are returned/displayed if as_tibble is used.
#'
#' @inheritSection pillar::`pillar-package` Package options
#' @section Package options:
#'
#' The following options are used by the tibble and pillar packages
#' to format and print `tbl_df` objects.
#' Used by the formatting workhorse `trunc_mat()` and, therefore,
#' indirectly, by `print.tbl()`.
#'
#' * `tibble.print_max`: Row number threshold: Maximum number of rows printed.
#'   Set to `Inf` to always print all rows.  Default: 20.
#' * `tibble.print_min`: Number of rows printed if row number threshold is
#'   exceeded. Default: 10.
#' * `tibble.width`: Output width. Default: `NULL` (use `width` option).
#' * `tibble.max_extra_cols`: Number of extra columns printed in reduced form.
#'   Default: 100.
#'
#' @importFrom rlang is_empty
#' @importFrom stringr str_replace
#'
#' @param x Object to format or print.
#' @param ... Other arguments passed on to individual methods.
#' @param n Number of rows to show. If `NULL`, the default, will print all rows
#'   if less than option `tibble.print_max`. Otherwise, will print
#'   `tibble.print_min` rows.
#' @param width Width of text output to generate. This defaults to `NULL`, which
#'   means use `getOption("tibble.width")` or (if also `NULL`)
#'   `getOption("width")`; the latter displays only the columns that fit on one
#'   screen. You can also set `options(tibble.width = Inf)` to override this
#'   default and always print all columns.
#' @param n_extra Number of extra columns to print abbreviated information for,
#'   if the width is too small for the entire tibble. If `NULL`, the default,
#'   will print information about at most `tibble.max_extra_cols` extra columns.
#'   
#' @return Nothing
#'   
#' @examples
#' library(dplyr)
#' data("pbmc_small")
#' pbmc_small  %>% print()
#' @name print
NULL

#' @rdname print
#' @importFrom cli cat_line
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Assays
#' @importFrom vctrs new_data_frame
#' 
#' @export
print.Seurat <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  
  x |>
    as_tibble(n_dimensions_to_return = 5) |>
    new_data_frame(class = c("tidySeurat", "tbl")) %>%
    add_attr( GetAssayData(x) %>% nrow,  "number_of_features") %>%
    add_attr( Assays(x) , "assay_names") %>%
    add_attr( x@active.assay , "active_assay") %>%
    print()
  
  invisible(x)

}

