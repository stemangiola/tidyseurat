# This file is a replacement of the unexported functions in the tibble package, in order to specify "tibble abstraction in the header"

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
#' pbmc_small  %>% print()
#' @name print
NULL

#' @rdname print
#' @importFrom cli cat_line
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Assays
#' 
#' @export
print.Seurat <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  
  x %>%
    as_tibble(n_dimensions_to_return = 5) %>%
    
    # Get formatting
    tidyseurat_format_tbl(..., n = n, width = width, n_extra = n_extra) %>%
  
    # Hijack the tibble header
    map_chr(~ .x %>% str_replace("A tibble:", "A Seurat-tibble abstraction:")) %>%
    
    # Insert more info
    append(sprintf(
      "\033[90m# Transcripts=%s | Active assay=%s | Assays=%s\033[39m", 
      GetAssayData(x) %>% nrow,
      x@active.assay,
      Assays(x) %>% paste(collapse=", ")
      ), 
      after = 1
    ) %>%
    
    # Output
    cat_line()
  
  invisible(x)
}

#' @importFrom tibble trunc_mat
tidyseurat_format_tbl <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  mat <- trunc_mat(x, n = n, width = width, n_extra = n_extra)
  tidyseurat_format_truncated_mat(mat)
}

tidyseurat_pluralise_n <- function(message, n) {
  stopifnot(n >= 0)
  
  
  # Don't strip parens if they have a space in between
  # (but not if the space comes before the closing paren)
  
  if (n == 1) {
    # strip [
    message <- gsub("\\[([^\\] ]* *)\\]", "\\1", message, perl = TRUE)
    # remove ( and its content
    message <- gsub("\\([^\\) ]* *\\)", "", message, perl = TRUE)
  } else {
    # strip (
    message <- gsub("\\(([^\\) ]* *)\\)", "\\1", message, perl = TRUE)
    # remove [ and its content
    message <- gsub("\\[[^\\] ]* *\\]", "", message, perl = TRUE)
  }
  
  message
}

tidyseurat_nchar_width <- function(x) {
  nchar(x, type = "width")
}

#' @importFrom pillar style_subtle
#' @importFrom rlang names2
#' @importFrom pillar squeeze
tidyseurat_format_truncated_mat <- function(x, width = NULL, ...) {
  if (is.null(width)) {
    width <- x$width
  }
  
  width <- tidyseurat_tibble_width(width)
  
  named_header <- tidyseurat_format_header(x)
  if (all(names2(named_header) == "")) {
    header <- named_header
  } else {
    header <- paste0(
      tidyseurat_justify(
        paste0(names2(named_header), ":"),
        right = FALSE, space = NBSP
      ),
      # We add a space after the NBSP inserted by tidyseurat_justify()
      # so that wrapping occurs at the right location for very narrow outputs
      " ",
      named_header
    )
  }
  
  comment <- tidyseurat_format_comment(header, width = width)
  squeezed <- squeeze(x$mcf, width = width)
  mcf <- tidyseurat_format_body(squeezed)
  
  # Splitting lines is important, otherwise subtle style may be lost
  # if column names contain spaces.
  footer <- tidyseurat_pre_dots(tidyseurat_format_footer(x, squeezed))
  footer_comment <- tidyseurat_split_lines(tidyseurat_format_comment(footer, width = width))
  
  c(style_subtle(comment), mcf, style_subtle(footer_comment))
}

tidyseurat_format_header <- function(x) {
  x$summary
}

tidyseurat_format_body <- function(x) {
  format(x)
}

#' @importFrom pillar extra_cols
tidyseurat_format_footer <- function(x, squeezed_colonnade) {
  extra_rows <- tidyseurat_format_footer_rows(x)
  extra_cols <- tidyseurat_format_footer_cols(x, extra_cols(squeezed_colonnade, n = x$n_extra))
  
  extra <- c(extra_rows, extra_cols)
  if (length(extra) >= 1) {
    extra[[1]] <- paste0("with ", extra[[1]])
    extra[-1] <- map_chr(extra[-1], function(ex) paste0("and ", ex))
    collapse(extra)
  } else {
    character()
  }
}

tidyseurat_format_footer_rows <- function(x) {
  if (length(x$mcf) != 0) {
    if (is.na(x$rows_missing)) {
      "more rows"
    } else if (x$rows_missing > 0) {
      paste0(tidyseurat_big_mark(x$rows_missing), tidyseurat_pluralise_n(" more row(s)", x$rows_missing))
    }
  } else if (is.na(x$rows_total) && x$rows_min > 0) {
    paste0("at least ", tidyseurat_big_mark(x$rows_min), tidyseurat_pluralise_n(" row(s) total", x$rows_min))
  }
}

tidyseurat_format_footer_cols <- function(x, extra_cols) {
  if (length(extra_cols) == 0) return(NULL)
  
  vars <- tidyseurat_format_extra_vars(extra_cols)
  paste0(
    tidyseurat_big_mark(length(extra_cols)), " ",
    if (!identical(x$rows_total, 0L) && x$rows_min > 0) "more ",
    pluralise("variable(s)", extra_cols), vars
  )
}

tidyseurat_format_extra_vars <- function(extra_cols) {
  # Also covers empty extra_cols vector!
  if (is.na(extra_cols[1])) return("")
  
  if (anyNA(extra_cols)) {
    extra_cols <- c(extra_cols[!is.na(extra_cols)], cli::symbol$ellipsis)
  }
  
  paste0(": ", collapse(extra_cols))
}

tidyseurat_format_comment <- function(x, width) {
  if (length(x) == 0L) return(character())
  map_chr(x, tidyseurat_wrap, prefix = "# ", width = min(width, getOption("width")))
}

tidyseurat_pre_dots <- function(x) {
  if (length(x) > 0) {
    paste0(cli::symbol$ellipsis, " ", x)
  } else {
    character()
  }
}

tidyseurat_justify <- function(x, right = TRUE, space = " ") {
  if (length(x) == 0L) return(character())
  width <- tidyseurat_nchar_width(x)
  max_width <- max(width)
  spaces_template <- paste(rep(space, max_width), collapse = "")
  spaces <- map_chr(max_width - width, substr, x = spaces_template, start = 1L)
  if (right) {
    paste0(spaces, x)
  } else {
    paste0(x, spaces)
  }
}

tidyseurat_split_lines <- function(x) {
  # Avoid .ptype argument to vec_c()
  if (is_empty(x)) return(character())
  
  unlist(strsplit(x, "\n", fixed = TRUE))
}

tidyseurat_big_mark <- function(x, ...) {
  # The thousand separator,
  # "," unless it's used for the decimal point, in which case "."
  mark <- if (identical(getOption("OutDec"), ",")) "." else ","
  ret <- formatC(x, big.mark = mark, format = "d", ...)
  ret[is.na(x)] <- "??"
  ret
}

collapse <- function(x) paste(x, collapse = ", ")

# tidyseurat_wrap --------------------------------------------------------------------

NBSP <- "\U00A0"

tidyseurat_wrap <- function(..., indent = 0, prefix = "", width) {
  x <- paste0(..., collapse = "")
  wrapped <- tidyseurat_strwrap2(x, width - tidyseurat_nchar_width(prefix), indent)
  wrapped <- paste0(prefix, wrapped)
  wrapped <- gsub(NBSP, " ", wrapped)
  
  paste0(wrapped, collapse = "\n")
}

#' @importFrom fansi strwrap_ctl
tidyseurat_strwrap2 <- function(x, width, indent) {
  strwrap_ctl(x, width = max(width, 0), indent = indent, exdent = indent + 2)
}


op.tibble <- list(
  tibble.print_max = 20L,
  tibble.print_min = 10L,
  tibble.width = NULL,
  tibble.max_extra_cols = 100L,
  tibble.view_max = 1000L
)

tidyseurat_tibble_opt <- function(x, dplyr = TRUE) {
  x_tibble <- paste0("tibble.", x)
  res <- getOption(x_tibble)
  if (!is.null(res)) {
    return(res)
  }
  
  if (dplyr) {
    x_dplyr <- paste0("dplyr.", x)
    res <- getOption(x_dplyr)
    if (!is.null(res)) {
      return(res)
    }
  }
  
  op.tibble[[x_tibble]]
}

tidyseurat_tibble_width <- function(width) {
  if (!is.null(width)) {
    return(width)
  }
  
  width <- tidyseurat_tibble_opt("width")
  if (!is.null(width)) {
    return(width)
  }
  
  getOption("width")
}

tidyseurat_tibble_glimpse_width <- function(width) {
  if (!is.null(width)) {
    return(width)
  }
  
  width <- tidyseurat_tibble_opt("width")
  if (!is.null(width) && is.finite(width)) {
    return(width)
  }
  
  getOption("width")
}

# Pluralise if there are many columns

pluralise <- function(message, objects) {
  pluralise_n(message, length(objects))
}

pluralise_n <- function(message, n) {
  stopifnot(n >= 0)
  
  
  # Don't strip parens if they have a space in between
  # (but not if the space comes before the closing paren)
  
  if (n == 1) {
    # strip [
    message <- gsub("\\[([^\\] ]* *)\\]", "\\1", message, perl = TRUE)
    # remove ( and its content
    message <- gsub("\\([^\\) ]* *\\)", "", message, perl = TRUE)
  } else {
    # strip (
    message <- gsub("\\(([^\\) ]* *)\\)", "\\1", message, perl = TRUE)
    # remove [ and its content
    message <- gsub("\\[[^\\] ]* *\\]", "", message, perl = TRUE)
  }
  
  message
}

