# This file is a replacement of the unexported functions in the tibble
# package, in order to specify "tibble abstraction in the header"

#' @name tbl_format_header
#' @rdname tbl_format_header
#' @inherit pillar::tbl_format_header
#' 
#' @examples
#' # TODO
#' 
#' @importFrom rlang names2
#' @importFrom pillar align
#' @importFrom pillar get_extent
#' @importFrom pillar style_subtle
#' @importFrom pillar tbl_format_header
#' @export
tbl_format_header.tidySeurat <- function(x, setup, ...){
  
    number_of_features <- x |> attr("number_of_features")
    assay_names <- x |> attr("assay_names")
    active_assay <- x |> attr("active_assay")
  
    named_header <- setup$tbl_sum
  
    # Change name
    names(named_header) <- "A Seurat-tibble abstraction"
  
    if (all(names2(named_header) == "")) {
        header <- named_header
    } else {
        header <- paste0(
        align(paste0(names2(named_header), ":"), space=NBSP),
        " ", named_header) %>%
        # Add further info single-cell
        append(sprintf(
            "\033[90m Features=%s | Cells=%s | Active assay=%s | Assays=%s\033[39m",
            number_of_features,
            nrow(x),
            active_assay,
            assay_names %>% paste(collapse=", ")
        ), after=1)
    }
    style_subtle(pillar___format_comment(header, width=setup$width))
}


#' @name formatting
#' @rdname formatting
#' @aliases print
#' @inherit tibble::formatting
#' @return Prints a message to the console describing
#'   the contents of the `tidyseurat`.
#'
#' @param ... Passed on to \code{\link[pillar:tbl_format_setup]{tbl_format_setup()}}.
#' @param n_extra Number of extra columns to print abbreviated information for,
#'   if the width is too small for the entire tibble. If `NULL`, the default,
#'   will print information about at most `tibble.max_extra_cols` extra columns.
#' @examples
#' data(pbmc_small)
#' print(pbmc_small)
#' 
#' @importFrom vctrs new_data_frame
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat Assays
#' @export
print.Seurat <- function(x, ..., n=NULL, width=NULL, n_extra=NULL) {

    x |>
        as_tibble(n_dimensions_to_return=5) |>
        new_data_frame(class=c("tidySeurat", "tbl")) %>%
        add_attr(GetAssayData(x) %>% nrow,  "number_of_features") %>%
        add_attr(Assays(x) , "assay_names") %>%
        add_attr(x@active.assay , "active_assay") %>%
        print()
    invisible(x)
}

