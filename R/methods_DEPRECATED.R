#' (DEPRECATED) Extract and join information for transcripts.
#'
#'
#' @description join_transcripts() extracts and joins information for specified transcripts
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_transcripts
#' @rdname join_transcripts
#'
#' @param .data A tidyseurat object
#' @param transcripts A vector of transcript identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract transcript abundance from
#'
#' @details DEPRECATED, please use join_features()
#'
#' @return A `tbl` containing the information.for the specified transcripts
#' 
#' @examples
#'
#' print("DEPRECATED")
#'
#'
#' @export
#'
join_transcripts <- function(.data,
                          transcripts = NULL,
                          all = FALSE,
                          exclude_zeros = FALSE,
                          shape = "wide", ...) {
  UseMethod("join_transcripts", .data)
}
#' @export
join_transcripts.default <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "wide", ...)
  {
    print("tidyseurat says: This function cannot be applied to this object")
  }
#' @export
join_transcripts.Seurat <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "wide", ...)
  {
    
    deprecate_warn("0.2.1", "join_transcripts()", "tidyseurat::join_features()")
    
    
    .data %>%
      join_features(features = transcripts,
                       all = all,
                       exclude_zeros = exclude_zeros,
                       shape = shape, ...)
    
  }