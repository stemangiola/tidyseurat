
setClass("tidyseurat", contains="Seurat")

#' tidy for seurat
#' 
#' @name tidy
#' 
#' @param object A Seurat object
#' 
#' @description 
#' 
#' DEPRECATED. Not needed any more.
#' 
#' @return A tidyseurat object
#' 
#' @export
tidy <- function(object) {  UseMethod("tidy", object) }

#' @importFrom lifecycle deprecate_warn
#' 
#' @param object A Seurat object
#' 
#' @export
tidy.Seurat <- function(object){ 
  
  # DEPRECATE
  deprecate_warn(
    when = "0.2.0",
    what = "tidy()",
    details = "tidyseurat says: tidy() is not needed anymore."
  )
  
  object
  
}

#' @importFrom methods getMethod
setMethod(
  f = "show",
  signature = "Seurat",
  definition = function(object) {
    if (isTRUE(x = getOption(x = "restore_Seurat_show", default = FALSE))) {
      f <- getMethod(
        f = "show",
        signature = "Seurat",
        where = asNamespace(ns = "SeuratObject")
      )
      f(object = object)
    } else {
      object %>%
        print()
    }
  }
)

#' Extract and join information for transcripts.
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
#' @details This function extracts information for specified transcripts and returns the information in either long or wide format.
#'
#' @return A `tbl` containing the information.for the specified transcripts
#' 
#' @examples
#'
#' pbmc_small %>% 
#'  
#' join_transcripts(transcripts = c("HLA-DRA", "LYZ"))
#'
#'
#' @export
#'
join_transcripts <- function(.data,
                              transcripts = NULL,
                              all = FALSE,
                              exclude_zeros = FALSE,
                              shape = "long", ...) {
  UseMethod("join_transcripts", .data)
}
#' @export
join_transcripts.default <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long", ...)
  {
    print("This function cannot be applied to this object")
  }
#' @export
join_transcripts.Seurat <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long", ...)
  {
    
    message("tidyseurat says: A data frame is returned for independent data analysis.")

    .data %>%
      as_tibble() %>%
      
      when(
        
      # Shape is long
      shape == "long" ~ (.) %>% left_join(
        get_abundance_sc_long(
          .data = .data,
          transcripts = transcripts,
          all = all,
          exclude_zeros = exclude_zeros
        ),
        by = "cell"
      ) %>%
        select(cell, transcript, contains("abundance"), everything()),
      
      # Shape if wide
      ~ (.) %>% left_join(
          get_abundance_sc_wide(
            .data = .data,
            transcripts = transcripts,
            all = all, ...
          ),
          by = "cell"
        ) 
      )
      

   
  }


