setClass("tidyseurat", contains="Seurat")

#' @importFrom methods show
#' @import Seurat
#' @importFrom magrittr %>%
setMethod(
  f = "show",
  signature = "tidyseurat",
  definition = function(object) {

      object %>%
      print()
    

  }
)

#' tidy for seurat
#' 
#' @name tidy
#' 
#' @param object A Seurat object
#' 
#' @return A tidyseurat object
#' 
#' @export
tidy <- function(object) {  UseMethod("tidy", object) }

#' @importFrom methods as
#' 
#' @param object A Seurat object
#' 
#' @export
tidy.Seurat <- function(object){  as(object, "tidyseurat") }



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
#'
#' @details This function extracts information for specified transcripts and returns the information in either long or wide format.
#'
#' @return A `tbl` containing the information.for the specified transcripts
#' 
#' @examples
#'
#' pbmc_small %>% 
#' tidy %>% 
#' join_transcripts(transcripts = c("HLA-DRA", "LYZ"))
#'
#'
#' @export
#'
join_transcripts <- function(.data,
                              transcripts = NULL,
                              all = FALSE,
                              exclude_zeros = FALSE,
                              shape = "long") {
  UseMethod("join_transcripts", .data)
}
#' @export
join_transcripts.default <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long")
  {
    print("This function cannot be applied to this object")
  }
#' @export
join_transcripts.tidyseurat <-
  function(.data,
           transcripts = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long")
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
            all = all
          ),
          by = "cell"
        ) 
      )
      

   
  }


