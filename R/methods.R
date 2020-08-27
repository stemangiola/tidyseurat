

setClass("tidyseurat", contains="Seurat")

#' @importFrom methods show
#' @import Seurat
#' @importFrom magrittr %>%
setMethod(
  f = "show",
  signature = "tidyseurat",
  definition = function(object) {

      object %>%
      as_tibble() %>%
      print()
    

  }
)

#' tidy for seurat
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



#' Add differential transcription information to a tbl using edgeR.
#'
#' \lifecycle{experimental}
#'
#' @description join_transcripts() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_transcripts
#' @rdname join_transcripts
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param transcripts A formula with no response variable, representing the desired linear model
#' @param all The name of the sample column
#' @param exclude_zeros The name of the transcript/gene column
#' @param shape The name of the transcript/gene abundance column
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'\donttest{
#'
#'
#'
#' 	join_transcripts(
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'}
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


