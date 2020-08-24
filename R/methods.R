

setClass("tidyseurat", contains="Seurat")

#' @importFrom methods show
#' @import Seurat
#' @importFrom magrittr %>%
setMethod(
  f = "show",
  signature = "tidyseurat",
  definition = function(object) {

      object %>%
      to_tib %>%
      print()
    

  }
)

#' tidy for seurat
#' @export
tidy <- function(object) {  UseMethod("tidy", object) }

#' @export
tidy.Seurat <- function(object){  as(object, "tidyseurat") }

as_meta_data = function(.data, seurat_object){
  
  col_to_exclude =  get_special_columns(seurat_object)
  
  .data %>% 
    select(-one_of(col_to_exclude)) %>%
    data.frame(row.names = "cell")
}

get_special_columns = function(seurat_object){
  get_special_datasets(seurat_object) %>%
    map_chr(~ .x@cell.embeddings %>% colnames  )
}

get_special_datasets = function(seurat_object){
  seurat_object@reductions 
}

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
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
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
                              all = F,
                              exclude_zeros = F,
                              shape = "long") {
  UseMethod("join_transcripts", .data)
}
#' @export
join_transcripts.default <-
  function(.data,
           transcripts = NULL,
           all = F,
           exclude_zeros = F,
           shape = "long")
  {
    print("This function cannot be applied to this object")
  }
#' @export
join_transcripts.tidyseurat <-
  function(.data,
           transcripts = NULL,
           all = F,
           exclude_zeros = F)
  {
    
    message("tidyseurat says: A data frame is returned for independent data analysis.")

    .data %>%
      as_tibble() %>%
      left_join(
        get_abundance_sc_long(
          .data = .data,
          transcripts = transcripts,
          all = all,
          exclude_zeros = exclude_zeros
        ),
        by = "cell"
      ) %>%
      select(cell, transcript, contains("abundance"), everything())

   
  }


