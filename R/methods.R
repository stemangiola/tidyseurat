
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

#' Extract and join information for features.
#'
#'
#' @description join_features() extracts and joins information for specified features
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_features
#' @rdname join_features
#'
#' @param .data A tidyseurat object
#' @param features A vector of feature identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from
#'
#' @details This function extracts information for specified features and returns the information in either long or wide format.
#'
#' @return A `tbl` containing the information.for the specified features
#' 
#' @examples
#'
#' data("pbmc_small")
#' pbmc_small %>% 
#' join_features(features = c("HLA-DRA", "LYZ"))
#'
#'
#' @export
#'
join_features <- function(.data,
                              features = NULL,
                              all = FALSE,
                              exclude_zeros = FALSE,
                              shape = "long", ...) {
  UseMethod("join_features", .data)
}
#' @export
join_features.default <-
  function(.data,
           features = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long", ...)
  {
    print("tidyseurat says: This function cannot be applied to this object")
  }
#' @export
join_features.Seurat <-
  function(.data,
           features = NULL,
           all = FALSE,
           exclude_zeros = FALSE,
           shape = "long", ...)
  {
    .data %>%
      
      
      when(
        
      # Shape is long
      shape == "long" ~ (.) %>% left_join(
        get_abundance_sc_long(
          .data = .data,
          features = features,
          all = all,
          exclude_zeros = exclude_zeros
        ),
        by = "cell"
      ) %>%
        select(cell, feature, contains("abundance"), everything()),
      
      # Shape if wide
      ~ (.) %>% left_join(
          get_abundance_sc_wide(
            .data = .data,
            features = features,
            all = all, ...
          ),
          by = "cell"
        ) 
      )
      

   
  }


