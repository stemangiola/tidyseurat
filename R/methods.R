
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
#' @importFrom ttservice join_features
#'
#' @name join_features
#' @rdname join_features
#'
#' @param .data A Seurat object
#' @param features A vector of feature identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from and gene prefix, for shape="wide"
#'
#' @details This function extracts information for specified features and returns the information in either long or wide format.
#'
#' @return An object containing the information.for the specified features
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
NULL

#' join_features
#'
#' @docType methods
#' @rdname join_features
#'
#' @return An object containing the information.for the specified features
#'
setMethod("join_features", "Seurat",  function(.data,
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
        by = c_(.data)$name
      ) %>%
        select(!!c_(.data)$symbol, .feature, contains(".abundance"), everything()),
      
      # Shape if wide
      ~ (.) %>% left_join(
        get_abundance_sc_wide(
          .data = .data,
          features = features,
          all = all, ...
        ),
        by = c_(.data)$name
      ) 
    )
  
  
  
})


#' Aggregate cells
#'
#' @description Combine cells into groups based on shared variables and aggregate feature counts.
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang enquo
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom purrr map_int
#' @importFrom ttservice aggregate_cells
#' 
#' @name aggregate_cells
#' @rdname aggregate_cells
#' 
#' @param .data A tidySingleCellExperiment object
#' @param .sample A vector of variables by which cells are aggregated
#' @param slot The slot to which the function is applied
#' @param assays The assay to which the function is applied
#' @param aggregation_function The method of cell-feature value aggregation
#' 
#' @return A tibble object
#' 
#' @examples 
#' data("pbmc_small")
#' pbmc_small |>
#'   aggregate_cells(c(groups, letter.idents), assays = "RNA")
#'
#' @export
#'
NULL

#' aggregate_cells
#'
#' @docType methods
#' @rdname aggregate_cells
#'
#' @return An object containing the information.for the specified features
#'
setMethod("aggregate_cells", "Seurat",  function(.data,
                                                 .sample = NULL, 
                                                 slot = "data",
                                                 assays = NULL, 
                                                 aggregation_function = Matrix::rowSums){
  
    .sample = enquo(.sample)

    # Subset only wanted assays
    if(!is.null(assays)){
      DefaultAssay(.data) = assays[1]
      .data@assays = .data@assays[assays]
    }

    .data %>%
      
      tidyseurat::nest(data = -!!.sample) %>%
      mutate(.aggregated_cells = map_int(data, ~ ncol(.x))) %>% 
      mutate(data = map(data, ~ 
                          
                          # loop over assays
                          map2(
                            .x@assays, names(.x@assays),
                            
                            # Get counts
                            ~ GetAssayData(.x, slot = slot) %>%
                              aggregation_function(na.rm = T) %>%
                              tibble::enframe(
                                name  = ".feature",
                                value = sprintf("%s", .y)
                              ) %>%
                              mutate(.feature = as.character(.feature)) 
                          ) %>%
                          Reduce(function(...) full_join(..., by=c(".feature")), .)
                        
      )) %>%
      left_join(.data %>% tidyseurat::as_tibble() %>% subset_tidyseurat(!!.sample)) %>%
      tidyseurat::unnest(data) %>%
      
      tidyr::unite(".sample", !!.sample,  sep="___", remove = FALSE) |> 
      select(.feature, .sample, names(.data@assays), everything()) |> 
      drop_class("tidyseurat_nested") 
    
    
  })




