#' @importFrom methods getMethod
setMethod(
    f="show",
    signature="Seurat",
    definition=function(object) {
        if (isTRUE(x=getOption(x="restore_Seurat_show", default=FALSE))) {
            f <- getMethod(
                f="show",
                signature="Seurat",
                where=asNamespace(ns="SeuratObject"))
            f(object=object)
        } else { print(object) }
    }
)

setClass("tidyseurat", contains="Seurat")

#' @name tidy
#' @rdname tidy
#' @title tidy for `Seurat`
#'
#' @param object A `Seurat` object.
#' @return A `tidyseurat` object.
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small
#'
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @rdname tidy
#' @importFrom lifecycle deprecate_warn
#' @export
tidy.Seurat <- function(object){ 
  
    # DEPRECATE
    deprecate_warn(
        when="0.2.0",
        what="tidy()",
        details="tidyseurat says: tidy() is not needed anymore."
    )
  
    return(object)
}



#' @name join_features
#' @rdname join_features
#' @inherit ttservice::join_features
#' @aliases join_features,Seurat-method
#'
#' @param .data A tidyseurat object
#' @param assay assay name to extract feature abundance
#' @param slot slot name to extract feature abundance
#' 
#' @return A `tidyseurat` object
#'   containing information for the specified features.
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small %>% join_features(
#'   features=c("HLA-DRA", "LYZ"))
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr contains
#' @importFrom dplyr everything
#' @importFrom ttservice join_features
#' @export
setMethod("join_features", "Seurat", function(.data,
    features=NULL, all=FALSE, exclude_zeros=FALSE, shape="long",
    assay=NULL, slot="data", ...) {
  
  .feature = NULL
  
  if(shape == "long")
    .data |> 
    left_join(
      get_abundance_sc_long(
        .data=.data,
        features=features,
        all=all,
        exclude_zeros=exclude_zeros,
        assay=assay,
        slot=slot,
        ...
      ),
      by=c_(.data)$name
    ) %>%
    select(!!c_(.data)$symbol, .feature,
           contains(".abundance"), everything())
  else
    .data |> 
    left_join(
      get_abundance_sc_wide(
        .data=.data,
        features=features,
        all=all,
        assay=assay,
        slot=slot,
        ...
      ),
      by=c_(.data)$name
    ) 

})


#' @name aggregate_cells
#' @rdname aggregate_cells
#' @inherit ttservice::aggregate_cells
#' @aliases aggregate_cells,Seurat-method
#' 
#' @param .data A tidyseurat object
#' 
#' @examples 
#' data(pbmc_small)
#' pbmc_small_pseudo_bulk <- pbmc_small |>
#'   aggregate_cells(c(groups, letter.idents), assays="RNA")
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom ttservice aggregate_cells
#' @importFrom SeuratObject DefaultAssay
#' @importFrom Seurat DietSeurat
#' @importFrom purrr map_int
#' @export
setMethod("aggregate_cells", "Seurat",  function(.data,
    .sample=NULL, slot="data", assays=NULL,
    aggregation_function=Matrix::rowSums, ...){
    # Solve NOTE  
    data <- NULL
    .feature <- NULL
  
    .sample <- enquo(.sample)

    # Subset only wanted assays
    if(!is.null(assays)){
        DefaultAssay(.data) <- assays[1]
        .data = .data |> DietSeurat(assays = assays)
    }

    .data %>%
        nest(data=-!!.sample) %>%
        mutate(.aggregated_cells=map_int(data, ~ ncol(.x))) %>% 
        mutate(
            data=map(data, ~ 
                # Loop over assays
                map2(.x@assays, names(.x@assays),
                    # Get counts
                    ~ GetAssayData(.x, layer=slot) %>%
                        aggregation_function(na.rm=T) %>%
                        tibble::enframe(
                            name=".feature",
                            value=sprintf("%s", .y)
                        ) %>%
                        mutate(.feature=as.character(.feature)) 
                ) %>%
                Reduce(function(...) full_join(..., by=c(".feature")), .), 
                .progress = TRUE          
        )) %>%
        left_join(
            .data %>%
                as_tibble() %>%
                subset_tidyseurat(!!.sample)) %>%
        unnest(data) %>%
        tidyr::unite(".sample", !!.sample,  sep="___", remove=FALSE) |> 
        select(.feature, .sample, names(.data@assays), everything()) |> 
        drop_class("tidyseurat_nested") 
})