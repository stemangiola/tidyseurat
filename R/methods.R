
#setClass("Seurat")
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
