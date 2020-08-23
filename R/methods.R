
#setClass("Seurat")
setClass("tidyseurat", contains="Seurat")

#' @importFrom methods show
#' @import Seurat
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
tidy <- function(object){  as(object, "tidyseurat") }


