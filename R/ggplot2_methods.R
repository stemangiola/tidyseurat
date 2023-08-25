#' @name ggplot
#' @rdname ggplot
#' @inherit ggplot2::ggplot
#' @title Create a new \code{ggplot} from a \code{tidyseurat}
#' @return `ggplot`
#'
#' @examples
#' library(ggplot2)
#' data(pbmc_small)
#' pbmc_small |> 
#'   ggplot(aes(groups, nCount_RNA)) +
#'   geom_boxplot()
#' 
#' @importFrom purrr map
#' @importFrom rlang quo_name
#' @importFrom ggplot2 aes ggplot
#' @export
ggplot.Seurat <- function(data=NULL, mapping=aes(),
    ..., environment=parent.frame()) {
  
    # Deprecation of special column names
    .cols <- mapping %>% 
        unlist() %>% map(~ quo_name(.x)) %>% 
        unlist() %>% as.character()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
  
    data %>%
        as_tibble() %>%
        ggplot2::ggplot(mapping=mapping)
}
