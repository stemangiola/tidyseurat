#' @importFrom tibble as_tibble
to_tib = function(.data){ .data@meta.data %>% as_tibble(rownames = "cell") }



