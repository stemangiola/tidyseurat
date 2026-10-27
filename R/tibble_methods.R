#' @name as_tibble
#' @rdname as_tibble
#' @inherit tibble::as_tibble
#' @return `tibble`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> as_tibble()
#' 
#' @importFrom tibble as_tibble
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom tidyr spread
#' @importFrom tibble enframe
#' @export
as_tibble.Seurat <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=NULL){

    out <- x[[]] %>%
        tibble::as_tibble(rownames=c_(x)$name)

    # Attach reduced dimensions only if present
    if (length(x@reductions) > 0) {
        special_datasets <- get_special_datasets(x, ...) %>%
            map(~ {
                # If row == 1 do a trick
                if (is.null(dim(.x))) {
                    .x %>% tibble::enframe() %>%
                        spread(name, value) %>%
                        mutate(!!c_(x)$symbol := rownames(x[[]]))
                } else {
                    as_tibble(.x, rownames=c_(x)$name)
                }
            }) %>%
            reduce(left_join, by=c_(x)$name)

        out <- out %>% left_join(special_datasets, by=c_(x)$name)
    }

    out
}

#' @name glimpse
#' @rdname glimpse
#' @inherit pillar::glimpse
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> glimpse()
#' 
#' @importFrom tibble glimpse
#' @export
glimpse.tidyseurat <- function(x, width=NULL, ...){
    x %>%
        as_tibble() %>%
        tibble::glimpse(width=width, ...)
}
