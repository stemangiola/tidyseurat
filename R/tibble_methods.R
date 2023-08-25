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

    x[[]] %>%
        tibble::as_tibble(rownames=c_(x)$name) %>%

        # Attach reduced dimensions
        when(
            # Only if I have reduced dimensions and special datasets
            length(x@reductions) > 0 ~ (.) %>%
                left_join(
                    get_special_datasets(x, ...) %>%
                    map(~ .x %>% when(
                        # If row == 1 do a trick
                        dim(.) %>% is.null ~ {
                            (.) %>% tibble::enframe() %>%
                                spread(name, value) %>%
                                mutate(!!c_(x)$symbol := rownames(x[[]]))
                        },

                        # Otherwise continue normally
                        ~ as_tibble(., rownames=c_(x)$name)
                    )) %>%
                    reduce(left_join, by=c_(x)$name),
                by=c_(x)$name
                ),
            # Otherwise skip
            ~ (.)
        )
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
