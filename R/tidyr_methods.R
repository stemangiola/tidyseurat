#' unnest
#'
#' @importFrom tidyr unnest
#'
#' @param .data A tbl. (See tidyr)
#' @param cols <[`tidy-select`][tidyr_tidy_select]> Columns to unnest.
#'   If you `unnest()` multiple columns, parallel entries must be of
#'   compatible sizes, i.e. they're either equal or length 1 (following the
#'   standard tidyverse recycling rules).
#' @param ... <[`tidy-select`][tidyr_tidy_select]> Columns to nest, specified
#'   using name-variable pairs of the form `new_col = c(col1, col2, col3)`.
#'   The right hand side can be any valid tidy select expression.
#'
#'   \Sexpr[results=rd, stage=render]{lifecycle::badge("deprecated")}:
#'   previously you could write `df %>% nest(x, y, z)` and `df %>%
#'   unnest(x, y, z)`. Convert to `df %>% nest(data = c(x, y, z))`.
#'   and `df %>% unnest(c(x, y, z))`.
#'
#'   If you previously created new variable in `unnest()` you'll now need to
#'   do it explicitly with `mutate()`. Convert `df %>% unnest(y = fun(x, y, z))`
#'   to `df %>% mutate(y = fun(x, y, z)) %>% unnest(y)`.
#' @param names_sep If `NULL`, the default, the names will be left
#'   as is. In `nest()`, inner names will come from the former outer names;
#'   in `unnest()`, the new outer names will come from the inner names.
#'
#'   If a string, the inner and outer names will be used together. In `nest()`,
#'   the names of the new outer columns will be formed by pasting together the
#'   outer and the inner column names, separated by `names_sep`. In `unnest()`,
#'   the new inner names will have the outer names (+ `names_sep`) automatically
#'   stripped. This makes `names_sep` roughly symmetric between nesting and unnesting.
#' @param keep_empty See tidyr::unnest
#' @param names_repair See tidyr::unnest
#' @param ptype See tidyr::unnest
#' 
#'
#' @return A tt object
#'
#' @examples
#'
#' library(dplyr)
#' 
#' nest(tidyseurat(tidyseurat::counts_mini, sample, transcript, count), data = -transcript) %>%
#' unnest(data)
#'
#' @rdname tidyr-methods
#'
#' @export
unnest <- function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
                    names_sep = NULL, names_repair = "check_unique")  {
  UseMethod("unnest")
}

#' @export
#' @rdname tidyr-methods
unnest.default <-  function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
                             names_sep = NULL, names_repair = "check_unique")
{
  cols <- enquo(cols)
  tidyr::unnest(.data, !!cols, ..., keep_empty = keep_empty, ptype = ptype, 
                names_sep = names_sep, names_repair = names_repair)
}

#' @importFrom rlang quo_name
#' 
#' @export
#' @rdname tidyr-methods
unnest.tidyseurat_nested <- function (.data, cols, ..., keep_empty = FALSE, ptype = NULL, 
                                    names_sep = NULL, names_repair = "check_unique")
{
  # Need this otherwise crashes map
  .data_ = .data
  
  cols <- enquo(cols)

  .data_ %>% 
    when(
      
      # If my only column to unnest is tidyseurat
      pull(., !!cols) %>% .[[1]] %>% class %>% as.character() %>% eq("tidyseurat") ~  
        
        # Do my trick to unnest
        mutate(., !!cols := imap(
          !!cols, ~ .x %>%
            bind_cols(
              .data_ %>% select(-!!cols) %>% slice(.y, each= ncol(.x))
            )
        )) %>%
        pull(!!cols) %>%
        reduce(bind_rows),
      
      # Else do normal stuff
      ~ tidyr::unnest(., !!cols, ..., keep_empty = keep_empty, ptype = ptype, names_sep = names_sep, names_repair = names_repair)
      
    )

}

#' nest
#'
#' @importFrom tidyr nest
#'
#' @param .data A tbl. (See tidyr)
#' @param ... Name-variable pairs of the form new_col = c(col1, col2, col3) (See tidyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' nest(tidyseurat(tidyseurat::counts_mini, sample, transcript, count), data = -transcript)
#'
#' @rdname tidyr-methods
#'
#' @export
nest <- function (.data, ...)  {
  UseMethod("nest")
}

#' @export
#' @rdname tidyr-methods
nest.default <-  function (.data, ...)
{
  tidyr::nest(.data, ...)
}

#' @export
#' @rdname tidyr-methods
nest.tidyseurat <- function (.data, ...)
{
  my_data__ = .data
  cols <- enquos(...)
  col_name_data  = names(cols)
  
  my_data__ %>%
    
    # This is needed otherwise nest goes into loop and fails
    to_tib %>%
    tidyr::nest(...) %>%
    
    mutate(
      !!as.symbol(col_name_data) := map(
        !!as.symbol(col_name_data),
        ~ my_data__ %>% 
          
          # Subset cells
          filter(cell %in% .x$cell) %>%
          
          # Subset columns
          select(colnames(.x))
      )) %>%
    
    # Coerce to tidyseurat_nested for unnesting
    add_class("tidyseurat_nested")
  
}

#' @export
extract <- function  (data, col, into, regex = "([[:alnum:]]+)", remove = TRUE, 
										 convert = FALSE, ...)   {
	UseMethod("extract")
}

#' @export
extract.default <-  function  (data, col, into, regex = "([[:alnum:]]+)", remove = TRUE, 
															convert = FALSE, ...) 
{
	col = enquo(col)
	tidyr::extract(col = !!col, into = into, regex = regex, remove = remove, 
								 convert = convert, ...) 
}

#' @export
extract.tidyseurat <- function  (data, col, into, regex = "([[:alnum:]]+)", remove = TRUE, 
														convert = FALSE, ...) 
{
	
	col = enquo(col)
	
	data@meta.data = 
	  data@meta.data %>%
		tidyr::extract(col = !!col, into = into, regex = regex, remove = remove, 
									 convert = convert, ...)  
	
	data
		
	
}