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
#' @importFrom purrr imap
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
      pull(., !!cols) %>% .[[1]] %>% class %>% as.character() %>% eq("tidyseurat") %>% any ~  
        
        # Do my trick to unnest
        mutate(., !!cols := imap(
          !!cols, ~ .x %>%
            bind_cols(
              
              # Attach back the columns used for nesting
              .data_ %>% select(-!!cols) %>% slice(rep(.y, ncol(.x)))
            )
        )) %>%
        pull(!!cols) %>%
        reduce(bind_rows),
      
      # Else do normal stuff
      ~ (.) %>% 
        drop_class("tidyseurat_nested") %>%
        tidyr::unnest( !!cols, ..., keep_empty = keep_empty, ptype = ptype, names_sep = names_sep, names_repair = names_repair) %>%
        add_class("tidyseurat_nested")
      
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

#' @importFrom rlang enquos
#' 
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

#' Pivot data from wide to long
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("maturing")}
#'
#' `pivot_longer()` "lengthens" data, increasing the number of rows and
#' decreasing the number of columns. The inverse transformation is
#' [pivot_wider()]
#'
#' Learn more in `vignette("pivot")`.
#'
#' @details
#' `pivot_longer()` is an updated approach to [gather()], designed to be both
#' simpler to use and to handle more use cases. We recommend you use
#' `pivot_longer()` for new code; `gather()` isn't going away but is no longer
#' under active development.
#'
#' @param data A data frame to pivot.
#' @param cols <[`tidy-select`][tidyr_tidy_select]> Columns to pivot into
#'   longer format.
#' @param names_to A string specifying the name of the column to create
#'   from the data stored in the column names of `data`.
#'
#'   Can be a character vector, creating multiple columns, if `names_sep`
#'   or `names_pattern` is provided. In this case, there are two special
#'   values you can take advantage of:
#'
#'   * `NA` will discard that component of the name.
#'   * `.value` indicates that component of the name defines the name of the
#'     column containing the cell values, overriding `values_to`.
#' @param names_prefix A regular expression used to remove matching text
#'   from the start of each variable name.
#' @param names_sep,names_pattern If `names_to` contains multiple values,
#'   these arguments control how the column name is broken up.
#'
#'   `names_sep` takes the same specification as [separate()], and can either
#'   be a numeric vector (specifying positions to break on), or a single string
#'   (specifying a regular expression to split on).
#'
#'   `names_pattern` takes the same specification as [extract()], a regular
#'   expression containing matching groups (`()`).
#'
#'   If these arguments do not give you enough control, use
#'   `pivot_longer_spec()` to create a spec object and process manually as
#'   needed.
#' @param names_repair What happens if the output has invalid column names?
#'   The default, `"check_unique"` is to error if the columns are duplicated.
#'   Use `"minimal"` to allow duplicates in the output, or `"unique"` to
#'   de-duplicated by adding numeric suffixes. See [vctrs::vec_as_names()]
#'   for more options.
#' @param values_to A string specifying the name of the column to create
#'   from the data stored in cell values. If `names_to` is a character
#'   containing the special `.value` sentinel, this value will be ignored,
#'   and the name of the value column will be derived from part of the
#'   existing column names.
#' @param values_drop_na If `TRUE`, will drop rows that contain only `NA`s
#'   in the `value_to` column. This effectively converts explicit missing values
#'   to implicit missing values, and should generally be used only when missing
#'   values in `data` were created by its structure.
#' @param names_transform,values_transform A list of column name-function pairs.
#'   Use these arguments if you need to change the type of specific columns.
#'   For example, `names_transform = list(week = as.integer)` would convert
#'   a character week variable to an integer.
#' @param names_ptypes,values_ptypes A list of column name-prototype pairs.
#'   A prototype (or ptype for short) is a zero-length vector (like `integer()`
#'   or `numeric()`) that defines the type, class, and attributes of a vector.
#'   Use these arguments to confirm that the created columns are the types that
#'   you expect.
#'
#'   If not specified, the type of the columns generated from `names_to` will
#'   be character, and the type of the variables generated from `values_to`
#'   will be the common type of the input columns used to generate them.
#' @param ... Additional arguments passed on to methods.
#' @export
#' @examples
#' # See vignette("pivot") for examples and explanation
pivot_longer <- function(data,
                         cols,
                         names_to = "name",
                         names_prefix = NULL,
                         names_sep = NULL,
                         names_pattern = NULL,
                         names_ptypes = list(),
                         names_transform = list(),
                         names_repair = "check_unique",
                         values_to = "value",
                         values_drop_na = FALSE,
                         values_ptypes = list(),
                         values_transform = list(),
                         ...
) {
  
  ellipsis::check_dots_used()
  UseMethod("pivot_longer")
}

#' @export
pivot_longer.default <- function(data,
                                 cols,
                                 names_to = "name",
                                 names_prefix = NULL,
                                 names_sep = NULL,
                                 names_pattern = NULL,
                                 names_ptypes = list(),
                                 names_transform = list(),
                                 names_repair = "check_unique",
                                 values_to = "value",
                                 values_drop_na = FALSE,
                                 values_ptypes = list(),
                                 values_transform = list(),
                                 ...
) {
  cols <- enquo(cols)
  tidyr::pivot_longer(data,
                      cols,
                      names_to = names_to,
                      names_prefix = names_prefix,
                      names_sep = names_sep,
                      names_pattern = names_pattern,
                      names_ptypes = names_ptypes,
                      names_transform = names_transform,
                      names_repair = names_repair,
                      values_to =values_to,
                      values_drop_na = values_drop_na,
                      values_ptypes = values_ptypes,
                      values_transform = values_transform,
                      ...
  )
  
}

#' @export
pivot_longer.tidyseurat <- function(data,
                                  cols,
                                  names_to = "name",
                                  names_prefix = NULL,
                                  names_sep = NULL,
                                  names_pattern = NULL,
                                  names_ptypes = list(),
                                  names_transform = list(),
                                  names_repair = "check_unique",
                                  values_to = "value",
                                  values_drop_na = FALSE,
                                  values_ptypes = list(),
                                  values_transform = list(),
                                  ...
) {
  cols <- enquo(cols) %>% quo_names()
  
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  data %>%
    as_tibble() %>%
    tidyr::pivot_longer(cols,
                        names_to = names_to,
                        names_prefix = names_prefix,
                        names_sep = names_sep,
                        names_pattern = names_pattern,
                        names_ptypes = names_ptypes,
                        names_transform = names_transform,
                        names_repair = names_repair,
                        values_to =values_to,
                        values_drop_na = values_drop_na,
                        values_ptypes = values_ptypes,
                        values_transform = values_transform,
                        ...
    )
  
}