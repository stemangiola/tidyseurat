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
#' @return A tidyseurat objector a tibble depending on input
#' 
#' @examples
#'
#' library(dplyr)
#' pbmc_small %>% tidy %>% nest(data = -groups) %>% unnest(data) 
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
#' @return A tidyseurat objector a tibble depending on input
#' 
#' @examples
#'
#' library(dplyr)
#' pbmc_small %>% tidy %>% nest(data = -groups) %>% unnest(data) 
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
#' @importFrom rlang :=
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

#' Extract a character column into multiple columns using regular
#' expression groups
#'
#' Given a regular expression with capturing groups, `extract()` turns
#' each group into a new column. If the groups don't match, or the input
#' is NA, the output will be NA.
#'
#' @inheritParams expand
#' @param data A tidyseurat object
#' @param col Column name or position. This is passed to
#'   [tidyselect::vars_pull()].
#'
#'   This argument is passed by expression and supports
#'   [quasiquotation][rlang::quasiquotation] (you can unquote column
#'   names or column positions).
#' @param into Names of new variables to create as character vector.
#'    Use `NA` to omit the variable in the output.
#' @param regex a regular expression used to extract the desired values.
#'   There should be one group (defined by `()`) for each element of `into`.
#' @param remove If `TRUE`, remove input column from output data frame.
#' @param convert If `TRUE`, will run [type.convert()] with
#'   `as.is = TRUE` on new columns. This is useful if the component
#'   columns are integer, numeric or logical.
#'
#'   NB: this will cause string `"NA"`s to be converted to `NA`s.
#' @param ... Additional arguments passed on to methods.
#' @seealso [separate()] to split up by a separator.
#' @export
#' @examples
#' 
#' pbmc_small %>% tidy %>% extract(groups, into = "g", regex = "g([0-9])", convert = TRUE)
#' 
#' @return A tidyseurat objector a tibble depending on input
#'  
#' @importFrom tidyr extract
#' 
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
	  data %>%
	  as_tibble() %>%
	  tidyr::extract(col = !!col, into = into, regex = regex, remove = remove, convert = convert, ...) %>%
	  as_meta_data(data)
	
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
#' 
#' @return A tidyseurat objector a tibble depending on input
#' 
#' @export
#' @examples
#' # See vignette("pivot") for examples and explanation
#' 
#' library(dplyr)
#' pbmc_small %>% tidy %>% pivot_longer(c(orig.ident, groups), names_to = "name", values_to = "value") 
#' 
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

#' Unite multiple columns into one by pasting strings together
#'
#' Convenience function to paste together multiple columns into one.
#'
#' @param data A data frame.
#' @param col The name of the new column, as a string or symbol.
#'
#'   This argument is passed by expression and supports
#'   [quasiquotation][rlang::quasiquotation] (you can unquote strings
#'   and symbols). The name is captured from the expression with
#'   [rlang::ensym()] (note that this kind of interface where
#'   symbols do not represent actual objects is now discouraged in the
#'   tidyverse; we support it here for backward compatibility).
#' @param ... <[`tidy-select`][tidyr_tidy_select]> Columns to unite
#' @param sep Separator to use between values.
#' @param na.rm If `TRUE`, missing values will be remove prior to uniting
#'   each value.
#' @param remove If `TRUE`, remove input columns from output data frame.
#' @seealso [separate()], the complement.
#' 
#' @return A tidyseurat objector a tibble depending on input
#' 
#' @export
#' @examples
#' 
#' pbmc_small %>% tidy %>% unite("new_col", c(orig.ident, groups)) 
#' 
#' 
unite <- function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE) {
  ellipsis::check_dots_unnamed()
  UseMethod("unite")
}
#' @export
unite.default <-  function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)
{
  cols = enquo(col) 
  tidyr::unite(data, !!cols, ..., sep = sep, remove = remove, na.rm = na.rm)
}

#' @export
unite.tidyseurat <- function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)
{
  
  # Check that we are not modifying a key column
  cols = enquo(col) 
  if(intersect(cols %>% quo_names, get_special_columns(data) %>% c(get_needed_columns())) %>% length %>% gt(0) & remove)
    stop(sprintf("tidyseurat says: you are trying to rename a column that is view only %s (it is not present in the meta.data). If you want to mutate a view-only column, make a copy and mutate that one.", get_special_columns(data) %>% c(get_needed_columns()) %>% paste(collapse=", ")))
  
  
  data@meta.data = data %>% as_tibble() %>% tidyr::unite(!!cols, ..., sep = sep, remove = remove, na.rm = na.rm) %>% as_meta_data(data)
  
  data
  
  
}

#' Separate a character column into multiple columns with a regular
#' expression or numeric locations
#'
#' Given either a regular expression or a vector of character positions,
#' `separate()` turns a single character column into multiple columns.
#'
#' @inheritParams extract
#' @param sep Separator between columns.
#'
#'   If character, `sep` is interpreted as a regular expression. The default
#'   value is a regular expression that matches any sequence of
#'   non-alphanumeric values.
#'
#'   If numeric, `sep` is interpreted as character positions to split at. Positive
#'   values start at 1 at the far-left of the string; negative value start at -1 at
#'   the far-right of the string. The length of `sep` should be one less than
#'   `into`.
#' @param extra If `sep` is a character vector, this controls what
#'   happens when there are too many pieces. There are three valid options:
#'
#'   * "warn" (the default): emit a warning and drop extra values.
#'   * "drop": drop any extra values without a warning.
#'   * "merge": only splits at most `length(into)` times
#' @param fill If `sep` is a character vector, this controls what
#'   happens when there are not enough pieces. There are three valid options:
#'
#'   * "warn" (the default): emit a warning and fill from the right
#'   * "right": fill with missing values on the right
#'   * "left": fill with missing values on the left
#' @seealso [unite()], the complement, [extract()] which uses regular
#'   expression capturing groups.
#' 
#' @return A tidyseurat objector a tibble depending on input
#' 
#' @export
#' @examples
#' 
#'  un = pbmc_small %>% tidy %>% unite("new_col", c(orig.ident, groups)) 
#'  un %>% separate(col = new_col, into= c("orig.ident", "groups"))
#' 
separate <- function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                     convert = FALSE, extra = "warn", fill = "warn", ...) {
  ellipsis::check_dots_used()
  UseMethod("separate")
}
#' @export
separate.default <-  function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                              convert = FALSE, extra = "warn", fill = "warn", ...)
{
  cols = enquo(col)
  tidyr::separate(data, !!cols, into = into, sep = sep, remove = remove,
                  convert = convert, extra = extra, fill = fill, ...)
}

#' @export
separate.tidyseurat <- function(data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
                                convert = FALSE, extra = "warn", fill = "warn", ...)
{
  
  # Check that we are not modifying a key column
  cols = enquo(col)
  if(intersect(cols  %>% quo_names, get_special_columns(data) %>% c(get_needed_columns())) %>% length %>% gt(0) & remove)
    stop(sprintf("tidyseurat says: you are trying to rename a column that is view only %s (it is not present in the meta.data). If you want to mutate a view-only column, make a copy and mutate that one.", get_special_columns(data) %>% c(get_needed_columns()) %>% paste(collapse=", ")))
  
   
  data@meta.data =
    data %>%
    as_tibble() %>% 
    tidyr::separate(!!cols, into = into, sep = sep, remove = remove,  convert = convert, extra = extra, fill = fill, ...) %>%
    as_meta_data(data)
  
  
  data
  
  
}