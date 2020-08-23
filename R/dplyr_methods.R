#' drplyr-methods
#' 
#' @rdname dplyr-methods
#' 
#' @return A tibble


#' Arrange rows by column values
#'
#'
#' @description
#' `arrange()` order the rows of a data frame rows by the values of selected
#' columns.
#'
#' Unlike other dplyr verbs, `arrange()` largely ignores grouping; you
#' need to explicit mention grouping variables (or use  `by_group = TRUE`)
#' in order to group by them, and functions of variables are evaluated
#' once per data frame, not once per group.
#'
#' @details
#' ## Locales
#' The sort order for character vectors will depend on the collating sequence
#' of the locale in use: see [locales()].
#'
#' ## Missing values
#' Unlike base sorting with `sort()`, `NA` are:
#' * always sorted to the end for local data, even when wrapped with `desc()`.
#' * treated differently for remote data, depending on the backend.
#'
#' @return
#' An object of the same type as `.data`.
#'
#' * All rows appear in the output, but (usually) in a different place.
#' * Columns are not modified.
#' * Groups are not modified.
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @export
#' @param .data A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr). See *Methods*, below, for
#'   more details.
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Variables, or functions or
#'   variables. Use [desc()] to sort a variable in descending order.
#' @param .by_group If TRUE, will sort first by grouping variable. Applies to grouped data frames only.
#' 
#' @family single table verbs
#' @examples
#' `%>%` = magrittr::`%>%`
#' arrange(mtcars, cyl, disp)
arrange <- function(.data, ..., .by_group = FALSE) {
  UseMethod("arrange")
}

#' @param .by_group If `TRUE`, will sort first by grouping variable. Applies to
#'   grouped data frames only.
#' @rdname dplyr-methods
#' @export
#'
############# START ADDED tidyseurat ###################################
#' @inheritParams arrange
arrange.default <- function(.data, ..., .by_group = FALSE) {
  
  dplyr::arrange(.data, ..., .by_group = .by_group)
  
}

#' @export
#' @inheritParams arrange
arrange.tidyseurat <- function(.data, ..., .by_group = FALSE) {
  
  .data@meta.data = dplyr::arrange( .data@meta.data, ..., .by_group = .by_group) 
  
  .data
  
}

############# END ADDED tidyseurat #####################################

#' Efficiently bind multiple data frames by row and column
#'
#' This is an efficient implementation of the common pattern of
#' `do.call(rbind, dfs)` or `do.call(cbind, dfs)` for binding many
#' data frames into one.
#'
#' The output of `bind_rows()` will contain a column if that column
#' appears in any of the inputs.
#'
#' @param ... Data frames to combine.
#'
#'   Each argument can either be a data frame, a list that could be a data
#'   frame, or a list of data frames.
#'
#'   When row-binding, columns are matched by name, and any missing
#'   columns will be filled with NA.
#'
#'   When column-binding, rows are matched by position, so all data
#'   frames must have the same number of rows. To match by value, not
#'   position, see [mutate-joins].
#' @param .id Data frame identifier.
#'
#'   When `.id` is supplied, a new column of identifiers is
#'   created to link each row to its original data frame. The labels
#'   are taken from the named arguments to `bind_rows()`. When a
#'   list of data frames is supplied, the labels are taken from the
#'   names of the list. If no names are found a numeric sequence is
#'   used instead.
#' @return `bind_rows()` and `bind_cols()` return the same type as
#'   the first input, either a data frame, `tbl_df`, or `grouped_df`.
#' @examples
#' `%>%` = magrittr::`%>%`
#' one <- mtcars[1:4, ]
#' two <- mtcars[11:14, ]
#'
#' # You can supply data frames as arguments:
#' bind_rows(one, two)
#'
#' @name bind
NULL

############# START ADDED tidyseurat #####################################

#' @rdname dplyr-methods
#' 
#' @inheritParams bind
#' 
#' @export
#'
bind_rows <- function(..., .id = NULL) {
  UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(..., .id = NULL)
{
  dplyr::bind_rows(..., .id = .id)
}

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
#' 
bind_rows.tidyseurat <- function(..., .id = NULL)
{
  

    
  tts = flatten_if(dots_values(...), is_spliced)

  # Check if cell with same name
  merge(  tts[[1]] ,  tts[[2]] ,  add.cell.ids = 1:2 ) %>% tidy
  
  
}

############# END ADDED tidyseurat #####################################


############# START ADDED tidyseurat #####################################
#' @export
#' 
#' @inheritParams bind
#' 
#' @rdname dplyr-methods
bind_cols <- function(..., .id = NULL) {
  UseMethod("bind_cols")
}

#' @export
bind_cols.default <-  function(..., .id = NULL)
{
  dplyr::bind_cols(..., .id = .id)
}

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
#' 
bind_cols.tidyseurat <- function(..., .id = NULL)
{
  
  tts = 	tts = flatten_if(dots_values(...), is_spliced) 
  
  tts[[1]]@meta.data = dplyr::bind_cols( tts[[1]]@meta.data, tts[[2]], .id = .id) 
  
  tts[[1]]
  
}

############# END ADDED tidyseurat #####################################
############# START ADDED tidyseurat #####################################

# #' @importFrom dplyr arrange_all
# #' @export
# dplyr::arrange_all
#
# #' @importFrom dplyr arrange_at
# #' @export
# dplyr::arrange_at
#
# #' @importFrom dplyr arrange_if
# #' @export
# dplyr::arrange_if

############# END ADDED tidyseurat #####################################
############# START ADDED tidyseurat #####################################

#' distinct
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#' @param .keep_all If TRUE, keep all variables in .data. If a combination of ... is not distinct, this keeps the first row of values. (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'
#' distinct(tidyseurat::counts_mini)
#'
#'
#' @export
distinct <- function (.data, ..., .keep_all = FALSE)  {
  UseMethod("distinct")
}

#' @export
distinct.default <-  function (.data, ..., .keep_all = FALSE)
{
  dplyr::distinct(.data, ..., .keep_all = FALSE)
}

#' @export
distinct.tidyseurat <- function (.data, ..., .keep_all = FALSE)
{
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  .data@meta.data %>%
    as_tibble(rownames="cell") %>%
    dplyr::distinct(..., .keep_all = .keep_all)
  
}
############# END ADDED tidyseurat #####################################

############# START ADDED tidyseurat #####################################

# #' @importFrom dplyr distinct_all
# #' @export
# dplyr::distinct_all
#
# #' @importFrom dplyr distinct_at
# #' @export
# dplyr::distinct_at
#
# #' @importFrom dplyr distinct_if
# #' @export
# dplyr::distinct_if

############# END ADDED tidyseurat #####################################

#' Subset rows using column values
#'
#' `filter()` retains the rows where the conditions you provide a `TRUE`. Note
#' that, unlike base subsetting with `[`, rows where the condition evaluates
#' to `NA` are dropped.
#'
#' dplyr is not yet smart enough to optimise filtering optimisation
#' on grouped datasets that don't need grouped calculations. For this reason,
#' filtering is often considerably faster on [ungroup()]ed data.
#'
#' @section Useful filter functions:
#'
#' * [`==`], [`>`], [`>=`] etc
#' * [`&`], [`|`], [`!`], [xor()]
#' * [is.na()]
#' * [between()], [near()]
#'
#' @section Grouped tibbles:
#'
#' Because filtering expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped filtering:
#'
#'
#' The former keeps rows with `mass` greater than the global average
#' whereas the latter keeps rows with `mass` greater than the gender
#'
#' average.
#' @family single table verbs
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Logical predicates defined in
#'   terms of the variables in `.data`.
#'   Multiple conditions are combined with `&`. Only rows where the
#'   condition evaluates to `TRUE` are kept.
#' @param .preserve when `FALSE` (the default), the grouping structure
#'   is recalculated based on the resulting data, otherwise it is kept as is.
#' @return
#' An object of the same type as `.data`.
#'
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The number of groups may be reduced (if `.preserve` is not `TRUE`).
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @seealso [filter_all()], [filter_if()] and [filter_at()].
#' @export
#' @examples
#'
#' # Learn more in ?dplyr_tidy_eval
############# START ADDED tidyseurat #####################################
#' @export
filter <- function (.data, ..., .preserve = FALSE)  {
  UseMethod("filter")
}

#' @export
filter.default <-  function (.data, ..., .preserve = FALSE)
{
  dplyr::filter(.data, ..., .preserve = .preserve)
}

#' @export
filter.tidyseurat <- function (.data, ..., .preserve = FALSE)
{
  new_meta = dplyr::filter(.data@meta.data, ..., .preserve = .preserve)
  new_obj = subset(.data,   cells = rownames(new_meta ))
  new_obj@meta.data = new_meta
  
  new_obj
                   
}
############# END ADDED tidyseurat #####################################


#' Group by one or more variables
#'
#' @description
#' Most data operations are done on groups defined by variables.
#' `group_by()` takes an existing tbl and converts it into a grouped tbl
#' where operations are performed "by group". `ungroup()` removes grouping.
#'
#' @family grouping functions
#' @inheritParams arrange
#' @param ... In `group_by()`, variables or computations to group by.
#'   In `ungroup()`, variables to remove from the grouping.
#' @param .add When `FALSE`, the default, `group_by()` will
#'   override existing groups. To add to the existing groups, use
#'   `.add = TRUE`.
#'
#'   This argument was previously called `add`, but that prevented
#'   creating a new grouping variable called `add`, and conflicts with
#'   our naming conventions.
#' @param .drop When `.drop = TRUE`, empty groups are dropped. See [group_by_drop_default()] for
#'   what the default value is for this argument.
#' @return A [grouped data frame][grouped_df()], unless the combination of `...` and `add`
#'   yields a non empty set of grouping columns, a regular (ungrouped) data frame
#'   otherwise.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' by_cyl <- mtcars %>% group_by(cyl)
#'

############# START ADDED tidyseurat #####################################
#' @export
group_by <- function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))  {
  UseMethod("group_by")
}

#' @export
group_by.default <-  function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))
{
  dplyr::group_by(.data, ...,  .drop = .drop)
}

#' @export
group_by.tidyseurat <- function (.data, ..., .add = FALSE, .drop = group_by_drop_default(.data))
{
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  .data@meta.data %>%
    
    as_tibble(rownames="cell")
    dplyr::group_by( ..., .drop = .drop) 
  
}
############# END ADDED tidyseurat #####################################

#' Summarise each group to fewer rows
#'
#' @description
#' `summarise()` creates a new data frame. It will have one (or more) rows for
#' each combination of grouping variables; if there are no grouping variables,
#' the output will have a single row summarising all observations in the input.
#' It will contain one column for each grouping variable and one column
#' for each of the summary statistics that you have specified.
#'
#' `summarise()` and `summarize()` are synonyms.
#'
#' @section Useful functions:
#'
#' * Center: [mean()], [median()]
#' * Spread: [sd()], [IQR()], [mad()]
#' * Range: [min()], [max()], [quantile()]
#' * Position: [first()], [last()], [nth()],
#' * Count: [n()], [n_distinct()]
#' * Logical: [any()], [all()]
#'
#' @section Backend variations:
#'
#' The data frame backend supports creating a variable and using it in the
#' same summary. This means that previously created summary variables can be
#' further transformed or combined within the summary, as in [mutate()].
#' However, it also means that summary variables with the same names as previous
#' variables overwrite them, making those variables unavailable to later summary
#' variables.
#'
#' This behaviour may not be supported in other backends. To avoid unexpected
#' results, consider using new names for your summary variables, especially when
#' creating multiple summaries.
#'
#' @export
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs of summary
#'   functions. The name will be the name of the variable in the result.
#'
#'   The value can be:
#'
#'   * A vector of length 1, e.g. `min(x)`, `n()`, or `sum(is.na(y))`.
#'   * A vector of length `n`, e.g. `quantile()`.
#'   * A data frame, to add multiple columns from a single expression.
#' @family single table verbs
#' @return
#' An object _usually_ of the same type as `.data`.
#'
#' * The rows come from the underlying `group_keys()`.
#' * The columns are a combination of the grouping keys and the summary
#'   expressions that you provide.
#' * If `x` is grouped by more than one variable, the output will be another
#'   [grouped_df] with the right-most group removed.
#' * If `x` is grouped by one variable, or is not grouped, the output will
#'   be a [tibble].
#' * Data frame attributes are **not** preserved, because `summarise()`
#'   fundamentally creates a new data frame.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @examples
#' `%>%` = magrittr::`%>%`
#' # A summary applied to ungrouped tbl returns a single row
#' mtcars %>%
#'   summarise(mean = mean(disp))
#'
############# START ADDED tidyseurat #####################################
#' @export
summarise <- function (.data, ...)  {
  UseMethod("summarise")
}

#' @export
summarise.default <-  function (.data, ...)
{
  dplyr::summarise(.data, ...)
}

#' @export
summarise.tidyseurat <- function (.data, ...)
{
  
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  .data@meta.data %>%
    
    as_tibble(rownames="cell") %>%
    dplyr::summarise( ...)
  
}
############# END ADDED tidyseurat #####################################

#' Create, modify, and delete columns
#'
#' `mutate()` adds new variables and preserves existing ones;
#' `transmute()` adds new variables and drops existing ones.
#' New variables overwrite existing variables of the same name.
#' Variables can be removed by setting their value to `NULL`.
#'
#' @section Useful mutate functions:
#'
#' * [`+`], [`-`], [log()], etc., for their usual mathematical meanings
#'
#' * [lead()], [lag()]
#'
#' * [dense_rank()], [min_rank()], [percent_rank()], [row_number()],
#'   [cume_dist()], [ntile()]
#'
#' * [cumsum()], [cummean()], [cummin()], [cummax()], [cumany()], [cumall()]
#'
#' * [na_if()], [coalesce()]
#'
#' * [if_else()], [recode()], [case_when()]
#'
#' @section Grouped tibbles:
#'
#' Because mutating expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped mutate:
#'

#' With the grouped equivalent:
#'
#' The former normalises `mass` by the global average whereas the
#' latter normalises by the averages within gender levels.
#'
#' @export
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Name-value pairs.
#'   The name gives the name of the column in the output.
#'
#'   The value can be:
#'
#'   * A vector of length 1, which will be recycled to the correct length.
#'   * A vector the same length as the current group (or the whole data frame
#'     if ungrouped).
#'   * `NULL`, to remove the column.
#'   * A data frame or tibble, to create multiple columns in the output.
#'
#' @family single table verbs
#' @return
#' An object of the same type as `.data`.
#'
#' For `mutate()`:
#'
#' * Rows are not affected.
#' * Existing columns will be preserved unless explicitly modified.
#' * New columns will be added to the right of existing columns.
#' * Columns given value `NULL` will be removed
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#'
#' For `transmute()`:
#'
#' * Rows are not affected.
#' * Apart from grouping variables, existing columns will be remove unless
#'   explicitly kept.
#' * Column order matches order of expressions.
#' * Groups will be recomputed if a grouping variable is mutated.
#' * Data frame attributes are preserved.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' # Newly created variables are available immediately
#' mtcars %>% as_tibble() %>% mutate(
#'   cyl2 = cyl * 2,
#'   cyl4 = cyl2 * 2
#' )
#'
############# START ADDED tidyseurat #####################################
#' @export
mutate <- function(.data, ...) {
  UseMethod("mutate")
}

#' @export
mutate.default <-  function(.data, ...)
{
  dplyr::mutate(.data, ...)
}

#' @export
mutate.tidyseurat <- function(.data, ...)
{

  .data@meta.data = dplyr::mutate(.data@meta.data, ...) 

  .data
}
#' @export
mutate.nested_tidyseurat <- function(.data, ...)
{
  .data %>%
    drop_class(c("nested_tidyseurat", "tt")) %>%
    dplyr::mutate(...) %>%
    
    # Attach attributes
    reattach_internals(.data) %>%
    
    # Add class
    add_class("tt") %>%
    add_class("nested_tidyseurat")
  
  
}
############# END ADDED tidyseurat #####################################

#' Rename columns
#'
#' Rename individual variables using `new_name = old_name` syntax.
#'
#' @section Scoped selection and renaming:
#'
#' Use the three scoped variants ([rename_all()], [rename_if()], [rename_at()])
#' to renaming a set of variables with a function.
#'
#' @inheritParams arrange
#' @param ... <[`tidy-select`][dplyr_tidy_select]> Use `new_name = old_name`
#'   to rename selected variables.
#' @return
#' An object of the same type as `.data`.
#' * Rows are not affected.
#' * Column names are changed; column order is preserved
#' * Data frame attributes are preserved.
#' * Groups are updated to reflect new names.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @family single table verbs
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' iris <- as_tibble(iris) # so it prints a little nicer
#' rename(iris, petal_length = Petal.Length)
############# START ADDED tidyseurat #####################################
#' @export
rename <- function(.data, ...) {
  UseMethod("rename")
}

#' @export
rename.default <-  function(.data, ...)
{
  dplyr::rename(.data, ...)
}

#' @export
rename.tidyseurat <- function(.data, ...)
{
  .data@meta.data = dplyr::rename( .data@meta.data,  ...)
  
  .data
 
  
}
############# END ADDED tidyseurat #####################################

#' Group input by rows
#'
#'
#' See [this repository](https://github.com/jennybc/row-oriented-workflows)
#' for alternative ways to perform row-wise operations.
#'
#' `rowwise()` is used for the results of [do()] when you
#' create list-variables. It is also useful to support arbitrary
#' complex operations that need to be applied to each row.
#'
#' Currently, rowwise grouping only works with data frames. Its
#' main impact is to allow you to work with list-variables in
#' [summarise()] and [mutate()] without having to
#' use \code{[[1]]}. This makes `summarise()` on a rowwise tbl
#' effectively equivalent to [plyr::ldply()].
#'
#' @param .data Input data frame.
#'
#' @return A `tbl`
#'
#'   A `tbl`
#'
#' @export
#' @examples
#' `%>%` = magrittr::`%>%`
#' df <- expand.grid(x = 1:3, y = 3:1)
#' df_done <- df %>% rowwise() %>% do(i = seq(.$x, .$y))
############# START ADDED tidyseurat #####################################
#' @export
rowwise <- function(.data) {
  UseMethod("rowwise")
}

#' @export
rowwise.default <-  function(.data)
{
  dplyr::rowwise(.data)
}

#' @export
rowwise.tidyseurat <- function(.data)
{
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  .data@meta.data %>%
    
    as_tibble(rownames="cell") %>%
    dplyr::rowwise()
  
}
############# END ADDED tidyseurat #####################################



#' Left join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @export
#'
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidyseurat::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidyseurat::counts %>% left_join(annotation)
#'
left_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
  UseMethod("left_join")
}

#' @export
left_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                ...)
{
  dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
left_join.tidyseurat <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                ...)
{
  
  new_meta = x %>% to_tib() %>% dplyr::left_join( y, by = by, copy = copy, suffix = suffix, ...) 
  
  new_meta %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., "cell") %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        new_meta
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj@meta.data = new_meta %>% data.frame(row.names = "cell")
        new_obj
      } 
    )
  
}

#' Inner join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidyseurat::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidyseurat::counts %>% inner_join(annotation)
#'
#' @export
inner_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
  UseMethod("inner_join")
}

#' @export
inner_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),				 ...)
{
  dplyr::inner_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
inner_join.tidyseurat <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)
{
  new_meta = x %>% to_tib() %>% dplyr::inner_join( y, by = by, copy = copy, suffix = suffix, ...) 
  
  new_meta %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., "cell") %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        new_meta
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj = subset(x,   cells = new_meta %>% pull("cell"))
        new_obj@meta.data = new_meta %>% data.frame(row.names = "cell")
        new_obj
      } 
    )
  
}

#' Right join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidyseurat::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidyseurat::counts %>% right_join(annotation)
#'
#' @export
right_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
  UseMethod("right_join")
}

#' @export
right_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                 ...)
{
  dplyr::right_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
right_join.tidyseurat <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                 ...)
{

  new_meta = x %>% to_tib() %>% dplyr::right_join( y, by = by, copy = copy, suffix = suffix, ...) 
  
  new_meta %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., "cell") %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        new_meta
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj = subset(x,   cells = new_meta %>% pull("cell"))
        new_obj@meta.data = new_meta %>% data.frame(row.names = "cell")
        new_obj
      } 
    )
  
}


#' Full join datasets
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tt object
#'
#' @examples
#'`%>%` = magrittr::`%>%`
#' annotation = tidyseurat::counts %>% distinct(sample) %>% mutate(source = "AU")
#' tidyseurat::counts %>% full_join(annotation)
#'
#' @export
full_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
  UseMethod("full_join")
}

#' @export
full_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                ...)
{
  dplyr::full_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}

#' @export
full_join.tidyseurat <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
                                ...)
{

  new_meta = x %>% to_tib() %>% dplyr::full_join( y, by = by, copy = copy, suffix = suffix, ...) 
  
  new_meta %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., "cell") %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        new_meta
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj@meta.data = new_meta %>% data.frame(row.names = "cell")
        new_obj
      } 
    )
  
}

#' Subset rows using their positions
#'
#' @description
#' `slice()` lets you index rows by their (integer) locations. It allows you
#' to select, remove, and duplicate rows. It is accompanied by a number of
#' helpers for common use cases:
#'
#' * `slice_head()` and `slice_tail()` select the first or last rows.
#' * `slice_sample()` randomly selects rows.
#' * `slice_min()` and `slice_max()` select rows with highest or lowest values
#'   of a variable.
#'
#' If `.data` is a [grouped_df], the operation will be performed on each group,
#' so that (e.g.) `slice_head(df, n = 5)` will select the first five rows in
#' each group.
#'
#' @details
#' Slice does not work with relational databases because they have no
#' intrinsic notion of row order. If you want to perform the equivalent
#' operation, use [filter()] and [row_number()].
#'
#' @family single table verbs
#' @inheritParams arrange
#' @inheritParams filter
#' @param ... For `slice()`: <[`data-masking`][dplyr_data_masking]> Integer row
#'   values.
#'
#'   Provide either positive values to keep, or negative values to drop.
#'   The values provided must be either all positive or all negative.
#'   Indices beyond the number of rows in the input are silently ignored.
#'
#'   For `slice_helpers()`, these arguments are passed on to methods.
#'
#' @param n,prop Provide either `n`, the number of rows, or `prop`, the
#'   proportion of rows to select. If neither are supplied, `n = 1` will be
#'   used.
#'
#'   If `n` is greater than the number of rows in the group (or `prop > 1`),
#'   the result will be silently truncated to the group size. If the
#'   `prop`ortion of a group size is not an integer, it is rounded down.
#' @return
#' An object of the same type as `.data`. The output has the following
#' properties:
#'
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * Data frame attributes are preserved.
#' @section Methods:
#' These function are **generic**s, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' Methods available in currently loaded packages:
#'
#' * `slice()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice")}.
#' * `slice_head()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_head")}.
#' * `slice_tail()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_tail")}.
#' * `slice_min()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_min")}.
#' * `slice_max()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_max")}.
#' * `slice_sample()`: \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("slice_sample")}.
#' @export
#' @examples
#' mtcars %>% slice(1L)
#' # Similar to tail(mtcars, 1):
#' mtcars %>% slice(n())
#' mtcars %>% slice(5:n())
#' # Rows can be dropped with negative indices:
#' slice(mtcars, -(1:4))
#'
#' # First and last rows based on existing order
#' mtcars %>% slice_head(n = 5)
#' mtcars %>% slice_tail(n = 5)
#'
#' # Rows with minimum and maximum values of a variable
#' mtcars %>% slice_min(mpg, n = 5)
#' mtcars %>% slice_max(mpg, n = 5)
#'
#' # slice_min() and slice_max() may return more rows than requested
#' # in the presence of ties. Use with_ties = FALSE to suppress
#' mtcars %>% slice_min(cyl, n = 1)
#' mtcars %>% slice_min(cyl, n = 1, with_ties = FALSE)
#'
#' # slice_sample() allows you to random select with or without replacement
#' mtcars %>% slice_sample(n = 5)
#' mtcars %>% slice_sample(n = 5, replace = TRUE)
#'
#' # you can optionally weight by a variable - this code weights by the
#' # physical weight of the cars, so heavy cars are more likely to get
#' # selected
#' mtcars %>% slice_sample(weight_by = wt, n = 5)
#'
#' # Group wise operation ----------------------------------------
#' df <- tibble(
#'   group = rep(c("a", "b", "c"), c(1, 2, 4)),
#'   x = runif(7)
#' )
#'
#' # All slice helpers operate per group, silently truncating to the group
#' # size, so the following code works without error
#' df %>% group_by(group) %>% slice_head(n = 2)
#'
#' # When specifying the proportion of rows to include non-integer sizes
#' # are rounded down, so group a gets 0 rows
#' df %>% group_by(group) %>% slice_head(prop = 0.5)
#'
#' # Filter equivalents --------------------------------------------
#' # slice() expressions can often be written to use `filter()` and
#' # `row_number()`, which can also be translated to SQL. For many databases,
#' # you'll need to supply an explicit variable to use to compute the row number.
#' filter(mtcars, row_number() == 1L)
#' filter(mtcars, row_number() == n())
#' filter(mtcars, between(row_number(), 5, n()))
slice <- function(.data, ..., .preserve = FALSE) {
  UseMethod("slice")
}
#' @export
slice.default <-  function (.data, ..., .preserve = FALSE)
{
  dplyr::filter(.data, ..., .preserve = .preserve)
}

#' @export
slice.tidyseurat <- function (.data, ..., .preserve = FALSE)
{
  new_meta = dplyr::slice(.data@meta.data, ..., .preserve = .preserve)
  new_obj = subset(.data,   cells = rownames(new_meta ))
  new_obj@meta.data = new_meta
  
  new_obj
  
}

#' Subset columns using their names and types
#'
#' @description
#'
#' Select (and optionally rename) variables in a data frame, using a concise
#' mini-language that makes it easy to refer to variables based on their name
#' (e.g. `a:f` selects all columns from `a` on the left to `f` on the
#' right). You can also use predicate functions like [is.numeric] to select
#' variables based on their properties.
#'
#'
#' ## Overview of selection features
#'
#' ```{r, child = "man/rmd/overview.Rmd"}
#' ```
#'
#' @inheritParams arrange
#' @param ... <[`tidy-select`][dplyr_tidy_select]> One or more unquoted
#'   expressions separated by commas. Variable names can be used as if they
#'   were positions in the data frame, so expressions like `x:y` can
#'   be used to select a range of variables.
#' @return
#' An object of the same type as `.data`. The output has the following
#' properties:
#'
#' * Rows are not affected.
#' * Output columns are a subset of input columns, potentially with a different
#'   order. Columns will be renamed if `new_name = old_name` form is used.
#' * Data frame attributes are preserved.
#' * Groups are maintained; you can't select off grouping variables.
#'
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("select")}.
#'
#' @section Examples:
#' ```{r, child = "man/rmd/setup.Rmd"}
#' ```
#'
#' Here we show the usage for the basic selection operators. See the
#' specific help pages to learn about helpers like [starts_with()].
#'
#' The selection language can be used in functions like
#' `dplyr::select()` or `tidyr::pivot_longer()`. Let's first attach
#' the tidyverse:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' library(tidyverse)
#'
#' # For better printing
#' iris <- as_tibble(iris)
#' ```
#'
#' Select variables by name:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(height)
#'
#' iris %>% pivot_longer(Sepal.Length)
#' ```
#'
#' Select multiple variables by separating them with commas. Note how
#' the order of columns is determined by the order of inputs:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(homeworld, height, mass)
#' ```
#'
#' Functions like `tidyr::pivot_longer()` don't take variables with
#' dots. In this case use `c()` to select multiple variables:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% pivot_longer(c(Sepal.Length, Petal.Length))
#' ```
#'
#' ## Operators:
#'
#' The `:` operator selects a range of consecutive variables:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(name:mass)
#' ```
#'
#' The `!` operator negates a selection:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(!(name:mass))
#'
#' iris %>% select(!c(Sepal.Length, Petal.Length))
#'
#' iris %>% select(!ends_with("Width"))
#' ```
#'
#' `&` and `|` take the intersection or the union of two selections:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% select(starts_with("Petal") & ends_with("Width"))
#'
#' iris %>% select(starts_with("Petal") | ends_with("Width"))
#' ```
#'
#' To take the difference between two selections, combine the `&` and
#' `!` operators:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% select(starts_with("Petal") & !ends_with("Width"))
#' ```
#'
#' @family single table verbs
#' @export
select <- function(.data, ...) {
  UseMethod("select")
}

#' @export
select.default <-  function (.data, ...)
{
  dplyr::select(.data, ...)
}

#' @export
select.tidyseurat <- function (.data, ...)
{
  
  .data %>%
    to_tib() %>%
    select_helper(...) %>%
    when(
      
      # If key columns are missing
      (c("cell",  "orig.ident", "nCount_RNA", "nFeature_RNA") %in% colnames(.)) %>% any %>% `!` ~ {
        message("tidyseurat says: Key columns are missing. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # If valid seurat meta data
      ~ {
        .data@meta.data = (.) %>% data.frame(row.names = "cell")
        .data
      }
    )
  
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
select_helper = function(.data, ...){
  
  loc <- tidyselect::eval_select(expr(c(...)), .data)

  dplyr::select( .data, loc) 
}
