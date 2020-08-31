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
#' pbmc_small %>% tidy %>% arrange(nFeature_RNA)
arrange <- function(.data, ..., .by_group = FALSE) {
  UseMethod("arrange")
}

#' @param .by_group If `TRUE`, will sort first by grouping variable. Applies to
#'   grouped data frames only.
#' @rdname dplyr-methods
#' @export
#'
 
#' @inheritParams arrange
arrange.default <- function(.data, ..., .by_group = FALSE) {
  
  dplyr::arrange(.data, ..., .by_group = .by_group)
  
}

#' @importFrom tibble as_tibble
#' 
#' @export
#' @inheritParams arrange
arrange.tidyseurat <- function(.data, ..., .by_group = FALSE) {
  

  .data@meta.data = 
    .data %>% 
    as_tibble() %>%
    dplyr::arrange(  ..., .by_group = .by_group  ) %>%
    as_meta_data(.data) 
  
  .data
  
}

 

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
#' @param add.cell.ids from Seurat 3.0 A character vector of length(x = c(x, y)). Appends the corresponding values to the start of each objects' cell names.
#' 
#' @return `bind_rows()` and `bind_cols()` return the same type as
#'   the first input, either a data frame, `tbl_df`, or `grouped_df`.
#' @examples
#' `%>%` = magrittr::`%>%`
#' tt = pbmc_small %>% tidy
#' bind_rows(    tt, tt  )
#' 
#' tt_bind = tt %>% select(nCount_RNA ,nFeature_RNA)
#' tt %>% bind_cols(tt_bind)
#'
#' @name bind
NULL


#' @rdname dplyr-methods
#' 
#' @inheritParams bind
#' 
#' @export
#'
bind_rows <- function(..., .id = NULL,  add.cell.ids = NULL) {
  UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(..., .id = NULL,  add.cell.ids = NULL)
{
  dplyr::bind_rows(..., .id = .id)
}

#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
#' 
bind_rows.tidyseurat <- function(..., .id = NULL,  add.cell.ids = NULL)
{
  

    
  tts = flatten_if(dots_values(...), is_spliced)

  # Check if cell with same name
  merge(  tts[[1]] ,  tts[[2]]  ,  add.cell.ids = add.cell.ids) %>% tidy
  
  
}


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


#' distinct
#' 
#' 
#' @param .data A tbl. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#' @param .keep_all If TRUE, keep all variables in .data. If a combination of ... is not distinct, this keeps the first row of values. (See dplyr)
#'
#' @return A tidyseurat object
#'
#' @examples
#'
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% distinct(groups)
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
  
  .data %>%
    as_tibble() %>%
    dplyr::distinct(..., .keep_all = .keep_all)
  
}
 

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
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% filter(groups == "g1")
#'
#' # Learn more in ?dplyr_tidy_eval
 
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
  new_meta = .data %>% as_tibble() %>% dplyr::filter( ..., .preserve = .preserve) %>% as_meta_data(.data)
  new_obj = subset(.data,   cells = rownames(new_meta ))
  #new_obj@meta.data = new_meta
  
  new_obj
                   
}
 

#' Group by one or more variables
#' 
#' @importFrom dplyr group_by_drop_default
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
#' pbmc_small %>% tidy %>% group_by(groups)
#'

 
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
  
  .data %>%
    as_tibble() %>%
    dplyr::group_by( ..., .add = .add, .drop = .drop) 
  
}
 

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
#' pbmc_small %>% tidy %>% summarise(mean(nCount_RNA))
#'
 
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
  
  .data %>%
    as_tibble() %>%
    dplyr::summarise( ...)
  
}
 

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
#' pbmc_small %>% tidy %>% mutate(nFeature_RNA = 1)
#'
 
#' @export
mutate <- function(.data, ...) {
  UseMethod("mutate")
}

#' @export
mutate.default <-  function(.data, ...)
{
  dplyr::mutate(.data, ...)
}


#' @importFrom dplyr mutate
#' @importFrom rlang enquos
#' 
#' @export
mutate.tidyseurat <- function(.data, ...)
{

  # Check that we are not modifying a key column
  cols = enquos(...) %>% names 
  if(intersect(cols, get_special_columns(.data) %>% c(get_needed_columns())) %>% length %>% gt(0))
    stop(sprintf("tidyseurat says: you are trying to mutate a column that is view only %s (it is not present in the meta.data). If you want to mutate a view-only column, make a copy and mutate that one.", get_special_columns(.data) %>% paste(collapse=", ")))

  .data@meta.data =
    .data %>% 
    as_tibble %>%
    dplyr::mutate( ...)  %>% 
    as_meta_data(.data)

  .data
}

 

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
#' pbmc_small %>% tidy %>% rename(s_score = nFeature_RNA) 
#' 
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
  
  # Check that we are not modifying a key column
  cols = tidyselect::eval_select(expr(c(...)), .data@meta.data) 
  if(intersect(cols %>% names, get_special_columns(.data) %>% c(get_needed_columns())) %>% length %>% gt(0))
    stop(sprintf("tidyseurat says: you are trying to rename a column that is view only %s (it is not present in the meta.data). If you want to mutate a view-only column, make a copy and mutate that one.", get_special_columns(.data) %>% paste(collapse=", ")))
  
  .data@meta.data = dplyr::rename( .data@meta.data,  ...)
  
  .data
 
  
}
 

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
#' 
 
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
  
  .data %>%
    as_tibble() %>%
    dplyr::rowwise()
  
}
 

#' Left join datasets
#' 
#' @importFrom dplyr count
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidyseurat object
#'
#' @export
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' 
#' tt = pbmc_small %>% tidy
#' tt %>% left_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2))
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
  
  x %>% 
    as_tibble() %>%
    dplyr::left_join( y, by = by, copy = copy, suffix = suffix, ...) %>%
    
    when(
      
      # If duplicated cells returns tibble
      dplyr::count(., cell) %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        x@meta.data = (.) %>% as_meta_data(x)
        x
      } 
    )
  
}

#' Inner join datasets
#' 
#' @importFrom dplyr pull
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidyseurat object
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' 
#' tt = pbmc_small %>% tidy
#' tt %>% inner_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1))
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
  x %>%
    as_tibble() %>%
    dplyr::inner_join( y, by = by, copy = copy, suffix = suffix, ...)  %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., cell) %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj = subset(x,   cells =  pull(., "cell"))
        new_obj@meta.data = (.) %>% as_meta_data(new_obj)
        new_obj
      } 
    )
  
}

#' Right join datasets
#' 
#' @importFrom dplyr pull
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidyseurat object
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' 
#' tt = pbmc_small %>% tidy
#' tt %>% right_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1))
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

  x %>% 
    as_tibble() %>%
    dplyr::right_join( y, by = by, copy = copy, suffix = suffix, ...) %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., cell) %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        new_obj = subset(x,   cells = (.) %>% pull("cell"))
        new_obj@meta.data = (.) %>% as_meta_data(new_obj)
        new_obj
      } 
    )
  
}


#' Full join datasets
#' 
#' @importFrom dplyr pull
#'
#' @param x tbls to join. (See dplyr)
#' @param y tbls to join. (See dplyr)
#' @param by A character vector of variables to join by. (See dplyr)
#' @param copy If x and y are not from the same data source, and copy is TRUE, then y will be copied into the same src as x. (See dplyr)
#' @param suffix If there are non-joined duplicate variables in x and y, these suffixes will be added to the output to disambiguate them. Should be a character vector of length 2. (See dplyr)
#' @param ... Data frames to combine (See dplyr)
#'
#' @return A tidyseurat object
#'
#' @examples
#' `%>%` = magrittr::`%>%`
#' 
#' tt = pbmc_small %>% tidy
#' tt %>% full_join(tibble::tibble(groups = "g1", other=1:4)) 
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

 x %>% 
    as_tibble() %>%
    dplyr::full_join( y, by = by, copy = copy, suffix = suffix, ...)  %>%
    
    when(
      
      # If duplicated cells returns tibble
      count(., cell) %>% filter(n>1) %>% nrow %>% gt(0) ~ {
        message("tidyseurat says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # Otherwise return updated tidyseurat
      ~ {
        x@meta.data = (.) %>% as_meta_data(x)
        x
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
#' 
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% slice(1)
slice <- function(.data, ..., .preserve = FALSE) {
  UseMethod("slice")
}
#' @export
slice.default <-  function (.data, ..., .preserve = FALSE)
{
  dplyr::slice(.data, ..., .preserve = .preserve)
}

#' @export
slice.tidyseurat <- function (.data, ..., .preserve = FALSE)
{
  new_meta = dplyr::slice(.data@meta.data, ..., .preserve = .preserve)
  new_obj = subset(.data,   cells = rownames(new_meta ))
  #new_obj@meta.data = new_meta
  
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
#' @examples
#' 
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% select(cell, orig.ident )
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
    as_tibble() %>%
    select_helper(...) %>%
    when(
      
      # If key columns are missing
      (get_needed_columns() %in% colnames(.)) %>% all %>% `!` ~ {
        message("tidyseurat says: Key columns are missing. A data frame is returned for independent data analysis.")
        (.)
      },
      
      # If valid seurat meta data
      ~ {
        .data@meta.data = (.) %>% as_meta_data(.data)
        .data
      }
    )
  
}


#' Sample n rows from a table
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge("superseded")}
#' `sample_n()` and `sample_frac()` have been superseded in favour of
#' [slice_sample()]. While they will not be deprecated in the near future,
#' retirement means that we will only perform critical bug fixes, so we recommend
#' moving to the newer alternative.
#'
#' These functions were superseded because we realised it was more convenient to
#' have two mutually exclusive arguments to one function, rather than two
#' separate functions. This also made it to clean up a few other smaller
#' design issues with `sample_n()`/`sample_frac`:
#'
#' * The connection to `slice()` was not obvious.
#' * The name of the first argument, `tbl`, is inconsistent with other
#'   single table verbs which use `.data`.
#' * The `size` argument uses tidy evaluation, which is surprising and
#'   undocumented.
#' * It was easier to remove the deprecated `.env` argument.
#' * `...` was in a suboptimal position.
#'
#' @keywords internal
#' @param tbl A data.frame.
#' @param size <[`tidy-select`][dplyr_tidy_select]>
#'   For `sample_n()`, the number of rows to select.
#'   For `sample_frac()`, the fraction of rows to select.
#'   If `tbl` is grouped, `size` applies to each group.
#' @param replace Sample with or without replacement?
#' @param weight <[`tidy-select`][dplyr_tidy_select]> Sampling weights.
#'   This must evaluate to a vector of non-negative numbers the same length as
#'   the input. Weights are automatically standardised to sum to 1.
#' @param .env DEPRECATED.
#' @param ... ignored
#' @examples
#' 
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% sample_n(50) 
#' pbmc_small %>% tidy %>% sample_frac(0.1)
#' 
#' @return A tidyseurat object
#'
#' @export
sample_n <- function(tbl, size, replace = FALSE, weight = NULL, .env = NULL, ...) {
  UseMethod("sample_n")
}

#' @export
sample_n.default <- function(tbl, size, replace = FALSE, weight = NULL,
                             .env = parent.frame(), ...) {
  tbl %>% sample_n(size, replace = replace, weight = weight, .env = .env, ...)
}

#' @export
sample_n.tidyseurat <- function(tbl, size, replace = FALSE,
                                weight = NULL, .env = NULL, ...) {
  
  lifecycle::signal_superseded("1.0.0", "sample_n()", "slice_sample()")
  
  new_meta = tbl@meta.data %>% dplyr::sample_n( size, replace = replace, weight = weight, .env = .env, ...) 
  new_obj = subset(tbl,   cells = rownames(new_meta ))
  new_obj@meta.data = new_meta
  
  new_obj
  
}

#' @rdname sample_n
#' @export
sample_frac <- function(tbl, size = 1, replace = FALSE, weight = NULL, .env = NULL, ...) {
  UseMethod("sample_frac")
}

#' @export
sample_frac.default <- function(tbl, size, replace = FALSE, weight = NULL,
                             .env = parent.frame(), ...) {
  tbl %>% dplyr::sample_frac(size, replace = replace, weight = weight, .env = .env, ...)
}

#' @export
sample_frac.tidyseurat <- function(tbl, size = 1, replace = FALSE,
                                   weight = NULL, .env = NULL, ...) {

  lifecycle::signal_superseded("1.0.0", "sample_frac()", "slice_sample()")
  
  new_meta = tbl@meta.data %>% dplyr::sample_frac( size, replace = replace, weight = weight, .env = .env, ...) 
  new_obj = subset(tbl,   cells = rownames(new_meta ))
  new_obj@meta.data = new_meta
  
  new_obj
  
}


#' Count observations by group
#'
#' @description
#' `count()` lets you quickly count the unique values of one or more variables:
#' `df %>% count(a, b)` is roughly equivalent to
#' `df %>% group_by(a, b) %>% summarise(n = n())`.
#' `count()` is paired with `tally()`, a lower-level helper that is equivalent
#' to `df %>% summarise(n = n())`. Supply `wt` to perform weighted counts,
#' switching the summary from `n = n()` to `n = sum(wt)`.
#'
#' `add_count()` are `add_tally()` are equivalents to `count()` and `tally()`
#' but use `mutate()` instead of `summarise()` so that they add a new column
#' with group-wise counts.
#'
#' @param x A data frame, data frame extension (e.g. a tibble), or a
#'   lazy data frame (e.g. from dbplyr or dtplyr).
#' @param ... <[`data-masking`][dplyr_data_masking]> Variables to group by.
#' @param wt <[`data-masking`][dplyr_data_masking]> Frequency weights.
#'   Can be `NULL` or a variable:
#'
#'   * If `NULL` (the default), counts the number of rows in each group.
#'   * If a variable, computes `sum(wt)` for each group.
#' @param sort If `TRUE`, will show the largest groups at the top.
#' @param name The name of the new column in the output.
#'
#'   If omitted, it will default to `n`. If there's already a column called `n`,
#'   it will error, and require you to specify the name.
#' @param .drop For `count()`: if `FALSE` will include counts for empty groups
#'   (i.e. for levels of factors that don't exist in the data). Deprecated in
#'   `add_count()` since it didn't actually affect the output.
#' @return
#' An object of the same type as `.data`. `count()` and `add_count()`
#' group transiently, so the output has the same groups as the input.
#' @export
#' @examples
#' 
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% count(groups)
#' 
#' 
count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL, .drop = group_by_drop_default(x)) {
UseMethod("count")
}

#' @export
count.default <- function(x, ..., wt = NULL, sort = FALSE, name = NULL, .drop = group_by_drop_default(x)) {
  if (!missing(...)) {
    out <- dplyr::group_by(x, ..., .add = TRUE, .drop = .drop)
  }
  else {
    out <- x
  }
  out <- dplyr::tally(out, wt = !!enquo(wt), sort = sort, name = name)
  if (is.data.frame(x)) {
    out <- dplyr::dplyr_reconstruct(out, x)
  }
  out}
#' @export
count.tidyseurat <- function(x, ..., wt = NULL, sort = FALSE, name = NULL, .drop = group_by_drop_default(x)) {
  
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  x %>%
    as_tibble() %>%
    dplyr::count(  ..., wt = !!enquo(wt), sort = sort, name = name, .drop = .drop)
  
}

#' Extract a single column
#'
#' `pull()` is similar to `$`. It's mostly useful because it looks a little
#' nicer in pipes, it also works with remote data frames, and it can optionally
#' name the output.
#'
#' @inheritParams arrange
#' @inheritParams tidyselect::vars_pull
#' @param name An optional parameter that specifies the column to be used
#'   as names for a named vector. Specified in a similar manner as \code{var}.
#' @param ... For use by methods.
#' @return A vector the same size as `.data`.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("pull")}.
#' @export
#' @examples
#' 
#' `%>%` = magrittr::`%>%`
#' pbmc_small %>% tidy %>% pull(groups)
#' 
pull <- function(.data, var = -1, name = NULL, ...) {
  ellipsis::check_dots_used()
  UseMethod("pull")
}
#' @export
pull.default <- function(.data, var = -1, name = NULL, ...) {
  var = enquo(var)
  name = enquo(name)
 .data %>% dplyr::pull( var = !!var, name = !!name, ...)
  }
#' @export
pull.tidyseurat <- function(.data, var = -1, name = NULL, ...) {
  var = enquo(var)
  name = enquo(name)
  
  message("tidyseurat says: A data frame is returned for independent data analysis.")
  
  .data %>%
    as_tibble() %>%
    dplyr::pull( var = !!var, name = !!name, ...)
  
  
}
