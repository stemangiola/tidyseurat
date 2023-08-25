#' @name unnest
#' @rdname unnest
#' @inherit tidyr::unnest
#' @aliases unnest_seurat
#' @return `tidyseurat`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> 
#'     nest(data=-groups) |> 
#'     unnest(data)
#'
#' @importFrom rlang quo_name
#' @importFrom purrr imap
#' @importFrom tidyr unnest
#' @export
unnest.tidyseurat_nested <- function(data, cols, ...,
    keep_empty=FALSE, ptype=NULL, names_sep=NULL, 
    names_repair="check_unique", .drop, .id, .sep, .preserve) {

    cols <- enquo(cols)

    unnest_seurat(data, !!cols, ...,
        keep_empty=keep_empty, ptype=ptype,
        names_sep=names_sep, names_repair=names_repair)

}

#' @rdname unnest
#' @importFrom tidyr unnest
#' @importFrom purrr when
#' @importFrom rlang quo_name
#' @importFrom purrr imap
#' @export
unnest_seurat  <-  function(data, cols, ...,
    keep_empty=FALSE, ptype=NULL,
    names_sep=NULL, names_repair="check_unique",
    .drop, .id, .sep, .preserve) {
    # Need this otherwise crashes map
    .data_ <- data
  
    cols <- enquo(cols)
  
    .data_ %>% 
        when(
      
            # If my only column to unnest is tidyseurat
            pull(., !!cols) %>% .[[1]] %>% is("Seurat") %>% any ~  
            {
                # Do my trick to unnest
                list_seurat <- mutate(.,
                    !!cols := imap(
                        !!cols, ~ .x %>%
                            bind_cols_(
                                .data_ %>%
                                    select(-!!cols) %>%
                                    slice(rep(.y, nrow(as_tibble(.x))))
                            )
                        )) %>%
                    pull(!!cols)
                list_seurat[[1]] %>%
                    # Bind only if length list > 1
                    when(
                        length(list_seurat)>1 ~ bind_rows(.,
                            list_seurat[2:length(list_seurat)]),
                        ~ (.)
                    )
            },
      
            # Else do normal stuff
            ~ (.) %>% 
                drop_class("tidyseurat_nested") %>%
                tidyr::unnest(!!cols, ..., keep_empty=keep_empty,
                    ptype=ptype, names_sep=names_sep,
                    names_repair=names_repair) %>%
                add_class("tidyseurat_nested"))
}


#' @name nest
#' @rdname nest
#' @inherit tidyr::nest
#' @return `tidyseurat_nested`
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> 
#'     nest(data=-groups) |> 
#'     unnest(data)
#' 
#' @importFrom tidyr nest
#' @importFrom magrittr equals
#' @importFrom rlang enquos
#' @importFrom Seurat SplitObject
#' @importFrom rlang :=
#' @export
nest.Seurat <- function (.data, ..., .names_sep=NULL)
{
    cols <- enquos(...)
    col_name_data  <- names(cols)

    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>%
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
  
    my_data__ <- .data
  
    # This is for getting the column names
    dummy_nested <- 
        my_data__[1,] |>
        to_tib() %>%
        tidyr::nest(...)
  
    split_by_column <- 
        dummy_nested |> 
        select(-col_name_data) |>
        colnames()
  
    # If nesting on one group use the fast split
    if (split_by_column |> length() |> identical(1L))
  
        my_data__ |> 
            SplitObject(split.by=split_by_column) |>
            map(~ .x |> select(-split_by_column)) |> 
            enframe(name=split_by_column, value=col_name_data) |>
            # Coerce to tidyseurat_nested for unnesting
            add_class("tidyseurat_nested")
  
    # If arbitrary nest is needed use the slow one
    else
        my_data__ %>%
        # This is needed otherwise nest goes into loop and fails
            to_tib %>%
            tidyr::nest(...) %>%
      
            mutate(
                !!as.symbol(col_name_data) := map(
                    !!as.symbol(col_name_data),
                    ~ my_data__ %>% 
                    # Subset cells
                    filter(!!c_(.data)$symbol %in% 
                        pull(.x, !!c_(.data)$symbol)) %>%
                    # Subset columns
                    select(colnames(.x))
                )) |>
      
            # Coerce to tidyseurat_nested for unnesting
            add_class("tidyseurat_nested")
}

#' @name extract
#' @rdname extract
#' @inherit tidyr::extract
#' @return `tidyseurat`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |>
#'   extract(groups, 
#'     into="g", 
#'     regex="g([0-9])", 
#'     convert=TRUE)
#' 
#' @importFrom tidyr extract
#' @export
extract.Seurat <- function  (data, col, into,
    regex="([[:alnum:]]+)", remove=TRUE, convert=FALSE, ...) {
	
    col <- enquo(col)
	
	# Deprecation of special column names
	if (is_sample_feature_deprecated_used(
	    data, 
	    c(quo_name(col), into)
	)) {
	    data= ping_old_special_column_into_metadata(data)
	}
	
	data@meta.data <- 
	    data %>%
	    as_tibble() %>%
	    tidyr::extract(col=!!col, into=into, regex=regex,
            remove=remove, convert=convert, ...) %>%
	    as_meta_data(data)

	data
}

#' @name pivot_longer
#' @rdname pivot_longer
#' @inherit tidyr::pivot_longer
#' @return `tidyseurat`
#' 
#' @export
#' @examples
#' data(pbmc_small)
#' pbmc_small |> pivot_longer(
#'   cols=c(orig.ident, groups),
#'   names_to="name", values_to="value")
#' 
#' @importFrom ellipsis check_dots_used
#' @importFrom tidyr pivot_longer
#' @export
pivot_longer.Seurat <- function(data,
    cols, names_to="name", names_prefix=NULL,
    names_sep=NULL, names_pattern=NULL, names_ptypes=NULL,
    names_transform=NULL, names_repair="check_unique",
    values_to="value", values_drop_na=FALSE,
    values_ptypes=NULL, values_transform=NULL, ...) {
    cols <- enquo(cols) 
  
    message(data_frame_returned_message)
  
    # Deprecation of special column names
    if (is_sample_feature_deprecated_used(
        data, 
        c(quo_names(cols))
    )) {
        data= ping_old_special_column_into_metadata(data)
    }
  
    data %>%
        as_tibble() %>%
        tidyr::pivot_longer(!!cols, names_to=names_to,
            names_prefix=names_prefix, names_sep=names_sep,
            names_pattern=names_pattern, names_ptypes=names_ptypes,
            names_transform=names_transform, names_repair=names_repair,
            values_to=values_to, values_drop_na=values_drop_na,
            values_ptypes=values_ptypes, values_transform=values_transform,
            ...)
}

#' @name unite
#' @rdname unite
#' @inherit tidyr::unite
#' @return `tidyseurat`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> unite(
#'   col="new_col", 
#'   c("orig.ident", "groups"))
#'     
#' @importFrom rlang enquo enquos quo_name
#' @importFrom tidyr unite
#' @export
unite.Seurat <- function(data, col,
    ..., sep="_", remove=TRUE, na.rm=FALSE) {
  
    # Check that we are not modifying a key column
    cols <- enquo(col) 
  
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }

    .view_only_cols <- c(
        get_special_columns(data),
        get_needed_columns(data))
    
    .test <- intersect(
        quo_names(cols), 
        .view_only_cols)

    if (remove && length(.test)) {
        stop("tidyseurat says:",
            " you are trying to rename a column",
            " that is view only ", 
            paste(.view_only_cols, collapse=", "),
            " (it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one.")
    }
  
    data@meta.data <- data %>%
        as_tibble() %>%
        tidyr::unite(!!cols, ..., sep=sep,
            remove=remove, na.rm=na.rm) %>%
        as_meta_data(data)
  
    data
}

#' @name separate
#' @rdname separate
#' @inherit tidyr::separate
#' @return `tidyseurat`
#' 
#' @examples
#' data(pbmc_small)
#' un <- pbmc_small |> unite("new_col", c(orig.ident, groups))
#' un |> separate(new_col, c("orig.ident", "groups"))
#' 
#' @importFrom tidyr separate
#' @export
separate.Seurat <- function(data, col, into,
    sep="[^[:alnum:]]+", remove=TRUE, convert=FALSE,
    extra="warn", fill="warn", ...) {
  
    # Check that we are not modifying a key column
    cols <- enquo(col)
  
    # Deprecation of special column names
    if(is_sample_feature_deprecated_used(
        data, 
        c(quo_names(cols))
    )) {
        data= ping_old_special_column_into_metadata(data)
    }

    .view_only_cols <- c(
        get_special_columns(data),
        get_needed_columns(data))
    
    .test <- intersect(
        quo_names(cols), 
        .view_only_cols)

    if (remove && length(.test)) {
        stop("tidyseurat says:",
            " you are trying to rename a column",
            " that is view only ",
            paste(.view_only_cols, collapse=", "),
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one.")
    }
   
    data@meta.data =
        data %>%
        as_tibble() %>% 
        tidyr::separate(!!cols, into=into, sep=sep, remove=remove, 
            convert=convert, extra=extra, fill=fill, ...) %>%
        as_meta_data(data)
    data
}