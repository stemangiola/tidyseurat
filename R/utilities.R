#' @importFrom tibble as_tibble
#'
#' @keywords internal
#'
#' @param .data A tidyseurat
#' 
#' @noRd
to_tib <- function(.data) {
    .data[[]] %>%
        as_tibble(rownames=c_(.data)$name)
}

# Greater than
gt <- function(a, b) {
    a > b
}

# Smaller than
st <- function(a, b) {
    a < b
}

# Negation
not <- function(is) {
    !is
}

# Raise to the power
pow <- function(a, b) {
    a^b
}

# Equals
eq <- function(a, b) {
    a == b
}

prepend <- function(x, values, before=1) {
    n <- length(x)
    stopifnot(before > 0 && before <= n)
    if (before == 1) {
        c(values, x)
    } else {
        c(x[seq_len(before-1)], values, x[seq(before, n)])
    }
}

#' Add class to abject
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class <- function(var, name) {
    if (!name %in% class(var)) 
        class(var) <- prepend(class(var), name)
    return(var)
}

#' Remove class to abject
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' @keywords internal
drop_class <- function(var, name) {
    class(var) <- class(var)[!class(var) %in% name]
    return(var)
}

#' get abundance wide
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat DietSeurat
#' @importFrom SeuratObject DefaultAssay<-
#' @importFrom stats setNames
#'
#' @param .data A tidyseurat
#' @param features A character
#' @param all A boolean
#' @param assay assay name to extract feature abundance
#' @param slot slot in the assay, e.g. `data` and `scale.data`
#' @param prefix prefix for the feature names
#'
#' @return A Seurat object
#' @examples
#' data(pbmc_small)
#' pbmc_small %>%
#'   get_abundance_sc_wide(features=c("HLA-DRA", "LYZ"))
#'
#' @export
get_abundance_sc_wide <- function(.data, features=NULL, all=FALSE,
    assay=.data@active.assay, slot="data", prefix="") {

    # Solve CRAN warnings
    . <- NULL
    assays <- NULL
    counts <- NULL
  
    if (is.null(assay)) {
  	    assay <- .data@active.assay
    }

    # Check if output would be too big without forcing
    if(
        length(VariableFeatures(.data)) == 0  &
        is.null(features) &
        all == FALSE
    ) {
        stop("Your object do not contain variable trancript labels,\n",
			 " feature argument is empty and all argument is set to FALSE.\n",
			 " Either:\n",
			 " 1. use detect_variable_features() to select variable feature\n",
			 " 2. pass an array of features names\n",
			 " 3. set all=TRUE (this will output a very large object;",
             " does your computer have enough RAM?)\n")
    }

    # Get variable features if existing
    if(
        length(VariableFeatures(.data)) > 0  &
        is.null(features) &
        all == FALSE
    ) variable_genes <- VariableFeatures(.data)
    # Else
    else variable_genes <- NULL

    # Eliminate unneeded assays.
    # This because if a gene is not in an assay I am not interested about
    # this could cause an unneeded error
    DefaultAssay(.data) <- assay
    .data = .data |> DietSeurat(assays = assay)

    # Just grub last assay
    .data |> 
      GetAssayData(assay = assay, layer=slot) %>% 
        when(
            variable_genes %>% is.null %>% `!` ~ 
                (.)[ toupper(rownames(.)) %in% toupper(variable_genes),,drop=FALSE],
            features %>% is.null %>% `!` ~ 
                (.)[ toupper(rownames(.)) %in% toupper(features),,drop=FALSE],
            ~ stop("tidyseurat says: It is not convenient to",
                " extract all genes, you should have either variable",
                " features or feature list to extract.")
        ) |> 
        as.matrix() |> 
        t() |> 
        as_tibble(rownames=c_(.data)$name) %>%

        # Add prefix
        setNames(c(c_(.data)$name, sprintf("%s%s", prefix, colnames(.)[-1])))
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom Seurat VariableFeatures
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom purrr map2
#'
#' @param .data A tidyseurat
#' @param features A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#' @param assay assay name to extract feature abundance
#' @param slot slot in the assay, e.g. `data` and `scale.data`
#'
#' @return A Seurat object
#' @examples
#' data(pbmc_small)
#' pbmc_small %>%
#'   get_abundance_sc_long(features=c("HLA-DRA", "LYZ"))
#'
#' @export
get_abundance_sc_long <- function(.data, features=NULL, all=FALSE,
    exclude_zeros=FALSE, assay=Assays(.data), slot="data"){

    # Solve CRAN warnings
    . <- NULL
  
    if (is.null(assay)) {
  	    assay <- Assays(.data)
    }
  
    # Check if output would be too big without forcing
    if (
        length(VariableFeatures(.data)) == 0  &
        is.null(features) &
        all == FALSE
    ) {
        stop("Your object do not contain variable trancript labels,\n",
        " feature argument is empty and all argument is set to FALSE.\n",
        " Either:\n",
        " 1. use detect_variable_features() to select variable feature\n",
        " 2. pass an array of features names\n",
        " 3. set all=TRUE (this will output a very large object;",
        " does your computer have enough RAM?)\n")
    }


    # Get variable features if existing
    if(
        length(VariableFeatures(.data)) > 0  &
        is.null(features) &
        all == FALSE
    ) variable_genes <- VariableFeatures(.data)
    # Else
    else variable_genes <- NULL

    .data@assays %>%
  	    .[assay] %>%
        # Take active assay
        map2(assay,
            ~ .x %>%
                GetAssayData(layer = slot) %>%
                when(
                    variable_genes %>% is.null %>% `!` ~
                        (.)[variable_genes,, drop=FALSE],
                    features %>% is.null %>% `!` ~ 
                        (.)[ toupper(rownames((.))) %in% 
                            toupper(features), , drop=FALSE],
                    all ~ (.),
                    ~ stop("tidyseurat says: It is not convenient to",
                        " extract all genes, you should have either variable",
                        " features or feature list to extract.")
                ) %>%

                # Replace 0 with NA
                when(exclude_zeros ~ 
                        (.) %>%
                        { x=(.); x[x == 0] <- NA; x }, ~ (.)) %>%
                        data.frame(check.names=FALSE) %>%
                        as_tibble(rownames=".feature") %>%
                        tidyr::pivot_longer(
                            cols= - .feature,
                            names_to=c_(.data)$name,
                            values_to=".abundance" %>% paste(.y, sep="_"),
                            values_drop_na=TRUE
                        ) #%>%
                #mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...)
            full_join(..., by=c(".feature", c_(.data)$name)), .)
}

#' @importFrom dplyr select_if
#' @importFrom tibble column_to_rownames
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param seurat_object A tidyseurat
#'
#' @noRd
as_meta_data <- function(.data, seurat_object){

    # Solve CRAN warnings
    . <- NULL

    col_to_exclude <- get_special_columns(seurat_object)

    .data %>%
        select_if(!colnames(.) %in% col_to_exclude) %>%
        #select(-one_of(col_to_exclude)) %>%
        column_to_rownames(c_(seurat_object)$name)
}

#' @importFrom purrr map_chr
#'
#' @keywords internal
#'
#' @param seurat_object A tidyseurat
#'
#' @noRd
get_special_columns <- function(seurat_object){
    get_special_datasets(seurat_object) %>%
        map(~ .x %>% colnames  ) %>%
        unlist %>%
        as.character
}

get_special_datasets <- function(seurat_object, n_dimensions_to_return=Inf){
    seurat_object@reductions %>%
        map(~ .x@cell.embeddings[,
            1:min(n_dimensions_to_return, ncol(.x@cell.embeddings)),
            drop=FALSE])
}

get_needed_columns <- function(.data){
    c(c_(.data)$name)
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {
    v <- quo_name(quo_squash(v))
    gsub('^c\\(|`|\\)$', '', v) %>%
        strsplit(', ') %>%
        unlist
}


#' returns variables from an expression
#' @param expression an expression
#' @importFrom rlang enexpr
#' @return list of symbols
return_arguments_of <- function(expression){
    variables <- enexpr(expression) |> as.list()
    if(length(variables) > 1) {
        variables <- variables[-1] # removes first element which is function
    }
    variables
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
select_helper <- function(.data, ...){
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    dplyr::select( .data, loc)
}

data_frame_returned_message <- paste(
    "tidyseurat says:",
    "A data frame is returned for independent data analysis.")

duplicated_cell_names <- paste(
    "tidyseurat says:",
    "This operation lead to duplicated cell names.",
    "A data frame is returned for independent data analysis.")

#' @importFrom methods .hasSlot
clean_seurat_object <- function(.data){

    . <- NULL

    if (.hasSlot(.data, "images"))
        .data@images <-
            map(.data@images,
                ~ .x %>% when((.)@coordinates %>% nrow() %>% gt(0) ~ (.))) %>%

                # Drop NULL
                Filter(Negate(is.null), .)

        .data@assays <- .data@assays %>%
            map(~ {
                my_assay=.x
                if (.hasSlot(., "SCTModel.list"))
                    my_assay@SCTModel.list  =
                    map(my_assay@SCTModel.list,
                        ~ .x %>%
                            when((.)@cell.attributes %>%
                                nrow() %>% gt(0) ~ (.))) %>%

                            # Drop NULL
                            Filter(Negate(is.null), .)
                my_assay
            })

    .data
}


# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and
# is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
is_sample_feature_deprecated_used <- function(.data, 
    user_columns, use_old_special_names=FALSE) {
    
    cell <- any(str_detect(user_columns, regex("\\bcell\\b")))
    .cell <- any(str_detect(user_columns, regex("\\W*(\\.cell)\\W*")))
    
    old_standard_is_used <- 
        !"cell" %in% colnames(.data@meta.data) &&
        ("cell" %in% user_columns || (cell && !.cell))
    
    if (old_standard_is_used) {
        warning("tidyseurat says:",
            " from version 1.3.1, the special columns including",
            " cell id (colnames(se)) has changed to \".cell\".",
            " This dataset is returned with the old-style vocabulary (cell),",
            " however, we suggest to update your workflow",
            " to reflect the new vocabulary (.cell).")
        use_old_special_names <- TRUE
    }
    use_old_special_names
}

get_special_column_name_symbol <- function(name){
    list(name=name, symbol=as.symbol(name))
}

# Key column names
ping_old_special_column_into_metadata <- function(.data){
    .data@misc$cell__ <- get_special_column_name_symbol("cell")
    .data
}

get_special_column_name_cell <- function(name){
    list(name=name, symbol=as.symbol(name))
}

cell__ <- get_special_column_name_symbol(".cell")

c_ <- function(x){
    # Check if old deprecated columns are used
    if("cell__" %in% names(x@misc)) cell__ <- x@misc$cell__
    return(cell__)
}

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr vars
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr <- function(var, attribute, name) {
    attr(var, name) <- attribute
    var
}

#' Get specific annotation columns
#'
#' @keywords internal
#' @noRd
#' 
#' @importFrom rlang enquo
#' @importFrom purrr map
#' @importFrom dplyr distinct_at
#' @importFrom magrittr equals
#' @importFrom dplyr vars
#' 
#' @param .data A tibble
#' @param .col A vector of column names
#' 
#' @return A character
get_specific_annotation_columns <- function(.data, .col) {
    
    # Comply with CRAN NOTES
    . <- NULL
    
    # Make col names
    .col <- enquo(.col)
    
    # x-annotation df
    n_x <- .data |> distinct_at(vars(!!.col)) |> nrow()
    
    # element wise columns
    .data |>
        select(-!!.col) |>
        colnames() |>
        map(~ {
            n_.x <- .data |> distinct_at(vars(!!.col, .x)) |> nrow()
            if (n_.x == n_x) .x else NULL
        }) %>%
        # Drop NULL
        { (.)[lengths((.)) != 0] } |>
        unlist()
}

subset_tidyseurat <- function(.data, .column) {
    # Make col names
    .column <- enquo(.column)

    # Check if column present
    if (.data |> select(!!.column) |> colnames() %in% colnames(.data) %>% all %>% `!`)
        stop("tidyseurat says: some of the .column specified",
            " do not exist in the input data frame.")


    .data %>%
    # Selecting the right columns
        select(!!.column, get_specific_annotation_columns(.data, !!.column)) %>%
        distinct()
}

#' @importFrom Seurat GetAssayData
#' @importFrom methods is
GetAssayData_robust = function(seurat_assay, layer = NULL){
  
  if(
    seurat_assay |> is("Assay5") & 
    seurat_assay |> ncol() == 1
  ){
    m = seurat_assay@layers[[layer]] |> as.matrix()
    rownames(m) = rownames(seurat_assay)
    colnames(m) = colnames(seurat_assay)
    m
  }
    
  else 
    GetAssayData(seurat_assay, layer=layer)
}
