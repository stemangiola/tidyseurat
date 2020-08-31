#' @importFrom tibble as_tibble
#' 
#' @keywords internal
#' 
#' @param .data A tidyseurat
to_tib = function(.data){ .data@meta.data %>% as_tibble(rownames = "cell") }

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

# Raise to the power
pow = function(a,b){	a^b }

# Equals
eq = function(a,b){	a==b }

prepend = function (x, values, before = 1) 
{
  n <- length(x)
  stopifnot(before > 0 && before <= n)
  if (before == 1) {
    c(values, x)
  }
  else {
    c(x[1:(before - 1)], values, x[before:n])
  }
}
#' Add class to abject
#' 
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {
  
  if(!name %in% class(var)) class(var) <- prepend(class(var),name)
  
  var
}

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' @keywords internal
drop_class = function(var, name) {
  class(var) <- class(var)[!class(var)%in%name]
  var
}

#' get abundance long
#' 
#' @keywords internal
#' 
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' 
#' @param .data A tidyseurat
#' @param transcripts A character
#' @param all A boolean
#' 
#' 
#' @return A tidyseurat object
#' 
#'
#' @export
get_abundance_sc_wide = function(.data, transcripts = NULL, all = FALSE){
  
  # Solve CRAN warnings
  . = NULL
  
  # Check if output would be too big without forcing
  if(
    length(VariableFeatures(.data)) == 0  &
    is.null(transcripts) &
    all == FALSE
  ) stop("
				 Your object do not contain variable trancript labels,
				 transcript argument is empty and all argument is set to FALSE.
				 Either:
				 1. use detect_variable_features() to select variable feature
				 2. pass an array of transcripts names
				 3. set all = TRUE (this will output a very large object, do you computer have enough RAM?)
				 ")
  
  # Get variable features if existing
  if(
    length(VariableFeatures(.data)) > 0  &
    is.null(transcripts) &
    all == FALSE
  ) variable_genes = VariableFeatures(.data)
  
  # Else
  else variable_genes = NULL
  
  # Just grub last assay
  .data@assays %>%
    tail(1) %>%
    .[[1]] %>%
    when(
      variable_genes %>% is.null %>% `!` ~ (.)@counts[variable_genes,, drop=FALSE],
      transcripts %>% is.null %>% `!` ~ (.)@counts[transcripts,, drop=FALSE],
      ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
    ) %>%
    as.matrix() %>%
    t %>%
    as_tibble(rownames = "cell") 
  
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
#' @param transcripts A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#' 
#' @return A tidyseurat object
#'
#' @export
get_abundance_sc_long = function(.data, transcripts = NULL, all = FALSE, exclude_zeros = FALSE){
  
  # Solve CRAN warnings
  . = NULL
  
  # Check if output would be too big without forcing
  if(
    length(VariableFeatures(.data)) == 0  &
    is.null(transcripts) &
    all == FALSE
  ) stop("
				 Your object do not contain variable trancript labels,
				 transcript argument is empty and all argument is set to FALSE.
				 Either:
				 1. use detect_variable_features() to select variable feature
				 2. pass an array of transcripts names
				 3. set all = TRUE (this will output a very large object, do you computer have enough RAM?)
				 ")
  
   
  # Get variable features if existing
  if(
    length(VariableFeatures(.data)) > 0  &
    is.null(transcripts) &
    all == FALSE
  ) variable_genes = VariableFeatures(.data)
  
  # Else
  else variable_genes = NULL
  
  assay_names = .data@assays %>% names
  
  
  .data@assays %>%
    
    # Take active assay
    map2(assay_names,
         
         ~ .x %>%
           when(
             variable_genes %>% is.null %>% `!` ~ .x@data[variable_genes,, drop=FALSE],
             transcripts %>% is.null %>% `!` ~ .x@data[ toupper(rownames(.x@data)) %in% toupper(transcripts),, drop=FALSE],
             all  ~ .x@data,
             ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
           ) %>%
           
           # Replace 0 with NA
           when(exclude_zeros ~ (.) %>% { x = (.); x[x == 0] <- NA; x }, ~ (.)) %>%
           
           data.frame() %>%
           as_tibble(rownames = "transcript") %>%
           tidyr::pivot_longer(
             cols = -transcript,
             names_to ="cell", 
             values_to = "abundance" %>% paste(.y, sep="_"),
             values_drop_na  = TRUE
           ) 
         #%>%
         #mutate_if(is.character, as.factor) %>%
         
         
    ) %>%
    Reduce(function(...) left_join(..., by=c("transcript", "cell")), .)
  
}

#' @importFrom dplyr select_if
#' 
#' @keywords internal
#' 
#' @param .data A tibble
#' @param seurat_object A tidyseurat
#' 
as_meta_data = function(.data, seurat_object){
  
  # Solve CRAN warnings
  . = NULL
  
  col_to_exclude =  get_special_columns(seurat_object)
  
  .data %>% 
    select_if(!colnames(.) %in% col_to_exclude) %>%
    #select(-one_of(col_to_exclude)) %>%
    data.frame(row.names = "cell")
}

#' @importFrom purrr map_chr
#' 
#' @keywords internal
#' 
#' @param seurat_object A tidyseurat
#' 
get_special_columns = function(seurat_object){
  get_special_datasets(seurat_object) %>%
    map(~ .x %>% colnames  ) %>%
    unlist %>% 
    as.character 
}

get_special_datasets = function(seurat_object){
  seurat_object@reductions %>%
    map(~ .x@cell.embeddings[,1:min(5, ncol(.x@cell.embeddings)), drop=FALSE] )
  
}

get_needed_columns = function(){
  #c("cell",  "orig.ident", "nCount_RNA", "nFeature_RNA")
  c("cell")
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
  
  v = quo_name(quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) %>% 
    strsplit(', ') %>% 
    unlist 
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
select_helper = function(.data, ...){
  
  loc <- tidyselect::eval_select(expr(c(...)), .data)
  
  dplyr::select( .data, loc) 
}