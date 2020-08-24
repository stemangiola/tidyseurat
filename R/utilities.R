#' @importFrom tibble as_tibble
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

#' Add class to abject
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
#' @importFrom magrittr "%$%"
#'
#' @export
get_abundance_sc_wide = function(.data, transcripts = NULL, all = F){
  
  
  
  # Check if output would be too big without forcing
  if(
    length(VariableFeatures(.data)) == 0  &
    is.null(transcripts) &
    all == F
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
    all == F
  ) variable_genes = VariableFeatures(.data)
  
  # Else
  else variable_genes = NULL
  
  # Just grub last assay
  my_assay = .data@assays %>% names  %>% tail(1)
  
  
  my_assay %>%
    when(
      variable_genes %>% is.null %>% `!` ~ .x@counts[variable_genes,],
      transcripts %>% is.null %>% `!` ~ .x@counts[transcripts,],
      ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
    ) %>%
    as.matrix() %>%
    t %>%
    as_tibble(rownames = quo_name(.transcript)) 
  
}

#' get abundance long
#' @importFrom magrittr "%$%"
#' @importFrom Seurat VariableFeatures
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom purrr map2
#' 
#'
#' @export
get_abundance_sc_long = function(.data, transcripts = NULL, all = F, exclude_zeros = F){
  
  
  
  # Check if output would be too big without forcing
  if(
    length(VariableFeatures(.data)) == 0  &
    is.null(transcripts) &
    all == F
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
    all == F
  ) variable_genes = VariableFeatures(.data)
  
  # Else
  else variable_genes = NULL
  
  assay_names = .data@assays %>% names
  
  
  .data@assays %>%
    
    # Take active assay
    map2(assay_names,
         
         ~ .x %>%
           when(
             variable_genes %>% is.null %>% `!` ~ .x@data[variable_genes,, drop=F],
             transcripts %>% is.null %>% `!` ~ .x@data[ toupper(rownames(.x@data)) %in% toupper(transcripts),, drop=F],
             all  ~ .x@data,
             ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
           ) %>%
           
           # Replace 0 with NA
           when(exclude_zeros ~ (.) %>% { x = (.); x[x == 0] <- NA; x }, ~ (.)) %>%
           
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
