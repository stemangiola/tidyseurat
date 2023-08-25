#' Cell types of 80 PBMC single cells
#' 
#' A dataset containing the barcodes and cell types of 80 PBMC single cells.
#'
#' @format A tibble containing 80 rows and 2 columns.
#'   Cells are a subsample of the Peripheral Blood Mononuclear Cells (PBMC) 
#'   dataset of 2,700 single cell. Cell types were identified with SingleR.
#' \describe{
#'   \item{cell}{cell identifier, barcode}
#'   \item{first.labels}{cell type}
#' }
#' @source \url{https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html}
#' @usage data(cell_type_df)
#' @return `tibble`
"cell_type_df"

#' Intercellular ligand-receptor interactions for 
#' 38 ligands from a single cell RNA-seq cluster.
#'
#' A dataset containing ligand-receptor interactions within a sample.
#' There are 38 ligands from a single cell cluster versus 35 receptors 
#' in 6 other clusters.
#'
#' @format A `tibble` containing 100 rows and 9 columns.
#'   Cells are a subsample of the PBMC dataset of 2,700 single cells. 
#'   Cell interactions were identified with `SingleCellSignalR`.
#' \describe{
#'   \item{sample}{sample identifier}
#'   \item{ligand}{cluster and ligand identifier}
#'   \item{receptor}{cluster and receptor identifier}
#'   \item{ligand.name}{ligand name}
#'   \item{receptor.name}{receptor name}
#'   \item{origin}{cluster containing ligand}
#'   \item{destination}{cluster containing receptor}
#'   \item{interaction.type}{type of interation, paracrine or autocrine}
#'   \item{LRscore}{interaction score}
#' }
#' @source \url{https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html}
#' @usage data(pbmc_small_nested_interactions)
#' @return `tibble`
"pbmc_small_nested_interactions"