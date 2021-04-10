# Article workflow
library(tidyverse)
library(Seurat)
library(SingleR)
library(plotly)
library(tidyHeatmap)
library(tidyseurat)
options(future.globals.maxSize = 50068 * 1024^2)

PBMC <- readRDS("dev/PBMC_integrated.rds")

# Polishing

PBMC_clean <-
  PBMC_

  # Clean groups
  mutate(Phase = Phase %>% str_remove("^phase_")) %>%

  # Extract sample
  extract(sample, "sample", "./data/seurat/outs/([a-zA-Z0-9]+)")

# PBMC_clean = PBMC_clean %>% nest(data = -sample) %>% mutate(data = map(data, ~ .x %>% sample_n(200))) %>% unnest(data)

# Scaling
# PBMC_clean_scaled <-
#   PBMC_clean %>%
#    SCTransform(verbose = FALSE) %>%
#   FindVariableFeatures(verbose = FALSE)


# Dimensionality reduction
PBMC_clean_scaled_UMAP <-
  PBMC_clean %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, n.components = 3L)

# Clustering
PBMC_clean_scaled_UMAP_cluster <-
  PBMC_clean_scaled_UMAP %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(method = "igraph", verbose = FALSE)

# Cell_type classification Manual
markers <-
  PBMC_clean_scaled_UMAP_cluster %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC)

# Cell_type classification Automatic

# Get cell type reference data
hpca <- HumanPrimaryCellAtlasData()

# Infer cell identities
cell_type_df <-

  # extracting counts from Seurat object
  GetAssayData(PBMC_clean_scaled_UMAP_cluster, slot = 'counts', assay = "SCT") %>%
  log1p() %>%

  # SingleR
  SingleR(
    ref = hpca,
    labels = hpca$label.main,
    method = "cluster",
    clusters = PBMC_clean_scaled_UMAP_cluster %>% pull(seurat_clusters)
  ) %>%

  # Formatting results
  as.data.frame() %>%
  as_tibble(rownames = "seurat_clusters") %>%
  select(seurat_clusters, first.labels)


# Infer cell identities - cell wise
cell_type_df_single <-

  # extracting counts from Seurat object
  GetAssayData(PBMC_clean_scaled_UMAP_cluster, slot = 'counts', assay = "SCT") %>%
  log1p() %>%

  # SingleR
  SingleR(
    ref = hpca,
    labels = hpca$label.main,
    method = "single"
  ) %>%

  # Formatting results
  as.data.frame() %>%
  as_tibble(rownames = "cell") %>%
  select(cell, first.labels_single = first.labels)

# Join UMAP and cell type info
PBMC_clean_scaled_UMAP_cluster_cell_type <-
  PBMC_clean_scaled_UMAP_cluster %>%
  left_join(
    cell_type_df,
    by = "seurat_clusters"
  ) %>%
  left_join(
    cell_type_df_single,
    by = "cell"
  )


# Markers
PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) %>%
  saveRDS("dev/PBMC_marker_df.rds")

# Nesting
PBMC_clean_scaled_UMAP_cluster_cell_type  %>%
  sample_n(1000) %>%

  # Label lymphoid and myeloid
  tidyseurat::filter(first.labels != "Platelets") %>%
  tidyseurat::mutate(cell_class =
                       if_else(
                         `first.labels` %in% c("Macrophage", "Monocyte"),
                         "myeloid",
                         "lymphoid"
                       )
  ) %>%

  # Nesting
  nest(data = -cell_class) %>%

  # Identification of variable gene features
  mutate(variable_genes = map_chr(
    data, ~ .x %>%
      FindVariableFeatures() %>%
      RunPCA(verbose = FALSE) %>%
      FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
      pull(gene) %>%
      head() %>%
      paste(collapse=", ")
  ))

# # Reorder columns
# PBMC_clean_scaled_UMAP_cluster_cell_type %>%
#   count(seurat_clusters, first.labels_cluster = first.labels)

saveRDS(PBMC_clean_scaled_UMAP_cluster_cell_type, "dev/PBMC_clean_scaled_UMAP_cluster_cell_type.rds")
