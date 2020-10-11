# Article workflow

library(tidyverse)
library(Seurat)
library(SingleR)
library(plotly)
library(tidyseurat)
friendly_cols <- dittoSeq::dittoColors()

# PBMC = PBMC %>% 
#   select(1:11, -old.ident) %>%
#   mutate(sample = sprintf("./data/seurat/outs/%s", sample)) %>% 
#   mutate(Phase = sprintf("phase_%s", Phase))

PBMC_tidy <- readRDS("dev/PBMC.rds")

# Polishing

PBMC_tidy_clean <-
  PBMC_tidy %>%
  
  # Clean groups
  mutate(Phase = Phase %>% str_remove("^phase_")) %>%
  
  # Extract sample
  extract(sample, "sample", "./data/seurat/outs/([a-zA-Z0-9]+)")

# PBMC_tidy_clean = PBMC_tidy_clean %>% nest(data = -sample) %>% mutate(data = map(data, ~ .x %>% sample_n(50))) %>% unnest(data)

# Scaling

PBMC_tidy_clean_scaled <-
  PBMC_tidy_clean %>%
  SCTransform(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE)

# Dimensionality reduction
PBMC_tidy_clean_scaled_UMAP <-
  PBMC_tidy_clean_scaled %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, n.components = 3L)

# Clustering

PBMC_tidy_clean_scaled_UMAP_cluster <-
  PBMC_tidy_clean_scaled_UMAP %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(method = "igraph", verbose = FALSE)

# PBMC_tidy_clean_scaled_UMAP_cluster %>%
#   plot_ly( 
#     x = ~`UMAP_1`,
#     y = ~`UMAP_2`,
#     z = ~`UMAP_3`,
#     color = ~seurat_clusters,
#     colors = "Dark2"
#   )

# Cell_type classification Manual

markers <-
  PBMC_tidy_clean_scaled_UMAP_cluster %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC)

# Plot heatmap
PBMC_tidy_clean_scaled_UMAP_cluster %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )

ggsave("dev/heatmap.pdf", units = "mm", width = 89, height = 100)

# Cell_type classification Automatic

# Get cell type reference data
hpca <- HumanPrimaryCellAtlasData()

# Infer cell identities
cell_type_df <-
  
  # extracting counts from Seurat object
  PBMC_tidy_clean_scaled_UMAP_cluster@assays[["SCT"]]@counts %>%
  log1p() %>%
  # SingleR
  SingleR(
    ref = hpca,
    labels = hpca$label.main,
    method = "cluster",
    clusters = PBMC_tidy_clean_scaled_UMAP_cluster %>% pull(seurat_clusters)
  ) %>%
  
  # Formatting results
  as.data.frame() %>%
  as_tibble(rownames = "seurat_clusters") %>%
  select(seurat_clusters, first.labels)

# Join UMAP and cell type info
PBMC_tidy_clean_scaled_UMAP_cluster_cell_type <-
  PBMC_tidy_clean_scaled_UMAP_cluster %>%
  left_join(
    cell_type_df, 
    by = "seurat_clusters"
  )
# Reorder columns
PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  count(seurat_clusters, first.labels)

saveRDS(PBMC_tidy_clean_scaled_UMAP_cluster_cell_type, "dev/PBMC_tidy_clean_scaled_UMAP_cluster_cell_type.rds")
