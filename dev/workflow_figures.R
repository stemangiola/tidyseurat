# Article workflow

library(tidyverse)
library(Seurat)
library(SingleR)
library(plotly)
library(tidyHeatmap)
library(ggalluvial)
library(ggplot2)
library(tidyseurat)
options(future.globals.maxSize = 50068 * 1024^2)

# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 9),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )



PBMC_tidy_clean_scaled_UMAP_cluster_cell_type <- readRDS("dev/PBMC_tidy_clean_scaled_UMAP_cluster_cell_type.rds")

p1 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  pivot_longer(
    c(mito.fraction, S.Score, G2M.Score), 
    names_to="property", 
    values_to="Value"
  ) %>%
  mutate(property =  factor(property, levels = c("mito.fraction", "G2M.Score", "S.Score"))) %>%
  ggplot(aes(sample, Value)) + 
  geom_boxplot(outlier.size = 0.5 ) + 
  facet_wrap(~property, scales = "free_y" ) +
  custom_theme +
  theme(aspect.ratio=1)

p2 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  ggplot(aes(UMAP_1, UMAP_2, color=seurat_clusters)) +
  geom_point() +
  custom_theme +
  theme(aspect.ratio=1)

# PBMC_tidy_clean_scaled_UMAP_cluster %>%
#   plot_ly(
#     x = ~`UMAP_1`,
#     y = ~`UMAP_2`,
#     z = ~`UMAP_3`,
#     color = ~seurat_clusters,
#     colors = friendly_cols[1:9],sizes = 50, size = 1
#   )

# Cell_type classification Manual

markers <-
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC)

p3 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  join_transcripts(transcripts=markers$gene %>% head) %>%
  ggplot(aes(seurat_clusters, abundance_SCT, fill=first.labels)) +
  geom_violin() +
  facet_wrap(~transcript) +
  custom_theme

# Plot heatmap
p4 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )

p5 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  join_transcripts(transcripts=markers$gene) %>%
  group_by(seurat_clusters) %>%
  
  # Plot heatmap
  heatmap(
    .row = transcript,
    .column = cell, 
    .value = abundance_SCT, 
    palette_grouping = list(rep("black",9)), 
    palette_value = circlize::colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
  ) %>%
    
  # Add annotation
  add_tile(sample, palette = friendly_cols[1:7]) %>%
  add_point(PC_1) 
  
p6 = 
  PBMC_tidy_clean_scaled_UMAP_cluster_cell_type %>%
  tidyseurat::unite("cluster_cell_type", c(first.labels, seurat_clusters)) %>%
  pivot_longer(
    c(cluster_cell_type, first.labels_single),
    names_to = "classification", values_to = "value"
  ) %>%
  
  ggplot(aes(x = classification, stratum = value, alluvium = cell,
           fill = value, label = value)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  guides(fill = FALSE) +
  scale_fill_manual(values = c(
    "B_cell" = friendly_cols[1],
    "B_cell_4" = friendly_cols[1],
    "T_cells" = friendly_cols[2],
    "T_cells_0" = friendly_cols[2],
    "T_cells_3" = friendly_cols[2],
    "T_cells_5" = friendly_cols[2],
    "T_cells_6" = friendly_cols[2],
    "NK_cell" = friendly_cols[3],
    "NK_cell_1" = friendly_cols[3],
    "NK_cell_7" = friendly_cols[3],
    "Monocyte" = friendly_cols[4],
    "Monocyte_2"  = friendly_cols[4],
    "Monocyte_8"  = friendly_cols[4],
    "Fibroblasts" = friendly_cols[5],
    "GMP" = friendly_cols[6],
    "Macrophage" = friendly_cols[7],
    "Pre-B_cell_CD34+" = friendly_cols[8],
    "Pre-B_cell_CD34-" = friendly_cols[8],
    "HSC_-G-CSF" = friendly_cols[9]
  )) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1),
    text = element_text(size = 9),
    legend.position = "none",
    strip.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
    axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

ggsave("dev/summary_statistics.pdf", p1,  device = "pdf", width = 183/3, height = 50, units = "mm", useDingbats=FALSE)
ggsave("dev/UMAP_2D.pdf", p2,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)
ggsave("dev/violin.pdf", p3,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)
save_pdf(p5, filename = "dev/UMAP_tidyheatmap.pdf", width = 183+50, height = 150, units = "mm")
ggsave("dev/alluvial.pdf", p6,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)


