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



PBMC_clean_scaled_UMAP_cluster_cell_type <- readRDS("dev/PBMC_tidy_clean_scaled_UMAP_cluster_cell_type.rds")

p1 = 
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
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
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  sample_n(20000) %>%
  ggplot(aes(UMAP_1, UMAP_2, color=seurat_clusters)) +
  geom_point(size=0.05, alpha=0.2) +
  custom_theme +
  theme(aspect.ratio=1)

PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  sample_n(20000) %>%
  plot_ly(
    x = ~`UMAP_1`,
    y = ~`UMAP_2`,
    z = ~`UMAP_3`,
    color = ~seurat_clusters,
    colors = friendly_cols[1:24],sizes = 50, size = 1
  )

markers = readRDS("dev/PBMC_marker_df.rds")

p3 = 
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  arrange(first.labels) %>%
  mutate(seurat_clusters = fct_inorder(seurat_clusters)) %>%
  join_features(features=c("CD3D", "HLA-DRB1")) %>%
  ggplot(aes(y=seurat_clusters , x=abundance_SCT, fill=first.labels)) +
  geom_density_ridges(bandwidth = 0.2) +
  facet_wrap(~feature, nrow = 2) +
  coord_flip() +
  custom_theme

# Plot heatmap
p4 = 
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  sample_n(2000) %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )

p5 = 
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  sample_n(1000) %>%
  join_features(features=markers$gene) %>%
  mutate(seurat_clusters = as.integer(seurat_clusters)) %>%
  filter(seurat_clusters<10) %>%
  group_by(seurat_clusters) %>%
  
  # Plot heatmap
  heatmap(
    .row = feature,
    .column = cell, 
    .value = abundance_SCT, 
    palette_grouping = list(rep("black",9)), 
    palette_value = circlize::colorRamp2(c(-1.5, 0, 1.5), c("purple", "black", "yellow")),
    
    # ComplexHeatmap parameters
    row_gap = unit(0.1, "mm"), column_gap = unit(0.1, "mm")
  ) %>%
    
  # Add annotation
  add_tile(sample, palette = friendly_cols[1:7]) %>%
  add_point(PC_1) 
  
p6 = 
  PBMC_clean_scaled_UMAP_cluster_cell_type %>%
  tidyseurat::unite("cluster_cell_type", c(first.labels, seurat_clusters), remove=FALSE) %>%
  pivot_longer(
    c(seurat_clusters, first.labels_single),
    names_to = "classification", values_to = "value"
  ) %>%
  
  ggplot(aes(x = classification, stratum = value, alluvium = cell,
           fill = first.labels, label = value)) +
  scale_x_discrete(expand = c(1, 1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  # geom_text(stat = "stratum", size = 3) +
  geom_text_repel(stat = "stratum", size = 3,
             nudge_x      = 0.05,
             direction    = "y",
             angle        = 0,
             vjust        = 0,
             segment.size = 0.2
         ) +
  scale_fill_manual(values = friendly_cols) +
  #guides(fill = FALSE) +
  coord_flip() +
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

ggsave("dev/summary_statistics.pdf", p1,  device = "pdf", width = 183/3, height = 50, units = "mm", useDingbats=FALSE)
ggsave("dev/UMAP_2D.pdf", p2,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)
ggsave("dev/violin.pdf", p3,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)
save_pdf(p5, filename = "dev/UMAPheatmap.pdf", width = 183+50, height = 150, units = "mm")
ggsave("dev/alluvial.pdf", p6,  device = "pdf", width = 89, height = 100, units = "mm", useDingbats=FALSE)


