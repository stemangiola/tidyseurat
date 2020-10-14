# Article workflow

library(tidyverse)
library(Seurat)
library(SingleR)
library(plotly)
# library(future)
# plan(multisession, workers=10)
options(future.globals.maxSize = 50068 * 1024^2)
library(tidyseurat)
friendly_cols <- dittoSeq::dittoColors()

# PBMC = PBMC %>% 
#   select(1:11, -old.ident) %>%
#   mutate(sample = sprintf("./data/seurat/outs/%s", sample)) %>% 
#   mutate(Phase = sprintf("phase_%s", Phase))

PBMC_tidy <- readRDS("dev/PBMC.rds")

PBMC_tidy_clean_scaled <-
  PBMC_tidy %>%
  mutate(grouping = sample) %>%
  nest(sample_df = -grouping) %>%
  mutate(sample_df = map( sample_df,~  SCTransform(.x)))

my_features =
  PBMC_tidy_clean_scaled$sample_df %>%
  SelectIntegrationFeatures(nfeatures = 2000)

PBMC_integrated = 
  
  PBMC_tidy_clean_scaled$sample_df %>%
  PrepSCTIntegration(anchor.features = my_features) %>%
  FindIntegrationAnchors(
    normalization.method = "SCT",
    anchor.features = my_features
  ) %>%
  IntegrateData(normalization.method = "SCT")

PBMC_integrated %>% saveRDS("dev/PBMC_integrated.rds")