library(tidyverse)
library(glue)

tibble(
  observation = glue("observation {1:100}"),
  variable_1 = rep("...", 100),
  variable_2 = rep("...", 100),
  variable_3 = map(1:100, ~ tibble(a = 1:10, b = 1:10)), 
  variable_4 = map(1:100, ~ ggplot()),
  variable_5 = map(1:100, ~ lm(y ~ x, data = data.frame(x=1:10, y=1:10))),
  variable_6 = map(1:100, ~ pbmc_small),
  variable_7 = map(1:100, ~ tidySingleCellExperiment::pbmc_small)
)


my_vector = seq(1, 20); 

# Imperative
my_vector_modified = c()
for(i in 1:length(my_vector)) {
  my_vector_modified[i] = my_vector[i] * 2L
}
# Functional
my_vector_modified = my_vector |> map_int(~ .x * 2L)

df = data.frame(a= rep("a", ncol(SeuratObject::pbmc_small)), b= rep("b", ncol(SeuratObject::pbmc_small)))
rownames(df)  = colnames(SeuratObject::pbmc_small)

info = rep(1, ncol(SeuratObject::pbmc_small))
SeuratObject::pbmc_small |>
  AddMetaData(info, "info")

colData(tidySingleCellExperiment::pbmc_small) |> cbind()


# Subsampling

single_cell_data |>
  add_count(sample, name = "tot_cells") |>
  mutate(median_cells = min(tot_cells)) |>
  nest(data = -c(sample, median_cells)) |>
  mutate(data = map2(data, median_cells, ~ sample_n(.x, .y, replace = TRUE))) |>
  unnest(data)

# Define cell categories for analysis plotting

single_cell_data |>
  
  mutate(cell_differentiation = 
           case_when(
             curated_cell_type_pretty %in% c("B immature", "B mem") ~ "B",
             curated_cell_type_pretty %in% c("pDC") ~ "pDC",
             cell_differentiation == "lymphoid" ~ "T+NK",
             cell_differentiation == "myeloid" ~ "Myeloid"
           )
  ) |> 
  
  mutate(
    curated_cell_type_pretty = if_else(
      curated_cell_type_pretty %in% c("T gd1",  "T gd2"), 
      "gamma_delta" , 
      curated_cell_type_pretty
    )
  ) 


# Quality control




# Gating gamma delta
seurat_obj_sig = seurat_obj |>
  
  
  join_features(
    features = c("CD3D", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B"),
    shape = "wide",
    assay = "SCT"
    
  ) |>
  
  mutate(signature_score =
           scales::rescale(CD3D + TRDC + TRGC1+ TRGC2, to=c(0,1)) -
           scales::rescale(CD8A + CD8B, to=c(0,1))
  ) |>
  
  Seurat::FeaturePlot(signature_score) |>
  mutate( gate = tidygate::gate_int(UMAP_1, UMAP_2) ) |> 
  
  filter(gate == 1) %>%
  
  NormalizeData() |> 
  FindVariableFeatures( nfeatures = 100)

  split_group(sample) %>% 
  RunFastMNN() |> 
  RunUMAP(reduction = "mnn", dims = 1:20) |> 
  FindNeighbors( dims = 1:20, reduction = "mnn") |> 
  FindClusters( resolution = 0.3) |>


# gamma_delta_df = 
#   readRDS("cancer_only_analyses/integrated_counts_curated.rds")  |> 
#   #	{.x = (.); DefaultAssay(.x) = "RNA"; .x} |> 
#   filter(curated_cell_type_pretty %in% c("T gd1",  "T gd2")) |> 
#   
#   {
#     .x= (.)
#     DefaultAssay(.x) = "RNA"
#     .x[["SCT"]] = NULL
#     .x[["integrated"]] = NULL
#     .x
#   } |> 
#   NormalizeData() |> 
#   FindVariableFeatures( nfeatures = 100) |> 
#   mutate(batch_to_eliminate = sample) |> 
#   nest(data = -batch_to_eliminate) |>
#   pull(data) |> 
#   RunFastMNN() |> 
#   RunUMAP(reduction = "mnn", dims = 1:20) |> 
#   FindNeighbors( dims = 1:20, reduction = "mnn") |> 
#   FindClusters( resolution = 0.3) |>
#   mutate(gate = tidygate::gate_int(UMAP_1, UMAP_2, how_many_gates = 2, gate_list = readRDS("file66175abbca44.rds"))) |> 
#   tidysc::adjust_abundance(~ 1) |>
#   mutate(gamma_delta = case_when(
#     gate == 0 ~ "T gd vd2",
#     gate == 1 ~ "T gd vd1 LGALS1",
#     gate == 2 ~ "T gd vd1",
#   ))