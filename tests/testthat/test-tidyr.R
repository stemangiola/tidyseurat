context('tidyr test')

tt = pbmc_small@assays$RNA@counts %>% CreateSeuratObject() %>% tidy %>% mutate(groups = rep(1:2, dplyr::n()/2))

test_that("nest_unnest",{
  
  col_names = colnames(tt@meta.data) %>% c("cell")
  
  x =     tt %>% nest(data = -groups) %>% unnest(data) %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA() 
  y =     tt %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()


  expect_equal( 
    x@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1),
    y@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1)
  )
  
  
})
