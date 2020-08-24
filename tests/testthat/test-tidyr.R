test_that("nest_unnest",{
  
  col_names = colnames(tt@meta.data) %>% c("cell")
  
  x =     tt %>% nest(data = -Phase) %>% unnest(data) %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA() 
  y =     tt %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()


  expect_equal( 
    x@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1),
    y@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1)
  )
  
  
})
