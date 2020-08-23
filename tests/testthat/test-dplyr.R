context('dplyr test')



test_that("arrange",{


  tt_pca_aranged = tt %>% arrange(nFeature_RNA) %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()
  tt_pca = tt %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()

  expect_equal(
    Seurat::VariableFeatures(tt_pca_aranged),
    Seurat::VariableFeatures(tt_pca)
  )
  
  expect_equal(
    tt_pca_aranged@reductions$pca@cell.embeddings ,
    tt_pca@reductions$pca@cell.embeddings
  )


})

test_that("bind_rows",{
  
  tt_bind = bind_rows(    tt, tt  )
  

  expect_equal(
    tt_bind %>% select(cell) %>% tidyseurat:::to_tib() %>% dplyr::count(cell) %>% dplyr::count(n) %>% nrow,
    1
  )
  

  
})