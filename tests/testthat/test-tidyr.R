context('tidyr test')

tt = pbmc_small@assays$RNA@counts %>% CreateSeuratObject() %>% tidy %>% mutate(groups = sprintf("g%s", rep(1:2, dplyr::n()/2)))

test_that("nest_unnest",{
  
  col_names = colnames(tt@meta.data) %>% c("cell")
  
  x =     tt %>% nest(data = -groups) %>% unnest(data) %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA() 
  y =     tt %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()


  expect_equal( 
    x@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1),
    y@reductions$pca@cell.embeddings %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC_1)
  )
  
  
})

test_that("unite separate",{
  
  un = tt %>% unite("new_col", c(orig.ident, groups)) 

  expect_equal( un %>% select(new_col) %>% slice(1) %>% pull(new_col),   "SeuratProject_g1")
  
  se = un %>% separate(col = new_col, into= c("orig.ident", "groups"))
  
  expect_equal( se %>% select(orig.ident) %>% ncol,   1)
  
  
})

test_that("extract",{
  
  expect_equal( 
    tt %>% extract(groups, into = "g", regex = "g([0-9])", convert = TRUE) %>% pull(g) %>% class , 
    "integer"
  )
  
  
})

test_that("pivot_longer",{
  
  expect_equal( 
    tt %>% pivot_longer(c(orig.ident, groups), names_to = "name", values_to = "value")  %>% class %>% .[1], 
    "tbl_df"
  )
  
  
})