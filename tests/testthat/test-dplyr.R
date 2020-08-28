context('dplyr test')

tt = pbmc_small %>% tidy

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

test_that("bind_cols",{
  
  tt_bind = tt %>% select(nCount_RNA ,nFeature_RNA)
  
  
  expect_equal(
    tt %>% bind_cols(tt_bind) %>% select(nCount_RNA...9 ,nFeature_RNA...10) %>% ncol,
    2
  )
  
})

test_that("distinct",{
  
  expect_equal(   tt %>% distinct(groups) %>% ncol,    1  )
  
})

test_that("filter",{
  
  expect_equal(   tt %>% filter(groups == "g1") %>% ncol,    44  )
  
})

test_that("group_by",{
  
  expect_equal(   tt %>% group_by(groups) %>% nrow,    80  )
  
})

test_that("summarise",{
  
  expect_equal(   tt %>% summarise(mean(nCount_RNA)) %>% nrow,    1  )
  
})

test_that("mutate",{
  
  expect_equal(   tt %>% mutate(nFeature_RNA = 1) %>% distinct(nFeature_RNA) %>% nrow,    1  )
  
})

test_that("rename",{
  
  expect_equal(   tt %>% rename(s_score = nFeature_RNA) %>% select(s_score) %>% ncol,    1  )
  
})

test_that("left_join",{
  
  expect_equal(   tt %>% left_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2)) %>% `@` (meta.data) %>% ncol,    9  )
  
})

test_that("inner_join",{
  
  expect_equal(   tt %>% inner_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1)) %>% ncol,    36  )
  
})

test_that("right_join",{
  
  expect_equal(   tt %>% right_join(tt %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1)) %>% ncol,    36  )
  
})

test_that("full_join",{
  
  expect_equal(   tt %>% full_join(tibble::tibble(groups = "g1", other=1:4)) %>% nrow,    212  )
  
})

test_that("slice",{
  
  expect_equal(   tt %>% slice(1) %>% ncol,    1  )
  
})

test_that("select",{
  
  expect_equal(   tt %>% select(cell, orig.ident ) %>% class %>% as.character,    "tidyseurat"  )

  expect_equal(   tt %>% select( orig.ident ) %>% class %>% as.character %>% .[1],    "tbl_df"  )
  
})

test_that("sample_n",{
  
  expect_equal(   tt %>% sample_n(50) %>% ncol,   50  )
  

})

test_that("sample_frac",{
  
  expect_equal(   tt %>% sample_frac(0.1) %>% ncol,   8  )
  
  
})

test_that("count",{
  
  expect_equal(   tt %>% count(groups) %>% nrow,   2  )
  
  
})