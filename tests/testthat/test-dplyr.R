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

test_that("bind_cols",{
  
  tt_bind = tt %>% select(batch, low_quality)
  
  
  expect_equal(
    tt %>% bind_cols(tt_bind) %>% select(batch...18, low_quality...19) %>% ncol,
    2
  )
  
})

test_that("distinct",{
  
  expect_equal(   tt %>% distinct(Phase) %>% ncol,    1  )
  
})

test_that("filter",{
  
  expect_equal(   tt %>% filter(Phase == "G2M") %>% ncol,    358  )
  
})

test_that("group_by",{
  
  expect_equal(   tt %>% group_by(Phase) %>% nrow,    1000  )
  
})

test_that("summarise",{
  
  expect_equal(   tt %>% summarise(mean(mito.tot)) %>% nrow,    1  )
  
})

test_that("mutate",{
  
  expect_equal(   tt %>% mutate(S.Score = 1) %>% distinct(S.Score) %>% nrow,    1  )
  
})

test_that("rename",{
  
  expect_equal(   tt %>% rename(s_score = S.Score) %>% select(s_score) %>% ncol,    1  )
  
})

test_that("left_join",{
  
  expect_equal(   tt %>% left_join(tt %>% distinct(Phase) %>% mutate(new_column = 1:3)) %>% `@` (meta.data) %>% ncol,    18  )
  
})

test_that("inner_join",{
  
  expect_equal(   tt %>% inner_join(tt %>% distinct(Phase) %>% mutate(new_column = 1:3) %>% slice(1)) %>% ncol,    358  )
  
})

test_that("right_join",{
  
  expect_equal(   tt %>% right_join(tt %>% distinct(Phase) %>% mutate(new_column = 1:3) %>% slice(1)) %>% ncol,    358  )
  
})

test_that("full_join",{
  
  expect_equal(   tt %>% full_join(tibble::tibble(Phase = "G2M", other=1:4)) %>% nrow,    2074  )
  
})

test_that("slice",{
  
  expect_equal(   tt %>% slice(1) %>% ncol,    1  )
  
})

test_that("select",{
  
  expect_equal(   tt %>% select(cell, S.Score ) %>% class %>% as.character,    "tidyseurat"  )

  expect_equal(   tt %>% select( S.Score ) %>% class %>% as.character %>% .[1],    "tbl_df"  )
  
})



