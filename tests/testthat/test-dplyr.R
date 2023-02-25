context('dplyr test')

library(Seurat)
data("pbmc_small")
set.seed(42)

# test_that("arrange",{
# 
# 
#   pbmc_small_pca_arranged = pbmc_small %>% arrange(nFeature_RNA) %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()
#   pbmc_small_pca = pbmc_small %>% Seurat::ScaleData() %>% Seurat::FindVariableFeatures() %>% Seurat::RunPCA()
# 
#   expect_equal(
#     Seurat::VariableFeatures(pbmc_small_pca_arranged),
#     Seurat::VariableFeatures(pbmc_small_pca)
#   )
# 
#   expect_equal(
#     pbmc_small_pca_arranged[["pca"]]@cell.embeddings ,
#     pbmc_small_pca[["pca"]]@cell.embeddings,
#     tolerance=0.1
#   )
# 
# 
# })

test_that("bind_rows",{

  pbmc_small_bind = bind_rows(    pbmc_small, pbmc_small  )


  expect_equal(
    pbmc_small_bind %>% select(cell) %>% as_tibble() %>% dplyr::count(cell) %>% dplyr::count(n) %>% nrow,
    1
  )



})

test_that("bind_cols",{

  pbmc_small_bind = pbmc_small %>% select(nCount_RNA ,nFeature_RNA)

  pbmc_small %>%
    bind_cols(pbmc_small_bind) %>%
    select(nCount_RNA...2 ,nFeature_RNA...3) %>%
    ncol %>%
  expect_equal( 2)

})

test_that("distinct",{

  expect_equal(   pbmc_small %>% distinct(groups) |> ncol(),    1  )

})

test_that("filter",{

  expect_equal(   pbmc_small %>% filter(groups == "g1") |> ncol(),    44  )

})

test_that("group_by",{

  expect_equal(   pbmc_small %>% group_by(groups) %>% nrow,    80  )

})

test_that("summarise",{

  expect_equal(   pbmc_small %>% summarise(mean(nCount_RNA)) %>% nrow,    1  )

})

test_that("mutate",{

  expect_equal(   pbmc_small %>% mutate(nFeature_RNA = 1) %>% distinct(nFeature_RNA) %>% nrow,    1  )

})

test_that("rename",{

  expect_equal(   pbmc_small %>% rename(s_score = nFeature_RNA) %>% select(s_score) |> ncol(),    1  )

})

test_that("left_join",{

  pbmc_small %>% left_join(pbmc_small %>% distinct(groups) %>% mutate(new_column = 1:2)) %>% `@` (meta.data) |> ncol() %>%
  expect_equal( 8  )

})

test_that("inner_join",{

  expect_equal(   pbmc_small %>% inner_join(pbmc_small %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1)) |> ncol(),    36  )

})

test_that("right_join",{

  expect_equal(   pbmc_small %>% right_join(pbmc_small %>% distinct(groups) %>% mutate(new_column = 1:2) %>% slice(1)) |> ncol(),    36  )

})

test_that("full_join",{

  expect_equal(   pbmc_small %>% full_join(tibble::tibble(groups = "g1", other=1:4)) %>% nrow,    212  )

})

test_that("slice",{

  expect_equal(   pbmc_small %>% slice(1) |> ncol(),    1  )

})

test_that("select",{

  expect_equal(   pbmc_small %>% select(cell, orig.ident ) %>% class %>% as.character,    "Seurat"  )

  expect_equal(   pbmc_small %>% select( orig.ident ) %>% class %>% as.character %>% .[1],    "tbl_df"  )

})

test_that("sample_n",{

  expect_equal(   pbmc_small %>% sample_n(50) |> ncol(),   50  )

  expect_equal(   pbmc_small %>% sample_n(500, replace = TRUE) |> ncol(),   29  )

})

test_that("slice_sample",{
  
  pbmc_small |> 
    slice_sample(n=50) |> 
    ncol() |>
    expect_equal( 50 )
  

})

test_that("sample_frac",{

  expect_equal(   pbmc_small %>% sample_frac(0.1) |> ncol(),   8  )

  expect_equal(   pbmc_small %>% sample_frac(10, replace = TRUE) |> ncol(),   29  )
  
})

test_that("count",{

  pbmc_small |>  
    count(groups) |> 
    nrow() |> 
  expect_equal(     2  )


})

test_that("count",{
  
  pbmc_small |>  
    add_count(groups) |> 
    nrow() |> 
    expect_equal(     230  )
  
  
})