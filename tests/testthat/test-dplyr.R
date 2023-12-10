context('dplyr test')

library(Seurat)
data("pbmc_small")
set.seed(42)

test_that("arrange", {
  
  pbmc_small |> 
    arrange(nFeature_RNA) |> 
    expect_warning(regexp = "`arrange\\(\\)` was deprecated in tidyseurat .*")
  
  # pbmc_small_pca_arranged <- pbmc_small |> arrange(nFeature_RNA) |> Seurat::ScaleData() |> Seurat::FindVariableFeatures() |> Seurat::RunPCA()
  # pbmc_small_pca <- pbmc_small |> Seurat::ScaleData() |> Seurat::FindVariableFeatures() |> Seurat::RunPCA()
  # expect_equal(
  #   Seurat::VariableFeatures(pbmc_small_pca_arranged),
  #   Seurat::VariableFeatures(pbmc_small_pca)
  # )
  
  # # Failing only for ATLAS CRAN, but succeding for the rest
  # expect_equal(
  #   pbmc_small_pca_arranged[["pca"]]@cell.embeddings,
  #   pbmc_small_pca[["pca"]]@cell.embeddings,
  #   tolerance=0.1
  # )
  # expect_equal(
  #   pbmc_small_pca_arranged |> as_tibble() |>dplyr::slice_head(n = 1),
  #   pbmc_small_pca |> as_tibble() |> dplyr::slice_min(nFeature_RNA, n = 1)
  # )

})

test_that("bind_cols", {
  pbmc_small_bind <- pbmc_small |> select(nCount_RNA, nFeature_RNA)
  pbmc_small |>
    ttservice::bind_cols(pbmc_small_bind) |>
    select(nCount_RNA...2, nFeature_RNA...3) |>
    ncol() |>
    expect_equal(2)
})

test_that("distinct", {
  expect_equal(pbmc_small |> distinct(groups) |> ncol(), 1)
})

test_that("filter", {
  expect_equal(
    pbmc_small |> filter(groups == "g1") |> ncol(),
    sum(pbmc_small[[]]$groups == "g1")
  )
})

test_that("group_by", {
  expect_equal(
    pbmc_small |> group_by(groups) |> nrow(),
    nrow(pbmc_small[[]])
  )
})

test_that("summarise", {
  expect_equal(pbmc_small |> summarise(mean(nCount_RNA)) |> nrow(), 1)
})

test_that("mutate", {
  expect_equal(pbmc_small |> mutate(nFeature_RNA = 1) |> distinct(nFeature_RNA) |> nrow(), 1)
})

test_that("rename", {
  expect_equal(pbmc_small |> rename(s_score = nFeature_RNA) |> select(s_score) |> ncol(), 1)
})

test_that("left_join", {
  expect_equal(
    pbmc_small |> left_join(pbmc_small |> distinct(groups) |> mutate(new_column = 1:2) |> slice(1)) |> ncol(),
    nrow(pbmc_small[[]])
  )
})

test_that("inner_join", {
  expect_equal(
    pbmc_small |> inner_join(pbmc_small |> distinct(groups) |> mutate(new_column = 1:2) |> slice(1)) |> ncol(),
    sum(pbmc_small[[]]$groups == "g2")
  )
})

test_that("right_join", {
  expect_equal(
    pbmc_small |> right_join(pbmc_small |> distinct(groups) |> mutate(new_column = 1:2) |> slice(1)) |> ncol(),
    sum(pbmc_small[[]]$groups == "g2")
  )
})

test_that("full_join", {
  expect_equal(
    pbmc_small |> full_join(tibble::tibble(groups = "g1", other = 1:4)) |> nrow(),
    sum(pbmc_small[[]]$groups == "g1") * 4 + sum(pbmc_small[[]]$groups == "g2")
  )
})

test_that("slice", {
  expect_equal(pbmc_small |> slice(1) |> ncol(), 1)
  expect_equal(
    pbmc_small |> slice(1:6) |> colnames(),
    colnames(pbmc_small) |> head(6))
})

test_that("sample_n", {
  expect_equal(pbmc_small |> sample_n(50) |> ncol(), 50)
  expect_equal(
    pbmc_small |> sample_n(500, replace = TRUE) |> ncol(),
    pbmc_small |> as_tibble() |> ncol()
  )
})

test_that("slice_sample", {
  pbmc_small |>
    slice_sample(n = 50) |>
    ncol() |>
    expect_equal(50)
})

test_that("slice_head", {
  pbmc_small |>
    slice_head(n = 50) |>
    ncol() |>
    expect_equal(50)
  expect_equal(
    colnames(pbmc_small) |> head(n = 50),
    pbmc_small |> slice_head(n = 50) |> colnames()
  )
})

test_that("slice_tail", {
  pbmc_small |>
    slice_tail(n = 50) |>
    ncol() |>
    expect_equal(50)
  expect_equal(
    colnames(pbmc_small) |> tail(n = 50),
    pbmc_small |> slice_tail(n = 50) |> colnames()
  )
})

test_that("slice_min", {
  pbmc_small |>
    slice_min(nFeature_RNA, n = 5) |>
    ncol() |>
    expect_equal(5)
  
  # Arrange is deprecated
  # expect_equal(
  #   pbmc_small |> as_tibble() |> arrange(nFeature_RNA) |> head(n = 5) %>% pull(.cell),
  #   pbmc_small |> slice_min(nFeature_RNA, n = 5) |> colnames()
  # )
})

test_that("slice_max", {
  pbmc_small |>
    slice_max(nFeature_RNA, n = 5) |>
    ncol() |>
    expect_equal(5)
  
  # Arrange is deprecated
  # expect_equal(
  #   pbmc_small |> as_tibble() |> arrange(desc(nFeature_RNA)) |> head(n = 5) %>% pull(.cell),
  #   pbmc_small |> slice_max(nFeature_RNA, n = 5) |> colnames()
  # )
})

test_that("slice_min slice_max tibble input for order_by", {
  pbmc_small |>
    slice_min(tibble::tibble(nFeature_RNA, nCount_RNA), n = 5) |>
    ncol() |>
    expect_equal(5)
  pbmc_small |>
    slice_max(tibble::tibble(nFeature_RNA, nCount_RNA), n = 5) |>
    ncol() |>
    expect_equal(5)
})

test_that("select", {
  expect_equal(pbmc_small |> select(cell, orig.ident) |> class() |> as.character(), "Seurat")
  expect_equal(pbmc_small |> select(orig.ident) |> class() |> as.character() |> purrr::pluck(1), "tbl_df")
})

test_that("sample_frac", {
  expect_equal(
    pbmc_small |> sample_frac(0.1) |> ncol(),
    nrow(pbmc_small[[]]) * 0.1
  )
  expect_equal(
    pbmc_small |> sample_frac(10, replace = TRUE) |> ncol(),
    pbmc_small |> as_tibble() |> ncol()
  )
})

test_that("count", {
  expect_equal(
    pbmc_small |> count(groups) |> nrow(),
    pbmc_small[[]]$groups |> unique() |> length()
  )
})

test_that("add_count", {
  expect_equal(
    pbmc_small |> add_count(groups) |> nrow(),
    pbmc_small |> rownames() |> length()
  )
})

test_that("rowwise", {
  expect_equal(
    pbmc_small |> rowwise() |> mutate(m = mean(c(nCount_RNA, nFeature_RNA))) |> purrr::pluck("m", 1),
    ((pbmc_small[, 1]$nCount_RNA + pbmc_small[, 1]$nFeature_RNA) / 2) |> unname()
  )
})

test_that("group_split() works for one variable", {
  fd <- pbmc_small |> 
    group_split(groups)
  expect_equal(length(fd), length(unique(pbmc_small$groups)))
})

test_that("group_split() works for combination of variables", {
  fd <- pbmc_small |> 
    group_split(groups, letter.idents)
  expect_equal(length(fd), length(unique(pbmc_small$groups)) *
                 length(unique(pbmc_small$letter.idents)))
})

test_that("group_split() works for one logical statement", {
  fd_log <- pbmc_small |> 
    group_split(groups=="g1")
  fd_var <- pbmc_small |> 
    group_split(groups=="g1")
  expect_equal(lapply(fd_var, count), lapply(fd_log, count))
})

test_that("group_split() works for two logical statements", {
  fd <- pbmc_small |>
    group_split(PC_1>0 & groups=="g1")
  fd_counts <- lapply(fd, count)
  expect_equal(c(fd_counts[[1]], fd_counts[[2]], use.names = FALSE), 
               list(75, 5))
})
