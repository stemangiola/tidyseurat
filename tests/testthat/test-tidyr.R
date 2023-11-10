context('tidyr test')

data("pbmc_small")
tt <- GetAssayData(pbmc_small, layer = 'counts', assay = "RNA") |> CreateSeuratObject() |> mutate(groups = sprintf("g%s", rep(1:2, dplyr::n()/2)))

test_that("nest_unnest", {
  col_names <- colnames(tt[[]]) |> c("cell")
  x <- tt |> nest(data = -groups) |> unnest(data) |> Seurat::NormalizeData() |>  Seurat::ScaleData() |> Seurat::FindVariableFeatures() |> Seurat::RunPCA()
  y <- tt |> Seurat::NormalizeData() |> Seurat::ScaleData() |> Seurat::FindVariableFeatures() |> Seurat::RunPCA()
  expect_equal(
    x[["pca"]]@cell.embeddings |> as_tibble(rownames = "cell") |> arrange(cell) |> pull(PC_1),
    y[["pca"]]@cell.embeddings |> as_tibble(rownames = "cell") |> arrange(cell) |> pull(PC_1)
  )
})

test_that("fast_vs_slow_nest", {
  expect_identical(
    tt |> mutate(groups2 = groups) |> nest(data = -c(groups, groups2)) |> select(-groups2),
    tt |> nest(data = -groups) 
  )
})

test_that("nest_unnest_slice_1", {
  expect_equal(
    tt |> nest(data = -groups) |> slice(1) |> unnest(data) |> ncol(),
    sum(tt[[]]$groups == "g1")
  )
})
test_that("unite separate", {
  un <- tt |> unite("new_col", c(orig.ident, groups))
  se <- un |> separate(col = new_col, into = c("orig.ident", "groups"))
  expect_equal(un |> select(new_col) |> slice(1) |> pull(new_col), "SeuratProject_g1")
  expect_equal(se |> select(orig.ident) |> ncol(), 1)
})

test_that("extract", {
  expect_equal(
    tt |> extract(groups, into = "g", regex = "g([0-9])", convert = TRUE) |> pull(g) |> class(),
    "integer"
  )
})

test_that("pivot_longer", {
  expect_equal(
    tt |> pivot_longer(c(orig.ident, groups), names_to = "name", values_to = "value")  |> class() |> magrittr::extract2(1),
    "tbl_df"
  )
})
