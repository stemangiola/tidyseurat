context('methods test')

data("pbmc_small")

test_that("join_features_long", {
  pbmc_small |> 
    join_features("CD3D", shape="long") |> 
    slice(1) |>
    pull(.abundance_RNA) |>
    expect_equal(6.35, tolerance = 0.1)
})

test_that("join_features_wide", {
  pbmc_small |> 
    join_features("CD3D", shape="wide") |> 
    slice(1) |>
    pull(CD3D) |>
    expect_equal(6.35, tolerance = 0.1)
})

test_that("aggregate_cells() returns expected values", {
  # Create pseudo-bulk object for testing
  pbmc_pseudo_bulk <-
    pbmc_small |>
    aggregate_cells(c(groups, letter.idents), assays = "RNA")
  
  # Check row length is unchanged
  pbmc_pseudo_bulk |>
    distinct(.feature) |> 
    nrow() |>
    expect_equal(pbmc_small |> nrow())
  
  # Check column length is correctly modified
  pbmc_pseudo_bulk |> 
    distinct(.sample) |> 
    nrow() |>
    expect_equal(pbmc_small |>
                   as_tibble() |>
                   select(groups, letter.idents) |>
                   unique() |>
                   nrow()
    )
  # Spot check for correctly aggregated count value of ACAP1 gene
  pbmc_pseudo_bulk |> 
    filter(.feature == "ACAP1" & .sample == "g1___A") |> 
    select(RNA) |> 
    as.numeric() |> 
    expect_equal(
      Seurat::DietSeurat(pbmc_small, assays = "RNA", features = "ACAP1")[, pbmc_small |>
                                  as_tibble() |>
                                  filter(groups == "g1", letter.idents == "A") |>
                                  pull(.cell)] |>
        LayerData() |> 
        sum())
})

test_that("get_abundance_sc_wide", {
  expect_equal(
    pbmc_small |> get_abundance_sc_wide() |> nrow(),
    pbmc_small[[]] |> nrow()
  )
  expect_equal(
    pbmc_small |> get_abundance_sc_wide() |> pull("S100A9") |> sum(),
    pbmc_small |> FetchData("S100A9") |> sum(), 
    tolerance = 0.1
 )
})

test_that("get_abundance_sc_long", {
  expect_equal(pbmc_small |> get_abundance_sc_long() |> ncol(), 3)
  expect_equal(
    pbmc_small |> get_abundance_sc_long() |> filter(.feature == "S100A9") |> pull(".abundance_RNA") |> sum(),
    pbmc_small |> FetchData("S100A9") |> sum(),
    tolerance = 0.1
  )
})



