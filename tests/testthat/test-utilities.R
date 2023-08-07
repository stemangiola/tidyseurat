context('utilities test')

data("pbmc_small")

test_that("get_abundance_sc_wide", {
  pbmc_abundance_wide <-
    pbmc_small |>
    get_abundance_sc_wide()
  pbmc_abundance_wide |> 
    nrow() |>
    expect_equal(80)
  pbmc_abundance_wide |> 
    ncol() |>
    expect_equal(21)
  pbmc_abundance_wide |> 
    pull("S100A9") |>
    sum() |>
    expect_equal(175.2576, tolerance = 0.1)
})

test_that("get_abundance_sc_long", {
  pbmc_abundance_long <-
    pbmc_small |>
    get_abundance_sc_long()
  pbmc_abundance_long |> 
    nrow() |>
    expect_equal(1600)
  pbmc_abundance_long |> 
    ncol() |>
    expect_equal(3)
  pbmc_abundance_long |> 
    filter(.feature == "S100A9") |>
    pull(".abundance_RNA") |>
    sum() |>
    expect_equal(175.2576, tolerance = 0.1)
})