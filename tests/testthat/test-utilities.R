context('utilities test')

data("pbmc_small")

test_that("get_special_column_name_symbol", {
  expect_equal(get_special_column_name_symbol(".cell")$symbol, rlang::sym(".cell"))
  expect_equal(get_special_column_name_symbol(".cell")$name, c(".cell"))
})

test_that("ping_old_special_column_into_metadata", {
  ping_old_special_column_into_metadata(pbmc_small) |>
    as_tibble() |>
    colnames() |>
    purrr::pluck(1) |>
    expect_equal("cell")
})

