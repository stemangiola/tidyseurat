context('pillar test')

test_string <- "A small string to test the function of pillar utilities."

test_that("pillar___format_comment", {
  test_string |>
    pillar___format_comment(width = 20) |>
    stringr::str_count("# ") |> 
    expect_equal(5)
})

test_that("pillar___strwrap2", {
  test_string |>
    pillar___strwrap2(width = 20, indent = 4) |>
    stringr::str_count("      ") |> 
    expect_equal(c(0, 1, 1, 1, 1))
})

test_that("pillar___wrap", {
  test_string |>
    pillar___wrap(width = 20) |>
    stringr::str_count("\n") |> 
    expect_equal(3)
})