context('print test')

data("pbmc_small")

test_that("print", {
  text <- capture.output(print(pbmc_small))
  expect_equal(grep("Seurat-tibble abstraction", text), 1)
  i <- grep(str <- ".*Features=([0-9]+).*", text)
  expect_equal(gsub(str, "\\1", text[i]), paste(nrow(pbmc_small)))
  i <- grep(str <- ".*Cells=([0-9]+).*", text)
  expect_equal(gsub(str, "\\1", text[i]), paste(ncol(pbmc_small)))
})

test_that("glimpse", {
  text <- capture.output(glimpse(pbmc_small))
  expect_equal(length(text), 37)
})
