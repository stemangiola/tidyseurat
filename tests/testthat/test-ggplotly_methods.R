context('ggplot test')

data("pbmc_small")

df <- pbmc_small
df$number <- rnorm(ncol(df))
df$factor <- sample(gl(3, 1, ncol(df)))

test_that("ggplot", {
  # cell metadata
  p <- ggplot(df, aes(factor, number)) 
  expect_silent(show(p))
  expect_s3_class(p, "ggplot")
  # assay data
  g <- sample(rownames(df), 1)
  fd <- join_features(df, g, shape="wide")
  p <- ggplot(fd, aes(factor, .data[[g]]))
  expect_silent(show(p))
  expect_s3_class(p, "ggplot")
  # reduced dimensions
  p <- ggplot(df, aes(PC_1, PC_2, col=factor))
  expect_silent(show(p))
  expect_s3_class(p, "ggplot")
})

test_that("plotly", {
  # cell metadata
  p <- plot_ly(df, x=~factor, y=~number, type="violin") 
  expect_silent(show(p))
  expect_s3_class(p, "plotly")
  # assay data
  g <- sample(rownames(df), 1)
  fd <- join_features(df, g, shape="wide")
  p <- plot_ly(fd, x=~factor, y=g, type="violin") 
  expect_silent(show(p))
  expect_s3_class(p, "plotly")
  # reduced dimensions
  p <- plot_ly(fd, x=~PC_1, y=~PC_2, type="scatter", mode="markers") 
  expect_silent(show(p))
  expect_s3_class(p, "plotly")
})
