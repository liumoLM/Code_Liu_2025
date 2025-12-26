library(testthat)
library(ggplot2)

source("code/newnewPlotKoh476Catalog.R")
source("code/Koh89_Koh476_Plotting_Functions.R")

test_that("newPlotKoh476Catalog handles normal input", {
  # Create a sample catalog with 476 values
  set.seed(123)
  catalog <- runif(476, min = 0, max = 100)

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog,
    plot_title = "Test Normal Input",
    num_x_labels = 5
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog handles zero values", {
  # All zeros - this was causing the viewport error
  catalog <- rep(0, 476)

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog,
    plot_title = "Test Zero Values",
    num_x_labels = 0 # No labels when all zeros
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog handles matrix input", {
  set.seed(456)
  catalog_matrix <- matrix(runif(476, 0, 50), ncol = 1)

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog_matrix,
    plot_title = "Test Matrix Input"
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog handles data.frame input", {
  set.seed(789)
  catalog_df <- data.frame(values = runif(476, 0, 50))

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog_df,
    plot_title = "Test DataFrame Input"
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog works with x_axis_label_skip", {
  set.seed(111)
  catalog <- runif(476, 0, 100)

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog,
    plot_title = "Test X-Axis Labels",
    x_axis_label_skip = 20
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog works with num_x_labels = NULL", {
  set.seed(222)
  catalog <- runif(476, 0, 100)

  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog,
    plot_title = "Test No Peak Labels",
    num_x_labels = NULL
  )

  expect_s3_class(p, "ggplot")
})

test_that("newPlotKoh476Catalog can be saved to PDF", {
  set.seed(333)
  catalog <- runif(476, 0, 100)

  temp_pdf <- tempfile(fileext = ".pdf")

  pdf(temp_pdf, width = 14, height = 6)
  p <- newPlotKoh476Catalog(
    Koh476.catalog = catalog,
    plot_title = "Test PDF Output"
  )
  print(p)
  dev.off()

  expect_true(file.exists(temp_pdf))
  expect_gt(file.size(temp_pdf), 0)

  unlink(temp_pdf)
})
