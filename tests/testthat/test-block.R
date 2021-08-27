test_that("Block", {
  features <- c(2, 3)
  block <- new("Block", features=features)
  block_ <- new("Block", features=features, explained_variance=0.2)

  expect_s4_class(block, "Block")
  expect_s4_class(block_, "Block")
  expect_equal(block@features, features)
  expect_equal(block_@features, features)
  expect_equal(block@explained_variance, 0)
  expect_equal(block_@explained_variance, 0.2)
})

test_that("Block - str", {
  features <- c(2, 3)
  block <- new("Block", features=features, explained_variance=0.2)
  result <- "Features (2, 3) explain 20% of the overall explained variance"

  expect_equal(str(block), result)
})