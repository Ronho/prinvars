test_that("Block", {
  features <- c(2, 3)
  blocks <- list()
  blocks[[1]] <- new("Block", features=features, is_valid=FALSE)
  blocks[[2]] <- new("Block", features=features, explained_variance=0.2)

  expect_s4_class(blocks[[1]], "Block")
  expect_s4_class(blocks[[2]], "Block")
  expect_equal(blocks[[1]]@features, features)
  expect_equal(blocks[[2]]@features, features)
  expect_equal(blocks[[1]]@explained_variance, 0)
  expect_equal(blocks[[2]]@explained_variance, 0.2)
  expect_equal(blocks[[1]]@is_valid, FALSE)
  expect_equal(blocks[[2]]@is_valid, TRUE)
})

test_that("Block - str", {
  features <- c(2, 3)
  block <- new("Block", features=features, explained_variance=0.2)
  result <- "Features (2, 3) explain 20% of the overall explained variance"

  expect_equal(str(block), result)
})