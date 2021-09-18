test_that("check_pla_equality", {

})

test_that("equal_block_elements", {
    
})

test_that("equal_features", {
    
})

test_that("equal_explained_variance", {
    
})

test_that("get_indices", {
    
})

test_that("conditional_matrix", {
    
})


test_that("get_feature_names", {
  data <- list()
  data[[1]] <- data.frame(
    a = c(1:3),
    b = c(4:6),
    c = c(7:9)
  )
  data[[2]] <- matrix(c(1:9), nrow=3, ncol=3)

  expect_equal(get_feature_names(x=data[[1]]), c("a", "b", "c"))
  expect_equal(get_feature_names(x=data[[2]]), c(1:3))
})

test_that("create_block", {
  selected_features <- c(2, 3)

  blocks <- list()
  blocks[[1]] <- create_block(
    feature_names=c("a", "b", "c", "d", "e"),
    selected_features=selected_features
  )
  blocks[[2]] <- create_block(
    feature_names=c(1:5),
    selected_features=selected_features
  )

  expect_s4_class(blocks[[1]], "Block")
  expect_s4_class(blocks[[2]], "Block")
  expect_equal(blocks[[1]]@features, c("b", "c"))
  expect_equal(blocks[[1]]@explained_variance, 0)
  expect_equal(blocks[[2]]@features, c(2, 3))
  expect_equal(blocks[[2]]@explained_variance, 0)
})

test_that("get_zero_count", {
  matrix <- matrix(c(
    0, 0, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow=3, ncol=3)

  expect_equal(get_zero_count(matrix), c(2, 2, 1))
  expect_error(get_zero_count(c()))
})