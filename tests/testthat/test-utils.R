test_that("get_feature_names", {
  data <- data.frame(
    a = c(1:3),
    b = c(4:6),
    c = c(7:9)
  )
  matrix <- matrix(c(1:9), nrow=3, ncol=3)

  expect_equal(get_feature_names(data), c("a", "b", "c"))
  expect_equal(get_feature_names(matrix), c(1:3))
})

test_that("manipulate_matrix", {
  matrix <- matrix(c(1:9), nrow=3, ncol=3)

  expect_equal(manipulate_matrix(matrix), cov(matrix))
  expect_equal(manipulate_matrix(matrix, manipulator="cor"), cor(matrix))
  expect_error(manipulate_matrix(matrix, manipulator="t"))
})

test_that("create_block", {
  feature_names <- c("a", "b", "c", "d", "e")
  feature_names_ <- c(1:5)
  selected_features <- c(2, 3)

  block <- create_block(feature_names, selected_features)
  block_ <- create_block(feature_names_, selected_features)

  expect_s4_class(block, "Block")
  expect_s4_class(block_, "Block")
  expect_equal(block@features, c("b", "c"))
  expect_equal(block@explained_variance, 0)
  expect_equal(block_@features, c(2, 3))
  expect_equal(block_@explained_variance, 0)
})

test_that("scale_eigen_vectors", {
  eigen_vectors <- matrix(c(1:9), nrow=3, ncol=3)
  expected_result <- matrix(c(
      1/3, 2/3, 1,
      4/6, 5/6, 1,
      7/9, 8/9, 1
      ), nrow=3, ncol=3)

  expect_equal(scale_eigen_vectors(eigen_vectors), expected_result)
})

test_that("get_threshold_matrix", {
  matrix <- matrix(c(
      0.1, -0.2, 0.3,
      -0.4, 0.5, 0.6,
      0.7, 0.8, -0.9
  ), nrow=3, ncol=3)
  expected_result <- matrix(c(
      0, 0, 0,
      1, 1, 1,
      1, 1, 1
      ), nrow=3, ncol=3)

  expect_warning(get_threshold_matrix(matrix, threshold=1.3))
  expect_equal(get_threshold_matrix(matrix, threshold=0.3), expected_result)
})

test_that("get_zero_count", {
  matrix <- matrix(c(
    0, 0, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow=3, ncol=3)

  expect_equal(get_zero_count(matrix), c(3, 2, 0))
  expect_error(get_zero_count(c()))
})

test_that("sum_vectors", {
  matrix <- matrix(c(
    1, 0, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow=3, ncol=3)

  expect_equal(sum_vectors(matrix, c(1,3)), c(1, 1, 2))
  expect_equal(sum_vectors(matrix, c(3)), c(0, 1, 1))
  expect_equal(sum_vectors(c(), c(4)), NULL)
  expect_error(sum_vectors(matrix, c(4)))
  expect_error(sum_vectors(matrix(c()), c(4)))
})