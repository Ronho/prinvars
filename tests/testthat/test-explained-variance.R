test_that("proportional_explained_variance", {
  eigen_vectors <- matrix(c(1:9), nrow=3, ncol=3)
  eigen_values <- c(0.5, 0.3, 0.2)
  variables <- c(1, 2)

  expect_equal(proportional_explained_variance(eigen_values, eigen_vectors, variables, type="exact"), 73.8)
  expect_equal(proportional_explained_variance(eigen_values, eigen_vectors, variables, type="approx"), 0.8)
  expect_equal(proportional_explained_variance(eigen_values, eigen_vectors, variables), 0.8)
  expect_error(expect_equal(proportional_explained_variance(eigen_values, eigen_vectors, variables, type="test")))
})

test_that("explained_variance.exact", {
  eigen_vectors <- matrix(c(1:9), nrow=3, ncol=3)
  eigen_values <- c(0.5, 0.3, 0.2)
  variables <- c(1, 2)

  expect_equal(explained_variance.exact(eigen_values, eigen_vectors, variables), 73.8)
  expect_equal(explained_variance.exact(eigen_values, c(), variables), 0)
})

test_that("vector_not_empty", {
  expect_equal(vector_not_empty(c(1, 2, 3)), TRUE)
  expect_equal(vector_not_empty(3), TRUE)
  expect_equal(vector_not_empty(c()), FALSE)
})

test_that("weighted_explained_variance", {
  eigen_vectors <- matrix(c(1:9), nrow=3, ncol=3)
  eigen_values <- c(0.5, 0.3, 0.2)
  variables <- c(1, 2)

  expect_equal(weighted_explained_variance(eigen_vectors, eigen_values, variables), 73.8)
})

test_that("explained_variance.approx", {
  eigen_values <- c(2, 3, 1, 0, 4)
  variables <- c(1,2)

  expect_equal(explained_variance.approx(eigen_values, variables), 5)
})