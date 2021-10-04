matrix <- matrix(data=seq(0.1, 0.9, by=0.1), nrow=3, ncol=3)
cutoff_expected_matrix <- matrix(
    data=c(0, 0, 0, 0, 0, 1, 1, 1, 1), nrow=3, ncol=3
)
percentage_expected_matrix <- matrix(data=c(0, 1, 1, 1, 1, 1, 1, 1, 1), nrow=3, ncol=3)
threshold <- 0.5

test_that("select_thresholding", {
  expect_equal(
    select_thresholding(
      eigen_vectors=matrix,
      threshold=threshold,
      mode="percentage"
    ),
    percentage_expected_matrix
  )
  expect_equal(
    select_thresholding(
      eigen_vectors=matrix,
      threshold=threshold,
      mode="cutoff"
    ),
    cutoff_expected_matrix
  )
  expect_error(
    select_thresholding(eigen_vectors=matrix, threshold=threshold, mode="Test")
  )
})

test_that("valid_threshold", {
  expect_silent(valid_threshold(threshold=threshold))
  expect_error(valid_threshold(threshold=1.1))
  expect_error(valid_threshold(threshold=-0.1))
})

test_that("cutoff", {
  expect_equal(cutoff(x=matrix, threshold=threshold), cutoff_expected_matrix)
})

test_that("percentage_per_eigen_vector", {
  expect_equal(
    percentage_per_eigen_vector(eigen_vectors=matrix, threshold=threshold),
    percentage_expected_matrix
  )
})

test_that("err_wrong_mode", {
  expect_error(err_wrong_mode(mode="Test"))
})