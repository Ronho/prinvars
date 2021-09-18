test_that("select_eigen_vector_scaling", {
  
})

test_that("scale_eigen_vectors", {
  eigen_vectors <- matrix(c(1:9), nrow=3, ncol=3)
  expected_result <- matrix(c(
      1/3, 2/3, 1,
      4/6, 5/6, 1,
      7/9, 8/9, 1
      ), nrow=3, ncol=3)

  expect_equal(
    scale_eigen_vectors(eigen_vectors=eigen_vectors),
    expected_result
  )
})