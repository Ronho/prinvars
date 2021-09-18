test_that("select_manipulator", {
  matrix <- matrix(c(1:9), nrow=3, ncol=3)

  expect_equal(select_manipulator(x=matrix, manipulator="cov"), cov(matrix))
  expect_equal(select_manipulator(x=matrix, manipulator="cor"), cor(matrix))
  expect_equal(select_manipulator(x=matrix, manipulator="none"), matrix)
  expect_error(select_manipulator(x=matrix, manipulator=""))  
})
test_that("validate_matrix_structure", {
  
})
test_that("is_quadratic_matrix", {
  
})
test_that("wrong_manipulator", {
  
})
