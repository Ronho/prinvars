matrix <- matrix(c(1:9), nrow=3, ncol=3)

test_that("select_manipulator", {
  expect_equal(select_manipulator(x=matrix, manipulator="cov"), cov(matrix))
  expect_equal(select_manipulator(x=matrix, manipulator="cor"), cor(matrix))
  expect_equal(select_manipulator(x=matrix, manipulator="none"), matrix)
  expect_error(select_manipulator(x=matrix, manipulator=""))  
})

test_that("validate_matrix_structure", {
  expect_silent(validate_matrix_structure(x=matrix))
  expect_equal(validate_matrix_structure(x=matrix), matrix)
  expect_silent(validate_matrix_structure(x=matrix()))
  expect_warning(validate_matrix_structure(x=matrix(c(1:12), nrow=3, ncol=4)))
})

test_that("is_quadratic_matrix", {
  expect_true(is_quadratic_matrix(x=matrix))
  expect_true(is_quadratic_matrix(x=matrix()))
  expect_false(is_quadratic_matrix(x=matrix(c(1:12), nrow=3, ncol=4)))
})

test_that("wrong_manipulator", {
  expect_error(wrong_manipulator(manipulator="Test"))
})
