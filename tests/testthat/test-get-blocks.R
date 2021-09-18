test_that("get_blocks", {
  
})

test_that("find_combination", {
  
})

test_that("are_enough_features", {
  
})

test_that("are_enough_eligible_features", {
  
})

test_that("get_eligible_features", {
  
})

test_that("check_cols", {
  
})

test_that("err_wrong_check", {
  
})


test_that("sum_vectors", {
  matrix <- matrix(c(
    1, 0, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow=3, ncol=3)

  expect_equal(sum_vectors(x=matrix, indices=c(1,3)), c(1, 1, 2))
  expect_equal(sum_vectors(x=matrix, indices=c(3)), c(0, 1, 1))
  expect_equal(sum_vectors(x=c(), indices=c(4)), NULL)
  expect_error(sum_vectors(x=matrix, indices=c(4)))
  expect_error(sum_vectors(x=matrix(c()), indices=c(4)))
})

test_that("are_exact_zeros", {
  vector <- c(1, 1, 0, 0, 0)
  zeros <- 3

  expect_true(are_exact_zeros(vector=vector, zeros=zeros))
  expect_false(are_exact_zeros(vector=vector, zeros=zeros+1))
})

test_that("warn_wrong_column_wise_combination", {
  combination <- c(1, 2, 3)
  expect_warning(
    warn_wrong_column_wise_combination(current_combination=combination)
  )
})

test_that("check_column_combination", {
  threshold_matrix <- matrix(c(
    1, 0, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow=3, ncol=3)
  zeros <- 1

  expect_silent(
    check_column_combination(
      threshold_matrix=threshold_matrix,
      row_combination=c(1, 1, 0),
      zeros=zeros,
      current_combination=c(2, 3)
    )
  )
  expect_warning(
    check_column_combination(
      threshold_matrix=threshold_matrix,
      row_combination=c(0, 1, 0),
      zeros=zeros,
      current_combination=c(2)
    )
  )
})

test_that("is_valid_combination", {
  threshold_matrix <- matrix(c(
    1, 1, 0, 0,
    0, 0, 1, 0,
    0, 1, 0, 0,
    0, 0, 1, 0
  ), nrow=4, ncol=4)
  check <- "rnc"
  zeros <- 2

  expect_true(
    is_valid_combination(
        threshold_matrix=threshold_matrix,
        current_combination=c(1, 2),
        zeros=zeros,
        check=check
    )
  )
  expect_warning(
    is_valid_combination(
      threshold_matrix=threshold_matrix, 
      current_combination=c(3, 4),
      zeros=zeros,
      check=check
    )
  )
  expect_false(
    is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=c(1),
      zeros=zeros,
      check=check
    )
  )
  expect_silent(
    is_valid_combination(
        threshold_matrix=threshold_matrix,
        current_combination=c(1, 2),
        zeros=zeros,
        check=check
    )
  )
})