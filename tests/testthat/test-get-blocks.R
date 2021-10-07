
test_that("get_blocks", {
  matrix <- matrix(c(
    1, 1, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow=3, ncol=3)

  expected_result <- list()
  expected_result[[1]] <- new("Block", features=c(1), is_valid=FALSE)
  expected_result[[2]] <- new("Block", features=c(3), is_valid=TRUE)
  expected_result[[3]] <- new("Block", features=c(2), is_valid=FALSE)
  
  expect_equal(
    get_blocks(threshold_matrix=matrix, feature_names=c(1:nrow(matrix)), check="rows"),
    expected_result
  )
})

test_that("find_combination", {
  matrix <- matrix(c(
    1, 1, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow=3, ncol=3)

  expect_equal(
    find_combination(
      threshold_matrix=matrix,
      eligible_features=c(1:3),
      ones=1,
      current_combination=c(),
      check="rows"
    ),
    list(combination=c(1), is_valid=FALSE)
  )
})

test_that("get_eligible_features", {
  expect_equal(
    get_eligible_features(
      zero_counts=c(1, 2, 3, 2, 4),
      zeros=2,
      untaken_features=c(1:3)
    ),
    c(2, 3)
  )
  expect_equal(
    get_eligible_features(
      zero_counts=c(1, 2, 3, 2, 4),
      zeros=2,
      untaken_features=c()
    ),
    vector("integer")
  )
})

test_that("are_enough_features", {
  expect_true(are_enough_features(untaken_features=c(1:5), ones=2))
  expect_false(are_enough_features(untaken_features=c(1:5), ones=3))
  expect_true(are_enough_features(untaken_features=c(), ones=0))
})

test_that("are_enough_eligible_features", {
  expect_true(
    are_enough_eligible_features(eligible_features=c(1:5), required_length=3)
  )
  expect_false(
    are_enough_eligible_features(eligible_features=c(), required_length=1)
  )
  expect_false(
    are_enough_eligible_features(eligible_features=c(), required_length=0)
  )
})

test_that("check_cols", {
  expect_equal(check_cols(check="rnc"), TRUE)
  expect_equal(check_cols(check="rows"), FALSE)
  expect_error(check_cols(check="Test"))
})

test_that("err_wrong_check", {
  expect_error(err_wrong_check(check="Test"))
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

  expect_equal(
    is_valid_combination(
        threshold_matrix=threshold_matrix,
        current_combination=c(1, 2),
        check=check,
        ones=2
    ),
    c(TRUE, TRUE)
  )
  expect_equal(
    is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=c(1),
      check=check,
      ones=1
    ),
    c(TRUE, FALSE)
  )
})