test_that("get_indices", {
  x <- matrix(c(1:16), nrow=4, ncol=4)
  blocks <- list(
    new("Block", features=c(1, 4)),
    new("Block", features=c(2)),
    new("Block", features=c(3))
  )
  object <- list(
    x=x,
    blocks=blocks
  )
  class(object) <- "pla"

  expect_error(get_indices(object=object, block_indices=c(4)))
  expect_error(get_indices(object=object, block_indices=c()))
  expect_equal(get_indices(object=object, block_indices=c(1, 3)), c(1, 4, 3))
})

test_that("check_indices", {
  expect_error(check_indices(indices=c(), max_length=0))
  expect_error(check_indices(indices=NULL, max_length=1))
  expect_error(check_indices(indices=c(2), max_length=1))
  expect_silent(check_indices(indices=c(1), max_length=1))
})

test_that("err_must_provide_indices", {
  expect_error(err_must_provide_indices())
})

test_that("err_index_out_of_bounds", {
  expect_error(err_index_out_of_bounds())
})

test_that("conditional_matrix", {
  data <- matrix(c(1:9), nrow=3, ncol=3)
  expected_result <- matrix(c(13, 18, 22, 30), nrow=2, ncol=2)
  
  expect_equal(conditional_matrix(x=data, indices=c(1)), expected_result)
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
    selected_features=selected_features,
    is_valid=FALSE,
    ev_influenced=c(1, 2)
  )
  blocks[[2]] <- create_block(
    feature_names=c(1:5),
    selected_features=selected_features,
    is_valid=TRUE,
    ev_influenced=c(4, 5)
  )

  expect_s4_class(blocks[[1]], "Block")
  expect_s4_class(blocks[[2]], "Block")
  expect_equal(blocks[[1]]@features, c("b", "c"))
  expect_equal(blocks[[1]]@explained_variance, 0)
  expect_equal(blocks[[1]]@is_valid, FALSE)
  expect_equal(blocks[[1]]@ev_influenced, c(1, 2))
  expect_equal(blocks[[2]]@features, c(2, 3))
  expect_equal(blocks[[2]]@explained_variance, 0)
  expect_equal(blocks[[2]]@is_valid, TRUE)
  expect_equal(blocks[[2]]@ev_influenced, c(4, 5))
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

test_that("str_loadings", {
  expect_type(
    str_loadings(
      loadings=matrix(c(1:16), nrow=4, ncol=4),
      threshold=0.5,
      threshold_mode="cutoff",
      feature_names=c("a", "b", "c", "d"),
      C=c(1:4)
    ),
    "character"
  )
})

test_that("select_sparse_type_orthogonal", {
  expect_true(
    select_sparse_type_orthogonal(type="data")
  )
  expect_false(
    select_sparse_type_orthogonal(type="dispersion")
  )
  expect_error(
    select_sparse_type_orthogonal(type="anything")
  )
})

test_that("select_sparse_type_not_orthogonal", {
  expect_equal(
    select_sparse_type_not_orthogonal(type="data"),
    "predictor"
  )
  expect_equal(
    select_sparse_type_not_orthogonal(type="dispersion"),
    "Gram"
  )
  expect_error(
    select_sparse_type_not_orthogonal(type="anything")
  )
})

test_that("err_wrong_sparse_type", {
  expect_error(err_index_out_of_bounds())
})