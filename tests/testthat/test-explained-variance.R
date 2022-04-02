eigen <- list(
  vectors=matrix(c(1:9), nrow=3, ncol=3),
  values=c(0.5, 0.3, 0.2)
)
feature_idxs <- c(1, 2)
threshold_matrix <- matrix(c(1,0,0,0,1,0,0,0,1), nrow=3, ncol=3)

test_that("calculate_explained_variance", {
  block <- new("Block", features=feature_idxs)
  blocks <- calculate_explained_variance(
    blocks=list(block),
    eigen=eigen,
    feature_names=c(1:3),
    type="approx",
    threshold_matrix=threshold_matrix
  )

  blocks_2 <- calculate_explained_variance(
    blocks=list(block),
    eigen=eigen,
    feature_names=c(1:3),
    type="approx",
    threshold_matrix=threshold_matrix,
    is_absolute=TRUE
  )

  expect_equal(blocks[[1]]@explained_variance, 0.8)
  expect_equal(blocks_2[[1]]@explained_variance, 0.8)
})

test_that("err_wrong_type", {
  expect_error(err_wrong_type(""))
})

test_that("proportional_explained_variance", {
  expect_equal(
      proportional_explained_variance(
        eigen=eigen,
        feature_idxs=feature_idxs,
        type="exact"
      ),
      37.4
  )
  expect_equal(
      proportional_explained_variance(
        eigen=eigen,
        feature_idxs=feature_idxs,
        type="approx",
        threshold_matrix=threshold_matrix
      ),
      0.8
  )
  expect_error(
    expect_equal(
      proportional_explained_variance(
        eigen=eigen,
        feature_idxs=feature_idxs,
        type="test"
      )
    )
  )
})

test_that("explained_variance.exact", {
  expect_equal(
    explained_variance.exact(eigen=eigen, feature_idxs=feature_idxs), 
    37.4
  )
})

test_that("vector_not_empty", {
  expect_equal(vector_not_empty(x=c(1, 2, 3)), TRUE)
  expect_equal(vector_not_empty(x=3), TRUE)
  expect_equal(vector_not_empty(x=c()), FALSE)
})

test_that("weighted_explained_variance", {
  expect_equal(
      weighted_explained_variance(eigen=eigen, feature_idxs=feature_idxs), 
      37.4
  )
})

test_that("explained_variance.approx", {
  eigen_values <- c(2, 3, 1, 0, 4)
  feature_idxs <- c(1,2)

  expect_equal(
    explained_variance.approx(
        eigen_values=eigen_values,
        feature_idxs=feature_idxs
    ),
    5
  )
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