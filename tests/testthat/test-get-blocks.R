# Tests for get_blocks.
test_get_blocks_helper <- function(method) {
  with_mocked_bindings(
    get_blocks_legacy = function(threshold_matrix, feature_names, check) {
      "legacy result"
    },
    get_blocks_clustering = function(threshold_matrix, feature_names, check) {
      "clustering result"
    },
    {
      threshold_matrix <- matrix(1, nrow = 5, ncol = 5)
      feature_names <- c("a", "b", "c", "d", "e")
      check <- TRUE

      result <- get_blocks(threshold_matrix, feature_names, check, method)

      return(result)
    }
  )
}

test_that("get_blocks uses the legacy method correctly", {
  result <- test_get_blocks_helper("legacy")
  expect_equal(result, "legacy result")
})

test_that("get_blocks uses the clustering method correctly", {
  result <- test_get_blocks_helper("clustering")
  expect_equal(result, "clustering result")
})

test_that("get_blocks throws an error for an invalid method", {
  expect_error(test_get_blocks_helper("invalid"))
})

# Tests for get_features_in_cluster.
test_that("get_features_in_cluster correctly extracts single row", {
  cluster_merge <- matrix(c(-1, -2, 3, -4, 5, -6, 7, -8, 9, -10),
                          nrow = 5,
                          ncol = 2,
                          byrow = TRUE)
  result <- get_features_in_cluster(cluster_merge, 1)
  expect_equal(result$features, c(1, 2))
  expect_equal(result$associated_rows_in_cluster, c(1))
})

test_that("get_features_in_cluster correctly iterates over multiple rows", {
  cluster_merge <- matrix(c(-1, -2, -3, -4, 1, -5, 2, -6, -7, -8,
                            3, 4, 5, -9, -10, 6, 7, 8, 9, -11),
                          ncol = 2,
                          byrow = TRUE)
  result <- get_features_in_cluster(cluster_merge, 4)
  expect_equal(result$features, c(3, 4, 6))
  expect_equal(result$associated_rows_in_cluster, c(4, 2))
})

test_that("get_features_in_cluster correctly iterates over all rows", {
  cluster_merge <- matrix(c(-1, -2, -3, -4, 1, -5, 2, -6, -7, -8,
                            3, 4, 5, -9, -10, 6, 7, 8, 9, -11),
                          ncol = 2,
                          byrow = TRUE)
  result <- get_features_in_cluster(cluster_merge, 10)
  expect_equal(result$features, c(7, 8, 9, 10, 1, 2, 5, 3, 4, 6, 11))
  expect_equal(result$associated_rows_in_cluster, c(10, 9, 7, 5, 8,
                                                    6, 3, 1, 4, 2))
})

test_that("get_features_in_cluster handles empty cluster_merge", {
  result <- get_features_in_cluster(matrix(ncol=2), 0)
  expect_equal(result$features, integer(0))
  expect_equal(result$associated_rows_in_cluster, c(0))
})

test_that("get_features_in_cluster throws error when wrong index provided", {
  expect_error(get_features_in_cluster(matrix(ncol=2), 1))
})

test_that("get_blocks throws an error for an invalid method", {
  threshold_matrix <- matrix(1, nrow = 5, ncol = 5)
  feature_names <- c("a", "b", "c", "d", "e")
  check <- TRUE
  method <- "invalid"

  expect_error(get_blocks(threshold_matrix, feature_names, check, method))
})


test_that("get_blocks_legacy", {
  matrix <- matrix(c(
    1, 1, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow=3, ncol=3)

  expected_result <- list()
  expected_result[[1]] <- new("Block",
    features=c(1),
    is_valid=TRUE,
    ev_influenced=c(1))
  expected_result[[2]] <- new("Block",
    features=c(3),
    is_valid=TRUE,
    ev_influenced=c(3))
  expected_result[[3]] <- new("Block",
    features=c(2),
    is_valid=FALSE,
    ev_influenced=c(2))

  expect_equal(
    get_blocks_legacy(
      threshold_matrix=matrix,
      feature_names=seq_len(nrow(matrix)),
      check="rows"),
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
      check="rows",
      taken_features=c()
    ),
    list(combination=c(1), is_valid=TRUE, ev_influenced=c(1))
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

test_that("check_cols", {
  expect_equal(check_cols(check="rnc"), TRUE)
  expect_equal(check_cols(check="rows"), FALSE)
  expect_error(check_cols(check="Test"))
})

test_that("err_wrong_check", {
  expect_error(err_wrong_check(check="Test"))
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
    list(is_valid=c(TRUE, TRUE), ev_influenced=c(1, 3))
  )
  expect_equal(
    is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=c(1),
      check=check,
      ones=1
    ),
    list(is_valid=c(TRUE, FALSE), ev_influenced=c(1))
  )
})