#' @title Decide which function to use to get blocks
#'
#' @description One can choose between "legacy" and "clustering". "legacy" is a
#' way of identifying blocks through a slow algorithm and "clustering" improves
#' this process by using clustering algorithm to find relevant combinations and
#' applying some optimizations. "legacy" can be removed at any time.
#'
#' @param threshold_matrix A m x m matrix containing 0 and 1 to represent
#' correlation between features.
#' @param feature_names A vector containing the names of the features in the
#' order that they are represented in the threshold_matrix.
#' @param check A string identifying how to prove a correct combination of
#' features.
#' @param method A string representing the method to use for retrieving block
#' structures out of the threshold_matrix. Can be "legacy" for the slower
#' variant or "clustering" for a faster implementation.
#'
#' @return
#' A list containing all found blocks.
#'
#' @keywords internal
get_blocks <- function(threshold_matrix, feature_names, check, method) {
  func <- switch(
    tolower(method),
    "legacy"=get_blocks_legacy,
    "clustering"=get_blocks_clustering,
    stop(
      sprintf(
        "Invalid argument: 'get_blocks_method' must be either 'legacy' or
        'clustering', not '%s'.",
        method
      )
    )
  )

  return(
    func(
      threshold_matrix=threshold_matrix,
      feature_names=feature_names,
      check=check
    )
  )
}

#' @title Get all features that correspond to a cluster.
#'
#' @description Finds all values feature indices that build a cluster together.
#' Per definition of the merge argument that is returned by hclust, value below
#' zero indicate the original variable and values above zero are a reference to
#' another cluster. Therefore, this function extracts the values building up a
#' cluster.
#'
#' @param cluster_merge A m x 2 matrix as retrieved by hclust that builds a tree
#' of clusters.
#' @param index The index of cluster_merge to find the features for.
#'
#' @return
#' A list containing
#' \item{features}{
#'   a list of integers identifying the features that correspond to the cluster.
#' }
#' \item{associated_rows_in_cluster}{
#'   rows of the cluster_merge that were used to get all features of the
#'   cluster.
#' }
#'
#' @keywords internal
get_features_in_cluster <- function(cluster_merge, index) {
  features <- integer(0)
  associated_rows <- c(index)
  for (value in cluster_merge[index, ]) {
    if (value < 0) {
      new_associations <- which(cluster_merge == value, arr.ind = TRUE)[, 1]
      features <- c(features, abs(value))
      associated_rows <- c(associated_rows, new_associations)
    } else {
      result <- get_features_in_cluster(cluster_merge=cluster_merge,
                                        index=abs(value))
      features <- c(features, result$features)
      associated_rows <- c(associated_rows, result$associated_rows_in_cluster)
    }
  }

  names(associated_rows) <- NULL

  return(
    list(
      features=unique(features),
      associated_rows_in_cluster=unique(associated_rows)
    )
  )
}

#' @title Get blocks using clustering.
#'
#' @description Finds all valid blocks using clustering algorithm. We apply some
#' optimization techniques such as searching for blocks of maximum size m/2.
#'
#' @param threshold_matrix A m x m matrix containing 0 and 1 to represent
#' correlation between features.
#' @param feature_names A vector containing the names of the features in the
#' order that they are represented in the threshold_matrix.
#' @param check A string identifying how to prove a correct combination of
#' features.
#'
#' @return
#' A list of block objects.
#'
#' @keywords internal
get_blocks_clustering <- function(threshold_matrix, feature_names, check) {
  num_features <- nrow(threshold_matrix)
  untaken_features <- 1:num_features
  blocks <- list()
  rows_removed <- integer(0)

  cluster <- hclust(dist(threshold_matrix))
  cluster_merge <- cluster$merge

  # Also include singleton clusters in the cluster merge structure.
  # Since singleton clusters are added, we need to shift the indices
  # respectively.
  singleton_clusters <- matrix(-1:-num_features, ncol=2, nrow=num_features)
  cluster_merge <- rbind(singleton_clusters, cluster_merge)
  cluster_merge[cluster_merge > 0] <- cluster_merge[cluster_merge > 0] +
    num_features

  for (index in seq_len(nrow(cluster_merge))) {
    # Skip row if already sorted out.
    if (index %in% rows_removed) next

    combination <- get_features_in_cluster(cluster_merge=cluster_merge,
                                           index=index)

    # Searching for clusters that have a size greater than half of the variables
    # is useless, because all remaining features must form a cluster.
    if (length(combination$features) > (num_features / 2)) {
      combination$features <- untaken_features
      all_indices <- seq_len(nrow(cluster_merge))
      combination$associated_rows_in_cluster <- all_indices[!all_indices %in%
                                                              rows_removed]
    }

    valid_combination <- is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=combination$features,
      ones=length(combination$features),
      check=check,
      taken_features=setdiff(c(1:num_features), untaken_features)
    )

    is_valid <- valid_combination$is_valid

    if (check_cols(check=check)) {
      is_valid <- is_valid[1] & is_valid[2]
    } else {
      is_valid <- is_valid[1]
    }

    if (is_valid) {
      blocks[[length(blocks) + 1]] <- create_block(
        feature_names=feature_names,
        selected_features=combination$features,
        is_valid=is_valid,
        ev_influenced=valid_combination$ev_influenced
      )
      untaken_features <- untaken_features[!untaken_features
                                           %in% combination$features]
      rows_removed <- c(rows_removed, combination$associated_rows_in_cluster)
      rows_removed <- unique(rows_removed)
    }
  }

  return(blocks)
}

get_blocks_legacy <- function(threshold_matrix, feature_names, check) {
  num_features <- nrow(threshold_matrix)
  untaken_features <- 1:num_features
  untaken_evs <- seq_len(ncol(threshold_matrix))
  zero_counts <- get_zero_count(threshold_matrix)
  ones <- 1
  blocks <- list()

  while (length(untaken_features) > 0 && ones <= length(untaken_features)) {
    eligible_features <- get_eligible_features(
      zero_counts=zero_counts,
      zeros=num_features - ones,
      untaken_features=untaken_features
    )

    while (length(eligible_features) >= ones) {
      combination <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        current_combination=c(),
        check=check,
        taken_features=setdiff(c(1:num_features), untaken_features)
      )

      if (!is.atomic(combination)) {
        is_valid <- combination$is_valid
        ev_influenced <- combination$ev_influenced
        combination <- combination$combination

        blocks[[length(blocks) + 1]] <- create_block(
          feature_names=feature_names,
          selected_features=combination,
          is_valid=is_valid,
          ev_influenced=ev_influenced
        )
        eligible_features <- eligible_features[
          !eligible_features %in% combination
        ]
        untaken_evs <- untaken_evs[!untaken_evs %in% ev_influenced]
        untaken_features <- untaken_features[!untaken_features %in% combination]
      } else {
        eligible_features <- eligible_features[-1]
      }
    }

    ones <- ones + 1
  }

  if (length(untaken_features) > 0) {
    blocks[[length(blocks) + 1]] <- create_block(
      feature_names=feature_names,
      selected_features=untaken_features,
      is_valid=FALSE,
      ev_influenced=untaken_evs
    )
  }

  return(blocks)
}

find_combination <- function(
  threshold_matrix,
  eligible_features,
  ones,
  current_combination,
  check,
  taken_features
) {
  num_remaining_features <- ones - length(current_combination)

  if (num_remaining_features == 0) {
    valid_combination <- is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=current_combination,
      ones=ones,
      check=check,
      taken_features=taken_features
    )
    is_valid <- valid_combination$is_valid

    if (check_cols(check=check)) {
      is_valid <- is_valid[1] & is_valid[2]
      if (is_valid) {
        result <- list(combination=current_combination, is_valid=is_valid,
          ev_influenced=valid_combination$ev_influenced)
        return(result)
      }
    } else {
      if (is_valid[1]) {
        result <- list(combination=current_combination, is_valid=is_valid[2],
          ev_influenced=valid_combination$ev_influenced)
        return(result)
      }
    }

    return(FALSE)
  } else {
    while (
      length(eligible_features) > 0 &&
      length(eligible_features) >= num_remaining_features
    ) {
      current_combination <- c(current_combination, eligible_features[1])
      eligible_features <- eligible_features[-1]

      result <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        current_combination=current_combination,
        check=check,
        taken_features=taken_features
      )

      if (!is.atomic(result)) {
        return(result)
      } else {
        current_combination <- current_combination[1:
          (length(current_combination) - 1)
        ]
      }
    }
    return(FALSE)
  }
}

get_eligible_features <- function(zero_counts, zeros, untaken_features) {
    eligible <- which(zero_counts >= zeros)
    eligible <- intersect(eligible, untaken_features)
    eligible <- as.vector(eligible, mode="integer")
    return(eligible)
}

is_valid_combination <- function(
  threshold_matrix,
  current_combination,
  ones,
  check,
  taken_features) {
  # check if rows are valid
  rows <- threshold_matrix[current_combination, ]
  if (length(current_combination) > 1) {
    rows <- colSums(rows)
  }
  ev_influenced <- which(rows >= 1)
  row_valid <- length(ev_influenced) == ones

  # check if columns are valid
  col_valid <- TRUE
  if (check_cols(check=check)) {
    cols <- t(threshold_matrix)[ev_influenced, ]
    if (length(ev_influenced) > 1) {
      cols <- colSums(cols)
    }
    col_valid <- length(which(cols >= 1)) == ones
  } else {
    if (length(taken_features) > 0) {
      cols <- t(threshold_matrix)[ev_influenced, taken_features]
      if (length(ev_influenced) > 1 && length(taken_features) > 1) {
        cols <- colSums(cols)
      }
      col_valid <- length(which(cols >= 1)) <= ones
    }
  }

  return(list(is_valid=c(row_valid, col_valid), ev_influenced=ev_influenced))
}

check_cols <- function(check) {
  result <- switch(
    tolower(check),
    "rnc"=TRUE,
    "rows"=FALSE,
    err_wrong_check(check=check)
  )

  return(result)
}

err_wrong_check <- function(check) {
  stop(
    paste(
      "'", check, "'", " is not a valid value for check.",
      sep=""
    )
  )
}
