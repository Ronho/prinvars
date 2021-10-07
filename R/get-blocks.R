get_blocks <- function(threshold_matrix, feature_names, check) {
  num_features <- nrow(threshold_matrix)
  untaken_features <- 1:num_features
  zero_counts <- get_zero_count(threshold_matrix)
  ones <- 1
  blocks <- list()

  while (length(untaken_features) > 0 && ones <= length(untaken_features)) {
    eligible_features <- get_eligible_features(
      zero_counts=zero_counts,
      zeros=num_features-ones,
      untaken_features=untaken_features
    )

    while (length(eligible_features) >= ones) {
      combination <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        current_combination=c(),
        check=check
      )

      if (!is.atomic(combination)) {
        is_valid <- combination$is_valid
        combination <- combination$combination
        blocks[[length(blocks)+1]] <- create_block(feature_names=feature_names, selected_features=combination, is_valid=is_valid)
        eligible_features <- eligible_features[!eligible_features %in% combination]
        untaken_features <- untaken_features[!untaken_features %in% combination]
      } else {
        eligible_features <- eligible_features[-1]
      }
    }

    ones <- ones+1
  }

  if(length(untaken_features) > 0) {
    blocks[[length(blocks)+1]] <- create_block(feature_names=feature_names, selected_features=untaken_features, is_valid=FALSE)
  }

  return(blocks)
}

find_combination <- function(
  threshold_matrix,
  eligible_features,
  ones,
  current_combination,
  check) {
  num_remaining_features <- ones - length(current_combination)

  if (num_remaining_features == 0) {
    is_valid <- is_valid_combination(threshold_matrix=threshold_matrix, current_combination=current_combination, ones=ones, check=check)

    if (check_cols(check=check)) {
      is_valid <- is_valid[1] & is_valid[2]
      if (is_valid) {
        result <- list(combination=current_combination, is_valid=is_valid)
        return(result)
      }
    } else {
      if (is_valid[1]) {
        result <- list(combination=current_combination, is_valid=is_valid[2])
        return(result)
      }
    }

    return(FALSE)
  } else {
    while(length(eligible_features) > 0 && length(eligible_features) >= num_remaining_features) {
      current_combination <- c(current_combination, eligible_features[1])
      eligible_features <- eligible_features[-1]

      result <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        current_combination=current_combination,
        check=check
      )

      if (!is.atomic(result)) {
        return(result)
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
  check) {
  # check if rows are valid
  rows <- threshold_matrix[current_combination, ]
  if (length(current_combination) > 1) {
    rows <- colSums(rows)
  }
  ev_influenced <- which(rows >= 1)
  row_valid <- length(ev_influenced) == ones

  # check if columns are valid
  cols <- t(threshold_matrix)[ev_influenced, ]
  if (length(ev_influenced) > 1) {
    cols <- colSums(cols)
  }
  col_valid <- length(which(cols >= 1)) == ones

  return(c(row_valid, col_valid))
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