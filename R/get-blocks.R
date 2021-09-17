# get block structure of the matrix
# x - transformed data, matrix
get_blocks <- function(threshold_matrix, feature_names, check) {
  ones <- 1
  number_features <- nrow(threshold_matrix)
  zero_counts <- get_zero_count(eigen_vectors=threshold_matrix)
  untaken_features <- 1:number_features
  blocks <- list()

  # iterate for each nxn combination
  while (ones <= number_features) {
    if (are_enough_features(untaken_features=untaken_features, ones=ones)) {
      zeros <- number_features-ones
      eligible_features <- get_eligible_features(
        zero_counts=zero_counts,
        zeros=zeros,
        untaken_features=untaken_features
      )

      combinations <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        zeros=zeros,
        current_combination=c(),
        check=check
      )

      for (combination in combinations) {
        untaken_features <- untaken_features[!untaken_features %in% combination]
        blocks[[length(blocks)+1]] <- create_block(
          feature_names=feature_names,
          selected_features=combination
        )
      }

      ones <- ones+1
    } else {
      blocks[[length(blocks)+1]] <- create_block(
        feature_names=feature_names,
        selected_features=untaken_features
      )
      ones <- number_features+1
    }
  }

  return(blocks)
}

# find combination that matches the ones
# x - data, matrix
# eligible_features - columns that match the required structure
# ones - number of 1s
# zeros - number of 0s
# current_combination - sequence of column indexes that match the required structure
find_combination <- function(threshold_matrix, eligible_features, ones, zeros, current_combination, check) {
  remaining_length <- ones - length(current_combination)

  if (remaining_length < 1) {
    valid_combination <- is_valid_combination(
      threshold_matrix=threshold_matrix,
      current_combination=current_combination,
      zeros=zeros,
      check=check
    )
    result <- if (valid_combination) current_combination else FALSE
    
    return(list(result))
  } else {
    combinations <- list()
    enough_features <- are_enough_eligible_features(
      eligible_features=eligible_features,
      remaining_length=remaining_length
    )
    while (enough_features) {
      current_combination <- c(current_combination, eligible_features[[1]])
      eligible_features <- eligible_features[-1]

      results <- find_combination(
        threshold_matrix=threshold_matrix,
        eligible_features=eligible_features,
        ones=ones,
        zeros=zeros,
        current_combination=current_combination,
        check=check
      )

      for (result in results) {
        if (result[1] != FALSE) {
          eligible_features <- eligible_features[!eligible_features %in% result]
          combinations[[length(combinations)+1]] <- result
          remaining_length <- ones
          current_combination <- c()
        } else {
          current_combination <- current_combination[
            -length(current_combination)
          ]
        }
      }
    }

    return(combinations)
  }
}

# check if there are enough features after finding a sequence that fits
are_enough_features <- function(untaken_features, ones) {
    return (length(untaken_features) >= 2*ones)
}

are_enough_eligible_features <- function(eligible_features, required_length) {
    return(
      !is.null(eligible_features) & 
      length(eligible_features) >= required_length
    )
}

get_eligible_features <- function(zero_counts, zeros, untaken_features) {
    eligible <- which(zero_counts >= zeros)
    eligible <- intersect(eligible, untaken_features)
    return(eligible)
}

is_valid_combination <- function(
  threshold_matrix,
  current_combination,
  zeros,
  check) {
  row_combination <- sum_vectors(
    x=threshold_matrix,
    indices=current_combination
  )

  if (are_exact_zeros(vector=row_combination, zeros=zeros)) {
    if (check_cols(check=check)) check_column_combination(
      threshold_matrix=threshold_matrix,
      row_combination=row_combination,
      zeros=zeros,
      current_combination=current_combination
    )

    return(TRUE)
  } else {

    return(FALSE)
  }
}

check_cols <- function(check) {
  result <- switch(
    check,
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

# tested
check_column_combination <- function(
  threshold_matrix,
  row_combination,
  zeros,
  current_combination) {
  column_combination <- sum_vectors(
    x=t(threshold_matrix),
    indices=which(row_combination >= 1)
  )
  
  if (!are_exact_zeros(vector=column_combination, zeros=zeros)) {
    warn_wrong_column_wise_combination(current_combination=current_combination)
  }
}

# tested
warn_wrong_column_wise_combination <- function(current_combination) {
  combination <- paste(unlist(current_combination), collapse=", ")
  warning(
    paste(
      "(",
      combination,
      ")",
      " is not a valid combination considering columns. 
      However, it is valid row-wise.",
      sep=""
    )
  )
}

# tested
are_exact_zeros <- function(vector, zeros) {
    return(length(which(vector == 0)) == zeros)
}

# sum vectors of matrix row-wise
# tested
sum_vectors <- function(x, indices) {
  if (length(indices) > 1) {
    result <- colSums(x[indices, ])
  } else (
    result <- x[indices, ]
  )

  return(result)
}