get_feature_names <- function(x) {
  feature_names <- colnames(x)
  feature_names <- if (length(feature_names) > 0) colnames(x) else 1:ncol(x)
  return(feature_names)
}

manipulate_matrix <- function(x, manipulator = "cov") {
  if (manipulator == "cor") return(cor(x, method = c("pearson")))
  if (manipulator == "cov") return(cov(x, method = c("pearson")))
  if (manipulator == "none") return(x)
  stop(paste("'", manipulator, "'", " is not a valid value for manipulator.", sep = ""))
}

create_block <- function(feature_names, selected_features) {
  if (length(feature_names) > 0) {
    selected_features <- feature_names[selected_features]
  }
  return <- new("Block", features=selected_features)
}

# scale eigen vectors (column wise) between -1 and 1
scale_eigen_vectors <- function(eigen_vectors) {
  scaled_eigen_vectors <- apply(
    eigen_vectors,
    MARGIN = 2,
    FUN = function(x) {
      x / max(abs(x))
    }
  )
  return(scaled_eigen_vectors)
}

get_threshold_matrix <- function(eigen_vectors, threshold) {
  if (threshold > 1 || threshold < 0)
    warning("Threshold should be between 0 and 1.")
  x <- ifelse(abs(eigen_vectors) > threshold, 1, 0)
  return(x)
}

# get number of zeros for each eigenvector
get_zero_count <- function(eigen_vectors) {
  zero_length <- apply(
    eigen_vectors,
    MARGIN = 2,
    FUN = function(x) {
      length(which(x == 0))
    }
  )
  return <- zero_length
}

# sum vectors of matrix rowwise
sum_vectors <- function(x, indices) {
  if (length(indices) > 1) {
    return <- colSums(x[indices, ])
  } else (
    return <- x[indices, ]
  )
}

check_pla_equality <- function(a, b) {
  is_equal <- TRUE
  for (a_block in a$blocks) {
    found <- FALSE
    for (b_block in b$blocks) {
      if ((a_block@explained_variance == b_block@explained_variance) && all(a_block@indices == b_block@indices)) {
        found <- TRUE
      }
    }
    is_equal <- if (found == FALSE) found else is_equal
  }
  return(is_equal)
}

apply_type <- function(eigen, scaled_ev, threshold, colnames, expvar, type) {
  if (type == "columns") {
    result <- transform(eigen, scaled_ev, threshold, colnames, expvar)
  } else if (type == "rows") {
    eigen$vectors <- t(eigen$vectors)
    result <- transform(eigen, scaled_ev, threshold, colnames, expvar)
  } else if (type == 'rnc') {
    pla_cols <- transform(eigen, scaled_ev, threshold, colnames, expvar)

    eigen$vectors <- t(eigen$vectors)
    pla_rows <- transform(eigen, scaled_ev, threshold, colnames, expvar)
    
    if (check_pla_equality(pla_cols, pla_rows)) {
      result <- pla_cols
    } else {
      stop("PLA for column and rows is not equal.")
    }
  } else {
    stop(paste("'", type, "'", " is not a valid value for type. It can be either 'columns' (DEFAULT), 'rows' or 'rnc'", sep = ""))
  }
  return(result)
}

transform <- function(eigen, scaled_ev, threshold, colnames, expvar) {
  eigen$vectors <- if (scaled_ev) scale_eigen_vectors(eigen$vectors) else eigen$vectors
  x <- get_threshold_matrix(eigen$vectors, threshold)
  blocks <- get_blocks(x, colnames)
  blocks <- calculate_explained_variance(blocks, eigen, colnames, expvar)

  result <- list(
    "eigen_vectors" = eigen$vectors,
    "threshold" = threshold,
    "blocks" = blocks
  )
  class(result) <- "pla"
  return(result)
}

calculate_explained_variance <- function(blocks, eigen, colnames, expvar) {
  blocks <- lapply(blocks, function(block) {
    col_idxs <- match(block@indices, colnames)
    block@explained_variance <- proportional_explained_variance(eigen$values,
                                                                eigen$vectors,
                                                                col_idxs,
                                                                expvar)
    return(block)
  })
  return(blocks)
}

# get block structure of the matrix
# x - transformed data, matrix
get_blocks <- function(x, colnames) {
  n <- 1 # number of 1s
  m <- nrow(x) # dimension of columns
  zero_length <- get_zeros(x) # array of integers indicating the amount of zeros in each column
  untaken_cols <- 1:ncol(x) # columns that are not already part of a block
  blocks <- list()

  # iterate for each nxn combination
  while (n < m) {
    if (length(untaken_cols) >= 2*n) { #check if there are enough columns after finding a sequence that fits
      dimension <- m-n # number of 0s
      eligible_cols <- which(zero_length >= dimension) # columns that have the required dimension
      eligible_cols <- intersect(eligible_cols, untaken_cols)

      # find a combination that works by using different starting points
      while (!is.null(eligible_cols) & length(eligible_cols) >= n) {
        element <- eligible_cols[[1]]
        eligible_cols <- eligible_cols[-1]

        value <- find_combination(x, possible_data = eligible_cols, required_length = n, dimension = dimension, current_combination = c(element))

        if (value[1] != FALSE) {
          blocks[[length(blocks)+1]] <- create_block(cols = value, colnames = colnames)
          untaken_cols <- untaken_cols[!untaken_cols %in% value]
          eligible_cols <- eligible_cols[!eligible_cols %in% value]
        }
      }

      n <- n+1
    } else {
      blocks[[length(blocks)+1]] <- create_block(cols = untaken_cols, colnames = colnames)
      n <- m
    }
  }

  return <- blocks
}

# find combination that matches the required_length
# x - data, matrix
# possible_data - columns that match the required structure
# required_length - number of 1s
# dimension - number of 0s
# current_combination - sequence of column indexes that match the required structure
find_combination <- function(x, possible_data, required_length, dimension, current_combination) {
  remaining_length <- required_length - length(current_combination)
  if (remaining_length < 1) {
    v <- sum_vectors(x, current_combination)

    if (length(which(v == 0)) == dimension) {
      return(current_combination)
    } else {
      return(FALSE)
    }
  } else {
    while (!is.null(possible_data) & length(possible_data) >= remaining_length) {
      element <- possible_data[[1]]
      possible_data <- possible_data[-1]

      result <- find_combination(x, possible_data = possible_data, required_length = required_length, 
      dimension = dimension, current_combination = c(current_combination, element))

      if (result[1] != FALSE) {
        return(result)
      }
    }
    return(FALSE)
  }
}