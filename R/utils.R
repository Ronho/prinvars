transform <- function(eigen, scaled_ev, threshold, colnames, expvar) {
  # scale eigen vectors (column wise) between -1 and 1
  if (scaled_ev)
    eigen$vectors = apply(
      eigen$vectors,
      MARGIN = 2,
      FUN = function(x) {
        x / max(abs(x))
      }
    )

  # create threshold matrix: each element is assigned a 0 or 1 depending
  # if it is underneath or above the threshold
  if (threshold > 1 || threshold < 0)
    warning("Threshold should be between 0 and 1.")
  x <- ifelse(abs(eigen$vectors) > threshold, 1, 0)

  blocks <- get_blocks(x, colnames=colnames)


  # get explained variance for each block
  blocks <- lapply(blocks, function(block) {
    col_idxs <- match(block@columns, colnames)
    block@explained_variance <- explained_variance(eigen$values,
                                                    eigen$vectors,
                                                    col_idxs,
                                                    expvar)
    return(block)
  })

  result <- list(
    "eigen_vectors" = eigen$vectors,
    "threshold" = threshold,
    "blocks" = blocks
  )
  class(result) <- "pla"
  return(result)
}

manipulate_matrix <- function(x, manipulator = "cov") {
  if (manipulator == "cor") return(cor(x, method = c("pearson")))
  if (manipulator == 'cov') return(cov(x, method = c("pearson")))
  stop(paste("'", manipulator, "'", " is not a valid value for manipulator.", sep = ""))
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
      dimension = m-n # number of 0s
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

# get number of zeros for each eigenvector
# x - data, matrix
get_zeros <- function(x) {
  zero_length <- c()
  for (idx in 1:ncol(x)) {
    zero_positions <- which(x[, idx] == 0)
    zero_length <- c(zero_length, length(zero_positions))
  }

  return <- zero_length
}

# find combination that matches the required_length
# x - data, matrix
# possible_data - columns that match the required structure
# required_length - number of 1s
# dimension - number of 0s
# current_combination - sequence of column indexes that match the required structure
find_combination <- function(x, possible_data, required_length, dimension, current_combination) {
  remaining_length = required_length - length(current_combination)
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

sum_vectors <- function(x, current_combination) {
  if (length(current_combination) > 1) {
    return <- rowSums(x[, current_combination])
  } else (
    return <- x[, current_combination]
  )
}

create_block <- function(cols, colnames) {
  if (length(colnames) > 0) {
    cols <- colnames[cols]
  }
  return <- new("Block", columns = cols)
}

explained_variance <- function(eigen_values, eigen_vectors, cols, type = "approx") {
  switch(type,
         "approx" = { return(explained_variance.approx(eigen_values, cols)) },
         "exact" = { return(explained_variance.exact(eigen_values, eigen_vectors, cols)) }, {
           stop(paste("'", type, "'", " is not a valid value for explained variance.", sep = ""))
         })
}

explained_variance.exact <- function(eigen_values, eigen_vectors, cols) {
  sum <- 0

  if (length(eigen_vectors) > 0) {
    for (col in 1:ncol(eigen_vectors)) {
      sum_vec <- 0
      for (row in cols) {
        sum_vec <- sum_vec + eigen_vectors[row, col]
      }
      sum <- sum + (eigen_values[col] * sum_vec ^ 2)
    }
  }

  return(sum / sum(eigen_values))
}

explained_variance.approx <- function(eigen_values, cols) {
  sum <- 0
  for (col in cols) {
    sum <- sum + eigen_values[col]
  }

  return(sum / sum(eigen_values))
}
