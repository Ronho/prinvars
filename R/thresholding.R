select_thresholding <- function(eigen_vectors, threshold, mode) {
  valid_threshold(threshold=threshold)
  threshold_matrix = switch(
    tolower(mode),
    "cutoff"=cutoff(x=eigen_vectors, threshold=threshold),
    "percentage"=percentage_per_eigen_vector(
      eigen_vectors=eigen_vectors,
      threshold=threshold
    ),
    err_wrong_mode(mode=mode)
  )

  return(threshold_matrix)
}

valid_threshold_matrix <- function(threshold_matrix) {
  row_zeros <- which(get_zero_count(threshold_matrix) == ncol(threshold_matrix))
  col_zeros <- which(
    get_zero_count(t(threshold_matrix)) == nrow(threshold_matrix)
  )

  if(length(row_zeros) > 0 || length(col_zeros) > 0) {
    err_invalid_pla()
  }
}

valid_threshold <- function(threshold) {
  if (threshold > 1 || threshold < 0) {
    stop("Threshold must be between 0 and 1.")
  }
}

cutoff <- function(x, threshold) {
  x <- ifelse(abs(x) > threshold, 1, 0)

  return(x)
}

percentage_per_eigen_vector <- function(eigen_vectors, threshold) {
  eigen_vectors = apply(
    eigen_vectors,
    MARGIN=2,
    FUN=function(eigen_vector) {
      boundary <- max(abs(eigen_vector)) * threshold
      eigen_vector <- cutoff(x=eigen_vector, threshold=boundary)
      return(eigen_vector)
    }
  )

  return(eigen_vectors)
}

err_invalid_pla <- function() {
  stop("PLA is not valid. Threshold might be too large.")
}

err_wrong_mode <- function(mode) {
  stop(
    paste(
      "'",
      mode,
      "'",
      " is not a valid value for threshold_mode.",
      sep=""
    )
  )
}