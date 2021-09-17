select_thresholding <- function(eigen_vectors, threshold, mode) {
  valid_threshold(threshold=threshold)
  eigen_vectors = switch(mode,
    "cutoff"=cutoff(x=eigen_vectors, threshold=threshold),
    "percentage"=percentage_per_eigen_vector(eigen_vectors=eigen_vectors, threshold=threshold),
    err_wrong_mode(mode=mode)
  )

  return(eigen_vectors)
}

valid_threshold <- function(threshold) {
  if (threshold > 1 || threshold < 0) {
    warning("Threshold should be between 0 and 1.")
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

err_wrong_mode <- function(mode) {
  stop(
    paste(
      "'", mode, "'", " is not a valid value for threshold_mode.",
      sep=""
    )
  )
}