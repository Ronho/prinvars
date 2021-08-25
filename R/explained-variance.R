proportional_explained_variance <- function(eigen_values, eigen_vectors, variables, type = "approx") {
    explained_variance <- 0
    switch(type,
            "approx" = { explained_variance <- explained_variance.approx(eigen_values, variables) },
            "exact" = { explained_variance <- explained_variance.exact(eigen_values, eigen_vectors, variables) }, {
                error <- paste("'", type, "'", " is not a valid value for explained variance.", sep = "")
                stop(error)
            })
    proportional_explained_variance = explained_variance / sum(eigen_values)
    return(proportional_explained_variance)
}

explained_variance.exact <- function(eigen_values, eigen_vectors, variables) {
  explained_variance <- 0
  if (vector_not_empty(eigen_vectors)) {
      explained_variance <- weighted_explained_variance(eigen_vectors, eigen_values, variables)
  }
  return(explained_variance)
}

vector_not_empty <- function(x) {
    return(length(x) > 0)
}

weighted_explained_variance <- function(eigen_vectors, eigen_values, variables) {
    explained_variance <- 0
    for (col in 1:ncol(eigen_vectors)) {
      summed_eigen_vector_elements <- sum(eigen_vectors[variables, col])
      explained_variance <- explained_variance + (eigen_values[col] * summed_eigen_vector_elements ^ 2)
    }
    return(explained_variance)
}

explained_variance.approx <- function(eigen_values, variables) {
  explained_variance <- sum(eigen_values[variables])
  return(explained_variance)
}