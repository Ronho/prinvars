calculate_explained_variance <- function(blocks, eigen, feature_names, type) {
  blocks <- lapply(blocks, function(block) {
    feature_idxs <- match(block@features, feature_names)
    block@explained_variance <- proportional_explained_variance(
      eigen=eigen,
      feature_idxs=feature_idxs,
      type=type
    )

    return(block)
  })

  return(blocks)
}

proportional_explained_variance <- function(eigen, feature_idxs, type) {
    explained_variance <- 0
    switch(
      tolower(type),
      "approx" = {
        explained_variance <- explained_variance.approx(
          eigen_values=eigen$values,
          feature_idxs=feature_idxs
      )},
      "exact" = {
        explained_variance <- explained_variance.exact(
          eigen=eigen,
          feature_idxs=feature_idxs
      )},
      {
        err_wrong_type(type)
      }
    )
    proportional_explained_variance <- explained_variance/sum(eigen$values)

    return(proportional_explained_variance)
}

err_wrong_type <- function(type) {
  stop(
    paste(
      "'", type, "'", " is not a valid value for explained variance.",
      sep=""
    )
  )
}

explained_variance.exact <- function(eigen, feature_idxs) {
  explained_variance <- 0
  if (vector_not_empty(x=eigen$vectors)) {
      explained_variance <- weighted_explained_variance(
        eigen=eigen,
        feature_idxs=feature_idxs
      )
  }

  return(explained_variance)
}

vector_not_empty <- function(x) {
    return(length(x) > 0)
}

weighted_explained_variance <- function(eigen, feature_idxs) {
    explained_variance <- 0
    for (col in 1:ncol(eigen$vectors)) {
      summed_eigen_vector_elements <- sum(eigen$vectors[feature_idxs, col]^2)
      explained_variance <- explained_variance + (
        eigen$values[col] * (summed_eigen_vector_elements) )
    }

    return(explained_variance)
}

explained_variance.approx <- function(eigen_values, feature_idxs) {
  explained_variance <- sum(eigen_values[feature_idxs])
  
  return(explained_variance)
}