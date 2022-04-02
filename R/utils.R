get_feature_names <- function(x) {
  feature_names <- colnames(x)
  feature_names <- if (length(feature_names) > 0) feature_names else 1:ncol(x)

  return(feature_names)
}

create_block <- function(feature_names, selected_features, is_valid) {
  if (length(feature_names) > 0) {
    selected_features <- feature_names[selected_features]
  }

  return(new("Block", features=selected_features, is_valid=is_valid))
}

# get number of zeros for each eigenvector
get_zero_count <- function(eigen_vectors) {
  zero_length <- apply(
    eigen_vectors,
    MARGIN=1,
    FUN=function(x) {
      length(which(x == 0))
    }
  )

  return(zero_length)
}

get_indices <- function(object, block_indices) {
  check_indices(indices=block_indices, max_length=length(object$blocks))
  colnames <- get_feature_names(x=object$x)
  indices <- vector()
  blocks <- object$blocks[block_indices]

  for (block in blocks) {
    indices <- c(indices, block@features)
  }

  col_idxs <- match(indices, colnames)

  return(col_idxs)
}

check_indices <- function(indices, max_length) {
  if (length(indices) == 0) {
    err_must_provide_indices()
  }
  if (max(indices) > max_length) {
    err_index_out_of_bounds()
  }
}

err_must_provide_indices <- function() {
  stop(paste("block_indices must have a value.", sep=""))
}

err_index_out_of_bounds <- function() {
  stop(paste("block_indices out of bounds.", sep=""))
}

conditional_matrix <- function(x, indices, drop=TRUE) {
  if (ncol(x) == length(indices)){
    return(FALSE)
  }

  if (drop==TRUE) {
    drop_indices <- indices
    keep_indices <- setdiff(1:ncol(x), drop_indices)
  } else {
    keep_indices <- indices
    drop_indices <- setdiff(1:ncol(x), keep_indices)
  }

  sigma_11 <- x[keep_indices, keep_indices]
  sigma_22 <- x[drop_indices, drop_indices]
  sigma_12 <- x[keep_indices, drop_indices]
  sigma_21 <- x[drop_indices, keep_indices]
  
  sigma_22.1 = sigma_11 + sigma_12 %*% solve(sigma_22) %*% sigma_21

  return(sigma_22.1)
}

select_threshold <- function(
  x,
  c,
  eigen,
  thresholds,
  threshold_mode,
  feature_names,
  check,
  expvar,
  helper) {
  if (length(thresholds) > 1) {
    result <- list()

    for (threshold in thresholds) {
      result[[length(result) + 1]] <- helper(
        x=x,
        c=c,
        eigen=eigen,
        threshold=threshold,
        threshold_mode=threshold_mode,
        feature_names=feature_names,
        check=check,
        expvar=expvar
      )
    }
  } else {
    result <- helper(
      x=x,
      c=c,
      eigen=eigen,
      threshold=thresholds,
      threshold_mode=threshold_mode,
      feature_names=feature_names,
      check=check,
      expvar=expvar
    )
  }

  return(result)
}

pla_helper <- function(
  x,
  c,
  eigen,
  threshold,
  threshold_mode,
  feature_names,
  check,
  expvar) {
  threshold_matrix <- select_thresholding(
    eigen_vectors=eigen$vectors,
    threshold=threshold,
    mode=threshold_mode
  )
  valid_threshold_matrix(threshold_matrix)
  blocks <- get_blocks(
    threshold_matrix=threshold_matrix,
    feature_names=feature_names,
    check=check
  )
  blocks <- calculate_explained_variance(
    blocks=blocks,
    eigen=eigen,
    feature_names=feature_names,
    type=expvar,
    threshold_matrix=threshold_matrix
  )

  result <- list(
    x=x,
    c=c,
    loadings=eigen$vectors,
    threshold=threshold,
    threshold_mode=threshold_mode,
    blocks=blocks
  )
  class(result) <- "pla"

  return(result)
}

spla_helper <- function(
  x,
  c,
  eigen,
  threshold,
  threshold_mode,
  feature_names,
  check,
  expvar) {
  threshold_matrix <- select_thresholding(
    eigen_vectors=eigen$vectors,
    threshold=threshold,
    mode=threshold_mode
  )
  threshold_matrix <- valid_threshold_matrix_spla(threshold_matrix)
  blocks <- get_blocks(
    threshold_matrix=threshold_matrix,
    feature_names=feature_names,
    check=check
  )

  blocks <- calculate_explained_variance(
    blocks=blocks,
    eigen=eigen,
    feature_names=feature_names,
    type=expvar,
    threshold_matrix=threshold_matrix,
    is_absolute=TRUE
  )

  result <- list(
    x=x,
    loadings=eigen$vectors,
    threshold=threshold,
    threshold_mode=threshold_mode,
    blocks=blocks
  )
  class(result) <- "pla"

  return(result)
}

str_loadings <- function(loadings, threshold, threshold_mode, feature_names) {
  loadings <- unclass(loadings)
  rownames(loadings) <- feature_names
  threshold_matrix <- select_thresholding(
    eigen_vectors=loadings,
    threshold=threshold,
    mode=threshold_mode
  )
  strrep <- format(round(loadings, digits=3L))
  nc <- nchar(strrep[1L], type="c")
  strrep[round(abs(threshold_matrix), digits=3L) == 0] <- strrep(" ", nc)

  return(strrep)
}

select_sparse_type_orthogonal <- function(type) {
  value = switch(
    tolower(type),
    "data"=TRUE,
    "dispersion"=FALSE,
    err_wrong_sparse_type(type=type, orthogonal=TRUE)
  )

  return(value)
}

select_sparse_type_not_orthogonal <- function(type) {
  value = switch(
    tolower(type),
    "data"="predictor",
    "dispersion"="Gram",
    err_wrong_sparse_type(type=type, orthogonal=FALSE)
  )

  return(value)
}

err_wrong_sparse_type <- function(type, orthogonal) {
  stop(
    paste(
      "'",
      type,
      "'",
      " is not a valid value for type using",
      " orthogonal = ",
      orthogonal,
      ".",
      sep=""
    )
  )
}