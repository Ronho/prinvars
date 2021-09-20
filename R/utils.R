get_feature_names <- function(x) {
  feature_names <- colnames(x)
  feature_names <- if (length(feature_names) > 0) feature_names else 1:ncol(x)

  return(feature_names)
}

create_block <- function(feature_names, selected_features) {
  if (length(feature_names) > 0) {
    selected_features <- feature_names[selected_features]
  }

  return(new("Block", features=selected_features))
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

  if (length(blocks) > 1) {
    for (block in blocks) {
      indices <- c(indices, block@features)
    }
  } else {
    indices <- c(indices, blocks@features)
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