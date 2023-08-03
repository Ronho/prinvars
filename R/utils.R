get_feature_names <- function(x) {
  feature_names <- colnames(x)
  if (length(feature_names) <= 0) {
    feature_names <- seq_len(ncol(x))
  }

  return(feature_names)
}

create_block <- function(
  feature_names,
  selected_features,
  is_valid,
  ev_influenced) {
  if (length(feature_names) > 0) {
    selected_features <- feature_names[selected_features]
  }

  return(new("Block", features=selected_features, is_valid=is_valid,
    ev_influenced=ev_influenced))
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
  if (ncol(x) == length(indices)) {
    return(FALSE)
  }

  if (drop == TRUE) {
    drop_indices <- indices
    keep_indices <- setdiff(seq_len(ncol(x)), drop_indices)
  } else {
    keep_indices <- indices
    drop_indices <- setdiff(seq_len(ncol(x)), keep_indices)
  }

  sigma_11 <- x[keep_indices, keep_indices]
  sigma_22 <- x[drop_indices, drop_indices]
  sigma_12 <- x[keep_indices, drop_indices]
  sigma_21 <- x[drop_indices, keep_indices]

  sigma_22.1 <- sigma_11 + sigma_12 %*% solve(sigma_22) %*% sigma_21

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
  expvar,
  orthogonal,
  criterion) {

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

  # Change order of rows to follow the block diagonal form.
  feature_idxs <- c()
  ev_idxs <- c()

  for (i in seq_along(blocks)) {
    feature_idxs <- c(feature_idxs, match(blocks[[i]]@features, feature_names))
    ev_idxs <- c(ev_idxs, blocks[[i]]@ev_influenced)
    blocks[[i]]@ev_influenced <- which(ev_idxs %in% blocks[[i]]@ev_influenced,
      arr.ind=TRUE, useNames=FALSE)
  }

  columns <- seq_len(ncol(eigen$vectors))
  I <- diag(1, nrow(eigen$vectors), ncol(eigen$vectors))
  P1 <- I[feature_idxs, ]
  P2 <- I[, ev_idxs]
  x_P1 <- x[, feature_idxs]
  eigen$vectors <- P1 %*% eigen$vectors %*% P2
  threshold_matrix <- P1 %*% threshold_matrix %*% P2
  feature_names <- feature_names[feature_idxs]
  colnames(eigen$vectors) <- sapply(columns[ev_idxs],
    function(num) paste("[,", num, "]", sep=""))

  if (orthogonal) {
    eigen$vectors[threshold_matrix == 0] <- 0
    svd <- svd(eigen$vectors)
    eigen$vectors <- svd$u %*% t(svd$v)
  }

  rownames(eigen$vectors) <- feature_names

 
  

  #### FITTING CRITERIA:
  # Calculate either RV, distcor, complete, or average
  if (criterion == "distcor"){
    #calculate fitting_criteria <- distcor
  }
  if (criterion == "RV"){
    #calculate fitting_criteria <- RV
  }
  if (criterion == "complete"){
    #calculate fitting_criteria <- complete
  }
  if (criterion == "average"){
    #calculate fitting_criteria <- average
  }
  #For now:
  #Replaces by random values to not cause errors
  fitting_criteria <- rep(1,ncol(eigen$vectors))

  

  eigen$values = vector(length = ncol(eigen$vectors))
  x_1 <- x_P1 %*% eigen$vectors[,1] %*% t(eigen$vectors[,1])
  eigen$values[1] <- sum(diag( t(x_1) %*% x_1 ))
  
  for (k in 2:ncol(eigen$vectors)){
    x_k <- x_P1 %*% eigen$vectors[,1:k] %*% t(eigen$vectors[,1:k])
    x_k_1 <- x_P1 %*% eigen$vectors[,1:(k-1)] %*% t(eigen$vectors[,1:(k-1)])
    eigen$values[k] <- sum(diag( t(x_k) %*% x_k )) - sum(diag( t(x_k_1) %*% x_k_1 ))
  }
  
  eigen$values <- eigen$values / sum(diag(t(x)%*%x))
  eigen$var.all = sum(diag(cov(x)))
  
  
  

  blocks <- calculate_explained_variance(
    blocks=blocks,
    eigen=eigen,
    feature_names=feature_names,
    type=expvar,
    threshold_matrix=threshold_matrix,
    is_absolute=TRUE
  )

  
 
  # Only first entry of each block should be depicted since the other entries
  # depend on this one
  for (block in blocks) {
    if (length(block@ev_influenced) > 1) {
      for (i in block@ev_influenced[-1]) {
        fitting_criteria[i] <- 0
      }
    }
  }

  fitting_criteria[1] <- 0 # First entry will not be used either

  result <- list(
    x=x,
    EC=fitting_criteria,
    loadings=eigen$vectors,
    threshold=threshold,
    threshold_mode=threshold_mode,
    blocks=blocks
  )
  class(result) <- "pla"

  return(result)
}

str_loadings <- function(
  loadings,
  threshold,
  threshold_mode,
  feature_names,
  C) {
  loadings <- unclass(loadings)
  threshold_matrix <- select_thresholding(
    eigen_vectors=loadings,
    threshold=threshold,
    mode=threshold_mode
  )

  # This should be TRUE for SPLA!
  if (!is.null(C)) {
    # Add fitting criteria for SPLA to output
    loadings <- rbind(loadings, rep.int(0, ncol(loadings)))
    loadings <- rbind(loadings, C)

    # Prevent threshold_matrix from overwriting new rows by columwraps
    threshold_matrix <- rbind(threshold_matrix, rep(1, ncol(loadings)))
  }

  strrep <- format(round(loadings, digits=3L))
  nc <- nchar(strrep[1L], type="c")

  # This should be TRUE for SPLA!
  if (!is.null(C)) {
    # Add fitting criteria for SPLA to output
    strrep[loadings == 0] <- strrep(" ", nc)
    feature_names <- c(feature_names, " ", "C:")
  } else {
    strrep[threshold_matrix == 0] <- strrep(" ", nc)
  }

  rownames(strrep) <- feature_names

  return(strrep)
}

select_sparse_type_orthogonal <- function(type) {
  value <- switch(
    tolower(type),
    "data"=TRUE,
    "dispersion"=FALSE,
    err_wrong_sparse_type(type=type, orthogonal=TRUE)
  )

  return(value)
}

select_sparse_type_not_orthogonal <- function(type) {
  value <- switch(
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
