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

  return(result)
}

get_threshold_matrix <- function(
    loadings,
    threshold){
  
  loadings[which(abs(loadings) <= threshold)] = 0
  loadings[which(abs(loadings) != 0)] = 1
  
  return(loadings)
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
  
  threshold_matrix <- get_threshold_matrix(
    loadings=eigen$vectors,
    threshold=threshold
  )

  blocks <- get_blocks(
    threshold_matrix=threshold_matrix,
    feature_names=feature_names
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
  Sigma,
  expvar,
  orthogonal,
  criterion) {
 

  threshold_matrix <- get_threshold_matrix(
    loadings=eigen$vectors,
    threshold=threshold
  )

  
  blocks <- get_blocks(
    threshold_matrix=threshold_matrix,
    feature_names=feature_names
  )
  
  if(length(blocks) == 1){
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    cat("No different blocks identified with these parameters.")
    stop()
  }


  # Change order of rows to follow the block diagonal form.
  feature_idxs <- c()
  ev_idxs <- c()

  for (i in seq_along(blocks)) {
    feature_idxs <- c(feature_idxs, match(blocks[[i]]@features, feature_names))
    ev_idxs <- c(ev_idxs, blocks[[i]]@ev_influenced)
    blocks[[i]]@ev_influenced <- which(ev_idxs %in% blocks[[i]]@ev_influenced,
      arr.ind=TRUE, useNames=FALSE)
  }
  
  x_P1 <- x[, feature_idxs]
  eigen$vectors <- eigen$vectors[feature_idxs, ev_idxs]
  threshold_matrix <- threshold_matrix[feature_idxs, ev_idxs] 
  
  feature_names <- feature_names[feature_idxs]
  colnames(eigen$vectors) <- as.character(ev_idxs)
  #colnames(eigen$vectors) <- sapply( (seq_len(ncol(eigen$vectors)))[ev_idxs],
  #  function(num) paste("[,", num, "]", sep=""))

  if (orthogonal) {
    eigen$vectors[threshold_matrix == 0] <- 0
    svd <- svd(eigen$vectors)
    eigen$vectors <- svd$u %*% t(svd$v)
  }

  rownames(eigen$vectors) <- feature_names

  Sigma_P1 <- Sigma[feature_idxs, feature_idxs]

  
  fitting_criteria <- list(
    "criterion"=criterion,
    "distcor"=numeric(length(ev_idxs)),
    "complete"=numeric(length(ev_idxs)),
    "average"=numeric(length(ev_idxs)),
    "rv"=numeric(length(ev_idxs))
  )
  
  
  for (block in blocks) {
    ev <- block@ev_influenced[1]
    feature_idxs <- match(block@features, feature_names)
    non_feature_idxs <- seq_along(feature_names)[-feature_idxs]
    

    
    fitting_criteria$distcor[ev] <- 1 - dcor(x_P1[, feature_idxs], x_P1[, non_feature_idxs])
    fitting_criteria$complete[ev] <- 1 - max(abs(Sigma_P1[feature_idxs, non_feature_idxs]))
    #fitting_criteria$average[ev] <- 1 - mean(as.vector(abs(Sigma_P1[feature_idxs, non_feature_idxs])))
    fitting_criteria$rv[ev] <- 1 - sum(diag((Sigma_P1[feature_idxs, non_feature_idxs] %*% t(Sigma_P1[feature_idxs, non_feature_idxs])))) /
      sqrt(sum(diag((Sigma_P1[feature_idxs, feature_idxs] %*% t(Sigma_P1[feature_idxs, feature_idxs])))) * sum(diag((Sigma_P1[non_feature_idxs, non_feature_idxs] %*% t(Sigma_P1[non_feature_idxs, non_feature_idxs])))))
    
    sorted.E <- sort(as.vector(abs(Sigma_P1[feature_idxs, non_feature_idxs])), decreasing = TRUE)
    n.mean <- min(length(feature_idxs), 5) #WHY 5? TUNING PARAEMTER - CHECK USING SIMULATIONS (PERHAPS 10?)
    fitting_criteria$average[ev] <- 1 - mean(sorted.E[1:n.mean])
    #1*length if considering only E upper right or lower left
    #2*length if considering both E
    
    
  } 
  
  

  eigen$values = vector(length = ncol(eigen$vectors))
  x_1 <- x_P1 %*% eigen$vectors[,1] %*% t(eigen$vectors[, 1])
  eigen$values[1] <- sum(diag(t(x_1) %*% x_1))
  
  for (k in 2:ncol(eigen$vectors)){
    x_k <- x_P1 %*% eigen$vectors[,1:k] %*% t(eigen$vectors[, 1:k])
    x_k_1 <- x_P1 %*% eigen$vectors[,1:(k-1)] %*% t(eigen$vectors[,1:(k - 1)])
    eigen$values[k] <- sum(diag( t(x_k) %*% x_k )) - sum(diag( t(x_k_1) %*% x_k_1 ))
  }

  eigen$values <- eigen$values / sum(diag(t(x) %*% x))
  eigen$var.all <- sum(diag(cov(x)))

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
  C,
  criterion) {
  loadings <- unclass(loadings)
  threshold_matrix <- get_threshold_matrix(
    loadings=loadings,
    threshold=threshold
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
    feature_names <- c(feature_names, " ", paste0(criterion, ":"))
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



run.spla <- function(tau, x, method, para, criterion, Sigma, cor, K, EC.min){
  
  test_error <- try(
    {
      invisible(capture.output(suppressWarnings({
        spla <- spla(x = x,
                     method = method,
                     para = para,
                     cor = cor,
                     criterion = criterion,
                     threshold = tau,
                     orthogonal = FALSE,
                     Sigma = Sigma,
                     K = K)
      })))
    },
    silent = TRUE
  )

  if(inherits(test_error, "try-error")) return(NA)
  
  n.blocks <- length(spla$blocks)
  if (n.blocks == 1) return(NA)
  
  EC <- spla[["EC"]][[criterion]]
  EC <- EC[EC != 0]
  if (max(EC) < EC.min) return(NA)
  
  if(missing(K)){K = ncol(spla$loadings)}
  
  return(c(max(EC), length(which(EC >= EC.min)), n.blocks, para, tau, K))
}






split <- function(x, EC.min, method, criterion, Sigma,
                  para.lim, para.steps, threshold.lim, threshold.steps, splits, greedy, K, n.cores){
  
  if(length(ncol(x)) == 0) return(list(x))
  
  para.lim[2] = min(para.lim[2], sqrt(ncol(x)), sqrt(nrow(x)))
  
  if(para.lim[2] < para.lim[1]){para.lim[1] <- para.lim[2]}
  
  ht <- spla.ht(x = x, method = method, criterion = criterion,
                EC.min = EC.min, para.lim = c(para.lim[1], para.lim[2]), para.steps = para.steps,
                threshold.lim = c(threshold.lim[1], threshold.lim[2]), threshold.steps = threshold.steps,
                splits = splits, Sigma = Sigma, greedy = greedy, K = K, n.cores = n.cores)
  
  
  if(length(ht)  == 0) return(list(colnames(x)))
  
  spla.obj <- spla(x = x, method = method, criterion = criterion,
                   para = ht$para[1], threshold = ht$threshold[1],
                   Sigma = Sigma, K = K)
  
  EC <- spla.obj[["EC"]][[criterion]]
  EC <- EC[EC!=0]
  
  if(splits == "multiple"){
    split.blocks <- which(EC >= EC.min)
    all.split.block.labels <- c()
    result.split <- list()
    
    cat("\n===== SPLITTING", length(split.blocks), "BLOCK(S) ===== with para =", ht$para[1], "and tau =", ht$threshold[1], "\n")
    
    counter <- 1
    for(split.block in split.blocks){
      split.block.labels <- spla.obj[["blocks"]][[split.block]]@features
      all.split.block.labels <- c(all.split.block.labels, split.block.labels)
      result.split[[counter]] <- colnames(x)[(colnames(x) %in% split.block.labels)]
      counter <- counter + 1
    }
    
    if(!all(colnames(x) %in% all.split.block.labels)){
      result.split[[counter]] <- colnames(x)[!(colnames(x) %in% all.split.block.labels)]
    }
    
    return(result.split)
  }
  
  if(splits == "single"){
    split.block  <- which(EC == max(EC))[1]
    split.labels <- spla.obj[["blocks"]][[split.block]]@features
    result.split <- list()
    
    cat("\n===== SPLITTING ===== with para =", ht$para[1], "and tau =", ht$threshold[1], "\n")
    
    result.split[[1]] <- colnames(x)[(colnames(x) %in% split.block.labels)]
    result.split[[2]] <- colnames(x)[!(colnames(x) %in% split.block.labels)]
    
    return(result.split)
  }
  
}


