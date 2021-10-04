#' @include block.R
#' @include utils.R
#' @include explained-variance.R
#' @include get-blocks.R
#' @include scale.R
#' @include thresholding.R
#' @include cor.R

#' @title Principle Loading Analysis
#'
#' This function performs principle loading analysis on a data set.
#'
#' @param x a numeric matrix or data frame which provides the data for the
#' principal loading analysis.
#' @param cor a logical value indicating whether the calculation should use the
#' correlation or the covariance matrix. The default is set to "TRUE" indicating
#' the use of the correlation matrix.
#' @param scaled_ev a logical value indicating whether the eigenvectors should
#' be scaled. The default is set to "FALSE".
#' @param thresholds a numeric or list of numeric used to determine "small"
#' values inside the eigenvectors. If multiple values are given, a list of pla
#' results will be returned. The default is set to 0.33.
#' @param threshold_mode a character string indicating how the threshold is
#' determined and used. "cutoff" indicates that the threshold value is used as
#' a general maximum for all elements. "percentage" indicates that the cutoff
#' value is determined by the maximum element of each vector multiplied with the
#' threshold value. The default is set to "cutoff".
#' @param expvar a character string indicating the method used for calculating
#' the explained variance."approx" indicates the use of an approximation.
#' "exact" indicates the exact calculation. The default is set to "approx".
#' @param check a character string indicating if only rows or rows and columns
#' are used."rows" checks if the rows fullfill the required structure. "rnc"
#' checks if rows and columns fullfill the required structure. The default is
#' set to "rnc".
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' single or list of pla class containing the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input of 'x'.
#' }
#' \item{c}{
#'   a numeric matrix or data frame which is the covariance or correlation
#'   matrix based on the input of 'cov'.
#' }
#' \item{loadings}{
#'   a matrix of variable loadings (i.e., a matrix whose columns contain the
#'   eigenvectors).
#' }
#' \item{threshold}{
#'   a numeric value which equals the input of thresholds at the corresponding
#'   position.
#' }
#' \item{threshold_mode}{
#'   a character string which equals the input of threshold_mode.
#' }
#' \item{blocks}{
#'   a list of blocks which are identified through the principle loading
#'   analysis.
#' }
#' See \url{https://arxiv.org/pdf/2007.05215.pdf},
#' \url{https://arxiv.org/pdf/2102.09912.pdf} for more information.
#' 
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9)
#' )
#' pla(data)
#' 
#' @export
pla <- function(x,
                cor = FALSE,
                scaled_ev = FALSE,
                thresholds = 0.33,
                threshold_mode = "cutoff",
                expvar = "approx",
                check = "rnc",
                ...) {
  chkDots(...)
  feature_names <- get_feature_names(x=x)
  c <- select_cor(x=x, cor=cor)
  eigen <- eigen(as.matrix(c))
  eigen$vectors <- select_eigen_vector_scaling(
    eigen_vectors=eigen$vectors,
    scale=scaled_ev
  )

  result <- select_threshold(
    x=x,
    c=c,
    eigen=eigen,
    thresholds=thresholds,
    threshold_mode=threshold_mode,
    feature_names=feature_names,
    check=check,
    expvar=expvar
  )

  return(result)
}

#' @title Print Function for pla S3
#'
#' Prints the blocks, threshold, threshold_mode and the loadings.
#'
#' @param x a pla object.
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9)
#' )
#' obj <- pla(data)
#' print(obj)
#' 
#' @export
print.pla <- function(x, ...) {
  chkDots(...)

  cat(
    "Explained Variances for each block with threshold",
    x$threshold,
    "and mode",
    x$threshold_mode,
    ":\n"
  )
  
  i = 1
  for (block in x$blocks) {
    cat("Block ", i, ": ", str(block), "\n", sep = "")
    i = i + 1
  }

  cat("\nLoadings:\n")
  print(
    str_loadings(
      loadings=x$loadings,
      threshold=x$threshold,
      threshold_mode=x$threshold_mode,
      feature_names=get_feature_names(x$x)
    ),
    quote=FALSE,
    ...
  )

  invisible(x)
}

#' @title Keep Blocks
#'
#' Used to only keep each variable of the original data set which is part of any
#' of the blocks according to the passed indices.
#'
#' @param object a pla object.
#' @param blocks a list of numeric values indicating the indices of the blocks
#' that should be kept.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input data for the pla
#'   object without any feature that is not part of the blocks that should be
#'   kept.
#' }
#' \item{cc_matrix}{
#'   a numeric matrix or data frame which contains either the conditional
#'   covariance or correlation matrix.
#' }
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9),
#' d = c(10:12)
#' )
#' obj <- pla(data)
#' data <- pla.keep_blocks(obj, c(1))
#' 
#' @export
pla.keep_blocks <- function(object, blocks, ...) {
  chkDots(...)
  col_idxs <- get_indices(object=object, block_indices=blocks)
  cc_matrix <- conditional_matrix(
    x=object$c,
    indices=col_idxs,
    drop=TRUE
  )
  x <- object$x[,col_idxs, drop = FALSE]

  result <- list(
    x=x,
    cc_matrix=cc_matrix
  )

  return(result)
}

#' @title Drop Blocks
#'
#' Used to remove each variable from the original data set which is part of any
#' of the blocks according to the passed indices.
#'
#' @param object a pla object.
#' @param blocks a list of numeric values indicating the indices of the blocks
#' that should be removed.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input data for the pla
#'   object without any feature that is not part of the blocks that should be
#'   removed.
#' }
#' \item{cc_matrix}{
#'   a numeric matrix or data frame which contains either the conditional
#'   covariance or correlation matrix.
#' }
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9),
#' d = c(10:12)
#' )
#' obj <- pla(data)
#' data <- pla.drop_blocks(obj, c(1))
#' 
#' @export
pla.drop_blocks <- function(object, blocks, ...) {
  chkDots(...)
  col_idxs <- get_indices(object=object, block_indices=blocks)
  conditional_matrix <- conditional_matrix(
    x=object$c,
    indices=col_idxs,
    drop=FALSE
  )
  x <- object$x[,-col_idxs, drop = FALSE]

  result <- list(
    x=x,
    conditional_matrix=conditional_matrix
  )

  return(result)
}
