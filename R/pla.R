#' @include block.R
#' @include utils.R
#' @include explained-variance.R
#' @include get-blocks.R
#' @include scale.R
#' @include thresholding.R
#' @include cov.R

#' @title Principle Loading Analysis
#'
#' This function performs principle loading analysis on a data set.
#'
#' @param x data.frame matrix, the raw data.
#' @param cov boolean, indicating the function which is used on
#' the data set; Options: TRUE for covariance, FALSE for correlation; defaults
#' to TRUE.
#' @param scaled_ev boolean, indicating whether the eigen vectors should be
#' scaled; defaults to FALSE.
#' @param thresholds numeric or list of numeric, used to determine "small"
#' values inside the eigenvectors; if multiple values are given, a list of pla
#' results will be returned; defaults to 0.33.
#' @param threshold_mode character string, indicating how the threshold is
#' determined and used; Options: "cutoff" the threshold value is used as an
#' maximum for all elements, "percentage" the cutoff value is determined by the
#' maximum element of each vector multiplied with the threshold value;
#' defaults to "cutoff".
#' @param expvar character string, determines the method used for calculating
#' the explained variance; Options: "approx" for an approximation, "exact" for
#' the exact calculation; defaults to "approx".
#' @param check character string, determines how eigen vectors are used;
#' Options: "rows" checks if the rows fullfill the required structure, "rnc"
#' checks if rows and columns fullfill the required structure; defaults to
#' "rnc".
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' pla class list of the following attributes:
#' \item{x}{
#'   matrix, transformed matrix
#' }
#' \item{eigen_vectors}{
#'   matrix, eigen vector matrix obtained from \code{eigen} and scaled if
#'   scaled_ev is TRUE
#' }
#' \item{thresholds}{
#'   numeric or list of numeric, equals input of thresholds
#' }
#' \item{threshold_mode}{
#'   numeric, equals input of threshold_mode
#' }
#' \item{blocks}{
#'   list of blocks, blocks describing the explained variance
#' }
#' \item{cov}{
#'   booleans, equals input of cov
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
                cov = TRUE,
                scaled_ev = FALSE,
                thresholds = 0.33,
                threshold_mode = "cutoff",
                expvar = "approx",
                check = "rnc",
                ...) {
  chkDots(...)
  feature_names <- get_feature_names(x=x)
  c <- select_cov(x=x, cov=cov)
  eigen <- eigen(as.matrix(c))
  eigen$vectors <- select_eigen_vector_scaling(
    eigen_vectors=eigen$vectors,
    scale=scaled_ev
  )

  result <- select_threshold(
    x=x,
    c=c,
    cov=cov,
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
#' Prints the blocks and threshold inside the object.
#'
#' @param x pla object.
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
  i = 1

  cat(
    "Explained Variances for each block with threshold",
    x$threshold,
    "and mode",
    x$threshold_mode,
    "\n"
  )

  for (block in x$blocks) {
    cat("Block ", i, ": ", str(block), "\n", sep = "")
    i = i + 1
  }

  invisible(x)
}

#' @title Keep Blocks
#'
#' Used to only keep each variable of the original data set which is part of any
#' of the blocks according to the passed indices.
#'
#' @param object data.frame matrix, the raw data; should be the same used to
#' obtain the blocks.
#' @param blocks list of numeric; indices of blocks that should be kept.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   matrix, transformed matrix with removed attributes
#' }
#' \item{conditional_matrix}{
#'   matrix, conditional matrix
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
#' @param object data.frame matrix, the raw data; should be the same used to
#' obtain the blocks.
#' @param blocks list of numeric; indices of blocks that should be
#' dropped.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   matrix, transformed matrix with removed attributes
#' }
#' \item{conditional_matrix}{
#'   matrix, conditional matrix
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
