#' @include block.R
#' @include utils.R
#' @include explained-variance.R
#' @include get-blocks.R
#' @include scale.R
#' @include thresholding.R
#' @include cor.R

#' @title Principal Loading Analysis
#'
#' @description This function performs a principal loading analysis on the given
#' data matrix and returns the results as an object of class \code{pla}.
#'
#' @param x a numeric matrix or data frame which provides the data for the
#' principal loading analysis.
#' @param cor a logical value indicating whether the calculation should use the
#' correlation or the covariance matrix.
#' @param scaled_ev a logical value indicating whether the eigenvectors should
#' be scaled.
#' @param thresholds a numeric value or list of numeric values used to determine
#' "small" values inside the eigenvectors. If multiple values are given, a list
#' of pla results will be returned.
#' @param threshold_mode a character string indicating how the threshold is
#' determined and used. \code{cutoff} indicates that the threshold value is used as
#' a general maximum for all elements. \codepercentage} indicates that the cutoff
#' value is determined by the maximum element of each vector multiplied with the
#' threshold value. The default is set to \code{cutoff}.
#' @param expvar a character string indicating the method used for calculating
#' the explained variance. \code{approx} uses the explained variances of each
#' eigenvectors i.e. their eigenvalues. \code{exact} uses the variance of each variable.
#' @param check a character string indicating if only rows, or if rows as well as columns
#' are used to detect the underlying block structure. \code{rows} checks if the rows fullfill
#' the required structure. \code{rnc} checks if rows and columns fullfill the required structure.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' single or list of pla class containing the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input of \code{x}.
#' }
#' \item{c}{
#'   a numeric matrix or data frame which is the covariance or correlation
#'   matrix based on the input of \code{cov}.
#' }
#' \item{loadings}{
#'   a matrix of variable loadings (i.e. a matrix containing the
#'   eigenvectors of the dispersion matrix).
#' }
#' \item{threshold}{
#'   a numeric value which equals the input of \code{thresholds}.
#' }
#' \item{threshold_mode}{
#'   a character string which equals the input of \code{threshold_mode}.
#' }
#' \item{blocks}{
#'   a list of blocks which are identified by principal loading analysis.
#' }
#' See \url{https://www.sciencedirect.com/science/article/pii/S0047259X21000324} and
#' \url{https://dl.acm.org/doi/10.1145/3475827.3475832} for more information.
#' 
#' @examples
#' require(AER)
#' data("OECDGrowth")
#'
#' ## the scales in OECDGrowth differ hence using the
#' ## correlation matrix is highly recommended
#'
#' pla(OECDGrowth,thresholds = 0.5) ## not recommended
#' pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#'
#' ## we obtain three blocks: (randd), (gdp85,gdp60) and 
#' ## (invest, school, popgrowth). Block 1, i.e. the 1x1 block 
#' ## (randd), explains only 5.76% of the overall variance.
#' ## Hence discarding this block seems appropriate.
#'
#' pla_obj = pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#' pla.drop_blocks(pla_obj, c(1)) ## drop block 1
#'
#' ## Sometimes, considering the blocks we keep rather than
#' ## the blocks we want to discard might be more convenient.
#'
#' pla.keep_blocks(pla_obj, c(2,3)) ## keep block 2 and block 3
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
    if (block@is_valid) {
      cat("Block ", i, ": ", str(block), "\n", sep = "")
    } else {
      cat(str(block), "\n")
    }
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
