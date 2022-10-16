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
#' data matrix.
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
#' determined and used. \code{cutoff} indicates the usage of a threshold value.
#' \code{percentage} indicates that the cutoff value is determined by the maximum
#' element of each vector multiplied with the threshold value.
#' @param expvar a character string indicating the method used for calculating
#' the explained variance. \code{approx} uses the explained variance of each
#' eigenvector i.e. its eigenvalue. \code{exact} uses the variance of each variable.
#' @param check a character string indicating if only rows or rows as well as columns
#' are used to detect the underlying block structure. \code{rows} checks if the rows fulfill
#' the required structure. \code{rnc} checks if rows and columns fulfill the required structure.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' single or list of \code{pla} class containing the following attributes:
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
#' See \insertRef{Bauer.2021}{prinvars} for more information.
#' 
#' @examples
#' if(requireNamespace("AER")){
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
#' ## Hence, discarding this block seems appropriate.
#'
#' pla_obj = pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#' pla.drop_blocks(pla_obj, c(1)) ## drop block 1
#'
#' ## Sometimes, considering the blocks we keep rather than
#' ## the blocks we want to discard might be more convenient.
#'
#' pla.keep_blocks(pla_obj, c(2,3)) ## keep block 2 and block 3
#' }
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
  x <- scale(x, center = TRUE, scale = FALSE)
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
    expvar=expvar,
    helper=pla_helper
  )

  return(result)
}

#' @title Print Function for pla S3
#'
#' @description Prints the blocks, threshold, threshold_mode and the loadings.
#'
#' @param x a pla object.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return
#' A pla object which equals the input of \code{x}.
#'
#' @examples
#' if(requireNamespace("AER")){
#' require(AER)
#' data("OECDGrowth")
#'
#' pla_obj = pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#' print(pla_obj)
#' }
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
  
  sum_expvar <- 0
  i <- 1
  for (block in x$blocks) {
    if (block@is_valid) {
      cat("Block ", i, ": ", str(block), "\n", sep = "")
    } else {
      cat(str(block), "\n")
    }

    sum_expvar <- sum_expvar + block@explained_variance
    i <- i + 1
  }

  cat(
    "\nAll blocks together explain ", 
    round(sum_expvar * 100, 2),
    "% of the total variance.\n",
    sep = ""
  )

  feature_names <- rownames(x$loadings)
  if (is.null(feature_names)) {
    feature_names <- get_feature_names(x$x)
  }

  cat("\nLoadings:\n")
  print(
    str_loadings(
      loadings=x$loadings,
      threshold=x$threshold,
      threshold_mode=x$threshold_mode,
      feature_names=feature_names,
      C=x$C
    ),
    quote=FALSE,
    ...
  )

  invisible(x)
}

#' @title Keep Blocks
#'
#' @description Used to pass the indices of the blocks we want to keep 
#' (i.e. which we do no want to be discarded).
#'
#' @param object a \code{pla} object.
#' @param blocks a list of numeric values indicating the indices of the blocks
#' that should be kept.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame containing the reduced set of original
#'   variables.
#' }
#' \item{cc_matrix}{
#'   a numeric matrix or data frame which contains the conditional dispersion
#'   matrix. Depending on the pla procedure, this is either the conditional
#'   covariance matrix or the conditional correlation matrix.
#' }
#'
#' @examples
#' if(requireNamespace("AER")){
#' require(AER)
#' data("OECDGrowth")
#'
#' pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#'
#' ## we obtain three blocks: (randd), (gdp85,gdp60) and 
#' ## (invest, school, popgrowth). Block 1, i.e. the 1x1 block 
#' ## (randd), explains only 5.76% of the overall variance.
#' ## Hence, discarding this block seems appropriate. Therefore,
#' ## we keep block 2 and block 3
#' 
#' pla_obj = pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#' pla.keep_blocks(pla_obj, c(2,3)) ## keep block 2 and block 3
#' }
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
#' @description Used to pass the indices of the blocks we want to discard.
#'
#' @param object a pla object.
#' @param blocks a list of numeric values indicating the indices of the blocks
#' that should be removed.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' list of the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame containing the reduced set of original
#'   variables.
#' }
#' \item{cc_matrix}{
#'   a numeric matrix or data frame which contains the conditional dispersion
#'   matrix. Depending on the pla procedure, this is either the conditional
#'   covariance matrix or the conditional correlation matrix.
#' }
#'
#' @examples
#' if(requireNamespace("AER")){
#' require(AER)
#' data("OECDGrowth")
#'
#' pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#'
#' ## we obtain three blocks: (randd), (gdp85,gdp60) and 
#' ## (invest, school, popgrowth). Block 1, i.e. the 1x1 block 
#' ## (randd), explains only 5.76% of the overall variance.
#' ## Hence, discarding this block seems appropriate.
#'
#' pla_obj = pla(OECDGrowth,cor=TRUE,thresholds = 0.5)
#' pla.drop_blocks(pla_obj, c(1)) ## drop block 1
#' }
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

#' @title Sparse Principal Loading Analysis
#'
#' @description This function performs a sparse principal loading analysis on the given
#' data matrix.
#'
#' @param x a numeric matrix or data frame which provides the data for the
#' principal loading analysis.
#' @param para a vector containing the penalization parameter for each variable
#' @param cor a logical value indicating whether the calculation should use the
#' correlation or the covariance matrix.
#' @param criterion a character string indicating if the weight-corrected evaluation criterion
#' (CEC) or the evaluation criterion (EC) is used. \code{corrected} changes the loadings to
#' weight all variables equally while \code{normal} does not change the loadings and therefore
#' weights the variables differently.
#' @param rho quadratic penalty parameter. We refer to REF(Zou2006) and REF(Bauer22) for
#' a more elaborate explanation.
#' @param max.iter maximum number of iterations.
#' @param trace a logical value indicating if the progress is printed.
#' @param eps.conv a numerical value as convergence criterion.
#' @param orthogonal a logical value indicating if the sparse loadings are orthogonalized.
#' @param check a character string indicating if only rows or rows as well as columns
#' are used to detect the underlying block structure. \code{rows} checks if the rows fulfill
#' the required structure. \code{rnc} checks if rows and columns fulfill the required structure.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' single or list of \code{pla} class containing the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input of \code{x}.
#' }
#' \item{EC}{
#'   a numeric vector that contains the weight-corrected evaluation criterion (CEC) if
#'   \code{criterion = "corrected"} and the evaluation criterion (EC) if \code{criterion = "normal"}.
#' }
#' \item{loadings}{
#'   a matrix of variable loadings (i.e. a matrix containing the sparse loadings).
#' }
#' \item{blocks}{
#'   a list of blocks which are identified by sparse principal loading analysis.
#' }
#' \item{W}{
#'   a matrix of variable loadings used to calculate the evaluation criterion. If \code{criterion = "corrected"},
#'   \code{W} contains an orthogonal matrix with equal weights in the first column of each loading-block. If
#'   \code{criterion = "normal"}, \code{W} equals \code{loadings}.
#' }
#' See REF(Bauer22) for more information.
#'
#' @examples
#' spla(USArrests, para = c(0.5,0.5,0.5,0.5), cor = TRUE)
#'
#' ## we obtain two blocks:
#' ## 1x1 (Urbanpop) and 3x3 (Murder, Aussault, Rape).
#' ## The large EC indicates that the given structure is reasonable.
#'
#' spla(USArrests, para = c(0.5,0.5,0.7,0.5), cor = TRUE)
#'
#' ## we obtain three blocks:
#' ## 1x1 (Urbanpop), 1x1 (Rape) and 2x2 (Murder, Aussault).
#' ## The mid-ish EC for (Murder, Aussault) indicates that the
#' ## found structure might not be adequate.
#'
#' @export
spla <- function(x,
                 para,
                 cor = FALSE,
                 criterion = "corrected",
                 rho = 1e-06,
                 max.iter = 200,
                 trace = FALSE,
                 eps.conv = 1e-3,
                 orthogonal = FALSE,
                 check = "rnc",
                 ...) {
  chkDots(...)

  feature_names <- get_feature_names(x = x)
  eigen <- list()
  x <- scale(x, center = TRUE, scale = FALSE)

  obj <- spca(
    x = x,
    K = ncol(x),
    para = para,
    type = "predictor",
    lambda = rho,
    use.corr = cor,
    max.iter = max.iter,
    trace = trace,
    eps.conv = eps.conv
  )

  eigen$vectors <- obj$loadings
  eigen$values <- obj$pev
  eigen$var.all <- obj$var.all

  result <- spla_helper(
    x = x,
    c = c(),
    eigen = eigen,
    threshold = 0,
    threshold_mode = "cutoff",
    feature_names = feature_names,
    check = check,
    expvar = "approx",
    orthogonal = orthogonal,
    criterion = criterion
  )

  return(result)
}
