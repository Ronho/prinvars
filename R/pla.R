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
#' \code{percentage} indicates that the cutoff value is determined by the
#' maximum element of each vector multiplied with the threshold value.
#' @param expvar a character string indicating the method used for calculating
#' the explained variance. \code{approx} uses the explained variance of each
#' eigenvector i.e. its eigenvalue. \code{exact} uses the variance of each
#' variable.
#' @param check a character string indicating if only rows or rows as well as
#' columns are used to detect the underlying block structure. \code{rows} checks
#' if the rows fulfill the required structure. \code{rnc} checks if rows and
#' columns fulfill the required structure.
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
#' See Bauer and Drabant (2021) for more information.
#'
#' @examples
#' if (requireNamespace("AER")) {
#'   require(AER)
#'   data("OECDGrowth")
#'
#'   ## The scales in OECDGrowth differ hence using the correlation matrix is
#'   ## highly recommended.
#'
#'   pla(OECDGrowth, thresholds = 0.5) ## not recommended
#'   pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'
#'   ## We obtain three blocks: (randd), (gdp85, gdp60) and (invest, school,
#'   ## popgrowth). Block 1, i.e. the 1x1 block (randd), explains only 5.76% of
#'   ## the overall variance. Hence, discarding this block seems appropriate.
#'
#'   pla_obj <- pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'   pla.drop_blocks(pla_obj, c(1)) ## drop block 1
#'
#'   ## Sometimes, considering the blocks we keep rather than the blocks we want
#'   ## to discard might be more convenient.
#'
#'   pla.keep_blocks(pla_obj, c(2, 3)) ## keep block 2 and block 3
#' }
#'
#' @references
#' \insertRef{Bauer.2021}{prinvars}
#'
#' @export
pla <- function(x,
                cor = FALSE,
                scaled_ev = FALSE,
                thresholds = 0.33,
                threshold_mode = c("cutoff", "percentage"),
                expvar = c("approx", "exact"),
                check = c("rnc", "rows"),
                ...) {
  chkDots(...)

  threshold_mode <- match.arg(threshold_mode)
  check <- match.arg(check)
  expvar <- match.arg(expvar)

  feature_names <- get_feature_names(x = x)
  x <- scale(x, center = TRUE, scale = FALSE)
  c <- select_cor(x = x, cor = cor)
  eigen <- eigen(as.matrix(c))
  eigen$vectors <- select_eigen_vector_scaling(
    eigen_vectors = eigen$vectors,
    scale = scaled_ev
  )

  result <- select_threshold(
    x = x,
    c = c,
    eigen = eigen,
    thresholds = thresholds,
    threshold_mode = threshold_mode,
    feature_names = feature_names,
    check = check,
    expvar = expvar,
    helper = pla_helper
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
#' if (requireNamespace("AER")) {
#'   require(AER)
#'   data("OECDGrowth")
#'
#'   pla_obj <- pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'   print(pla_obj)
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
      loadings = x$loadings,
      threshold = x$threshold,
      threshold_mode = x$threshold_mode,
      feature_names = feature_names,
      C = x$EC[[x$EC$criterion]],
      criterion = x$EC$criterion
    ),
    quote = FALSE,
    ...
  )

  invisible(x)
}

#' @title Keep Blocks
#'
#' @description Used to pass the indices of the blocks we want to keep (i.e.
#' which we do no want to be discarded).
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
#' if (requireNamespace("AER")) {
#'   require(AER)
#'   data("OECDGrowth")
#'
#'   pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'
#'   ## we obtain three blocks: (randd), (gdp85,gdp60) and (invest, school,
#'   ## popgrowth). Block 1, i.e. the 1x1 block (randd), explains only 5.76% of
#'   ## the overall variance. Hence, discarding this block seems appropriate.
#'   ## Therefore, we keep block 2 and block 3.
#'
#'   pla_obj <- pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'   pla.keep_blocks(pla_obj, c(2, 3)) ## keep block 2 and block 3
#' }
#'
#' @export
pla.keep_blocks <- function(object, blocks, ...) {
  chkDots(...)
  col_idxs <- get_indices(object = object, block_indices = blocks)
  cc_matrix <- conditional_matrix(
    x = object$c,
    indices = col_idxs,
    drop = TRUE
  )
  x <- object$x[, col_idxs, drop = FALSE]

  result <- list(
    x = x,
    cc_matrix = cc_matrix
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
#' if (requireNamespace("AER")) {
#'   require(AER)
#'   data("OECDGrowth")
#'
#'   pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'
#'   ## we obtain three blocks: (randd), (gdp85,gdp60) and (invest, school,
#'   ## popgrowth). Block 1, i.e. the 1x1 block (randd), explains only 5.76% of
#'   ## the overall variance. Hence, discarding this block seems appropriate.
#'
#'   pla_obj <- pla(OECDGrowth, cor = TRUE, thresholds = 0.5)
#'   pla.drop_blocks(pla_obj, c(1)) ## drop block 1
#' }
#'
#' @export
pla.drop_blocks <- function(object, blocks, ...) {
  chkDots(...)
  col_idxs <- get_indices(object = object, block_indices = blocks)
  conditional_matrix <- conditional_matrix(
    x = object$c,
    indices = col_idxs,
    drop = FALSE
  )
  x <- object$x[, -col_idxs, drop = FALSE]

  result <- list(
    x = x,
    conditional_matrix = conditional_matrix
  )

  return(result)
}

#' @title Sparse Principal Loading Analysis
#'
#' @description This function performs sparse principal loading analysis
#' on the given data matrix. We refer to Bauer (2022) for more information.
#' The corresponding sparse loadings are calculated either using \code{PMD} from
#' the \code{PMA} package or using \code{spca} from the \code{elasticnet}
#' package. The respective methods are given by Zou et al. (2006) and Witten et
#' al. (2009) respectively.
#' @param x a numeric matrix or data frame which provides the data for the
#' sparse principal loading analysis.
#' @param method chooses the methods to calculate the sparse loadings.
#' \code{PMD} uses the method introduced by Witten et al. (2009) provided in the \code{PMA} package,
#' \code{rSVD} uses the method introduced by Shen and Huang (2006) provided in the \code{irlba}
#' package, and \code{PP} uses the method introduced by Croux et al. (2013) provided in the \code{rrcovHD}
#' package.
#' @param para an intiger indicating the sparseness. When \code{method="PMD"}: \code{para}
#' gives the bound for the L1 regularization. When \code{method="rSVD"}: \code{para}
#' specifies the number of nonzero components in each loading.
#' @param cor a logical value indicating whether the calculation should use the
#' correlation or the covariance matrix.
#' @param criterion Which of the evaluation criteria (EC) should be used?
#' @param threshold a numeric value used to determine zero elements in the
#' loadings.
#' @param max.iter maximum number of iterations.
#' @param orthogonal a logical value indicating if the sparse loadings are
#' orthogonalized.
#' @param K The number of sparse loadings to be computed. Default is min(N,P)
#' for the methods \code{PMD} and \code{PP}, and min(N,P)-1 for the method \code{rSVD}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' single or list of \code{pla} class containing the following attributes:
#' \item{x}{
#'   a numeric matrix or data frame which equals the input of \code{x}.
#' }
#' \item{EC}{
#'   a numeric vector that contains the evaluation criteria.
#' }
#' \item{loadings}{
#'   a matrix of variable loadings (i.e. a matrix containing the sparse
#'   loadings).
#' }
#' \item{blocks}{
#'   a list of blocks which are identified by sparse principal loading analysis.
#' }
#' \item{W}{
#'   a matrix of variable loadings used to calculate the evaluation criterion.
#'   If \code{criterion="corrected"}, \code{W} contains an orthogonal matrix
#'   with equal weights in the first column of each loading-block. If
#'   \code{criterion="normal"}, \code{W} are the \code{loadings}.
#' }
#'
#' @references
#' \insertRef{Bauer.2022}{prinvars}
#' \insertRef{Witten.2009}{prinvars}
#' \insertRef{Zou.2006}{prinvars}
#'
#' @examples
#' #############
#' ## We replicate a synthetic example similar to Bauer (2022)
#' #############
#'
#' set.seed(1)
#' N <- 500
#' V1 <- rnorm(N, 0, 10)
#' V2 <- rnorm(N, 0, 11)
#'
#' ## Create the blocks (X_1,...,X_4) and (X_5,...,X_8) synthetically
#'
#' X1 <- V1 + rnorm(N, 0, 1) # X_j = V_1 + N(0,1) for j =1,...,4
#' X2 <- V1 + rnorm(N, 0, 1)
#' X3 <- V1 + rnorm(N, 0, 1)
#' X4 <- V1 + rnorm(N, 0, 1)
#'
#' X5 <- V2 + rnorm(N, 0, 1) # X_j = V_1 + N(0,1) for j =5,...9
#' X6 <- V2 + rnorm(N, 0, 1)
#' X7 <- V2 + rnorm(N, 0, 1)
#' X8 <- V2 + rnorm(N, 0, 1)
#'
#' X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8)
#'
#' ## Conduct SPLA to obtain the blocks (X_1,...,X_4) and (X_5,...,X_8)
#'
#' spla(X, para = 1.4, method = "PMD)
#'
#' @export
spla <- function(x,
                 method = c("PMD", "rSVD", "PP"),
                 para,
                 cor,
                 criterion = c("average", "distcor", "RV", "complete"),
                 threshold = 0.05,
                 max.iter,
                 orthogonal,
                 K,
                 Sigma,
                 ...) {
  chkDots(...)
  
  method <- match.arg(method)
  criterion <- match.arg(criterion)

  
  P <- ncol(x)
  N <- nrow(x)

  if(missing(para)){stop("Enter a sparseness parameter")}
  
  if(length(para) != 1) {stop("Enter a single sparseness parameter.")}
  
  if(missing(cor)){cor <- FALSE}
  
  if(missing(max.iter)){
    if(method == "PMD"){ max.iter <- 20}
    if(method == "rSVD"){ max.iter <- 10000}
  } 
  
  if(missing(orthogonal)){orthogonal <- FALSE}
  
  if(missing(K)){
    K <- min(N,P)
    if(method == "rSVD"){K <- min(N,P)-1}
  }else{
    if(K > min(N,P)){
      if("PMD" == method){
        warning("K must be less than or equal to min(nrow(x), ncol(x)).\nMethod continues with K <- min(nrow(x), ncol(x))")
        K <- min(N,P)
      }
      
      if(K >= min(N,P) & any(c("rSVD", "PP") == method)){
        warning("K must be strictly less than min(nrow(x), ncol(x)) when method = rSVD or method = PP.\nMethod continues with K <- min(nrow(x), ncol(x)) - 1")
        K <- min(N,P) - 1
      }
    }
  }
  
  if(missing(Sigma)){
    if(N > P){  
      Sigma <- cov(x)
    } else{  
      Sigma <- scadEst(x, 0.2)  #or 0.1 ??
      #Sigma <- poetEst(dat = x, k = 2, lambda = 0.1)
    }
  } 
  
  x <- scale(x, center = TRUE, scale = cor)
  feature_names <- get_feature_names(x = x)
  eigen <- list()


  switch(method,
    "PMD" = {
      obj <- PMD(
        x = x,
        K = K,
        type = "standard",
        sumabsv = para,
        sumabsu = sqrt(N),
        niter = max.iter,
        trace = FALSE,
        center = FALSE)
   
      eigen$vectors <- obj$v
    },
    "rSVD" = { 
      obj <- ssvd(
        x = x,
        k = K,
        n = para,
        maxit = max.iter)
      
      eigen$vectors <- obj$v
    },
    "PP" = {
      obj <- SPcaGrid(
        x = x,
        k = K,
        kmax = K,
        lambda = para)
      
      eigen$vectors <- obj$loadings
    },
    stop("Method unknown. Select one of PMD, rSVD, or PP")
  )
  

  result <- spla_helper(
    x = x,
    c = c(),
    eigen = eigen,
    threshold = threshold,
    threshold_mode = "cutoff",
    feature_names = feature_names,
    Sigma = Sigma,
    expvar = "approx",
    orthogonal = orthogonal,
    criterion = criterion
  )

  return(result)
}




#' @title Block Detection
#'
#' @description This function returns the blocks for a given loading matrix.
#' @param V a numeric matrix which contains the loadings.
#' @param threshold a numeric value used to determine zero elements in the
#' loadings.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' \item{blocks}{
#'   a list of blocks which are identified by sparse principal loading analysis.
#' }
#'
#' @references
#' \insertRef{Bauer.2022}{prinvars}
#'
#' @examples
#' #############
#'
#' @export
blocks <- function(V,
                 threshold = 0,
                 ...) {
  chkDots(...)
  
  
  if(missing(V)){
    stop("V is required.")
  }
  
  if(length(colnames(V)) == 0){ colnames(V) == as.character(1:ncol(V))}

  threshold_matrx <- get_threshold_matrix(V, threshold)
  
  
  result <- get_blocks(threshold_matrx, colnames(V))
  
  return(result)
}




#' @title Hyperparameter Tuning For Sparse Principal Loading Analysis
#'
#' @description With set.seed(1), the seed for random search can be determined.
#'
#' @export

spla.ht <- function(x,
                    method = c("PMD", "rSVD", "PP"),
                    criterion = c("distcor", "RV", "complete", "average"),
                    EC.min = 0.8,
                    para.lim = c(1, 3),
                    para.steps = 0.1,
                    threshold.lim = c(0, 0.5),
                    threshold.steps = 0.05,
                    K,
                    Sigma,
                    splits = c("multiple", "single"),
                    cor,
                    max.iter,
                    greedy = TRUE,
                    n.cores,
                    ...) {
  chkDots(...)
  
  N <- nrow(x)
  P <- ncol(x)
  
  
  if (EC.min < 0 || EC.min > 1) {
    stop("EC.min not between 0 and 1")
  }
  
  if(missing(cor)){cor = FALSE}
  
  if(missing(K)){
    K <- min(N,P)
    if(method == "rSVD"){K <- min(N,P)-1}
  }else{
    if(K > min(N,P)){
      if("PMD" == method){
        warning("K must be less than or equal to min(nrow(x), ncol(x)).\nMethod continues with K <- min(nrow(x), ncol(x))")
        K <- min(nrow(x), ncol(x))
      }
      
      if(K >= min(N,P) & any(c("rSVD", "PP") == method)){
        warning("K must be strictly less than min(nrow(x), ncol(x)) when method = rSVD or method = PP.\nMethod continues with K <- min(nrow(x), ncol(x)) - 1")
        K <- min(nrow(x), ncol(x)) - 1
      }
    }
  }
  
  para.grid <- seq(para.lim[2], para.lim[1], -para.steps)
  tau.grid <- seq(threshold.lim[1], threshold.lim[2], threshold.steps)
  
  if(missing(n.cores)){n.cores <- length(tau.grid)}
  
  output <- data.frame()

  for(para in para.grid){
    cat("\rpara", para, "till", para.lim[1], "remain.")

    output_list <- mclapply(tau.grid, function(tau) run.spla(tau=tau,x = x, method = method, para = para, cor = cor,
                                                             criterion = criterion,Sigma = Sigma, K = K, EC.min = EC.min),
                            mc.cores = n.cores)
    
    result <- na.omit(as.data.frame(do.call(rbind, output_list)))
    if(!nrow(result) == 0){
      colnames(result) <- c("largest EC", "n.splits", "n.blocks", "para", "threshold", "K")
      output <- rbind(output, result)
      output <- output[order(-output$n.splits, -output$`largest EC`),]
      
      if(greedy){
        cat("\rSplit with para =", output[1,]$para, "and threshold =", output[1,]$threshold, ".\n")
        return(output[1,])
      }
    }
  }

  #cat("\n")
  
  if (length(output)==0) {
    cat(" No (further) blocks identified on this grid.\n")
  } else {
    
    splits <- match.arg(splits)
    switch(splits,
           "multiple" = {
             output <- output[order(-output$n.splits, -output$EC),]
           },
           "single" = { 
             output <- output[order(-output$EC, output$n.blocks),]
           }
    )
    
    rownames(output) <- seq_len(nrow(output))
    return(output)
    
  }
  
}










  

#' @title Complete Split
#'
#' @description With set.seed(1), the seed for random search can be determined.
#'
#' @export

complete.split <- function(x,
                           EC.min = 0.8,
                           method = c("PMD", "rSVD", "PP"),
                           criterion = c("average", "distcor", "RV", "complete"),
                           para.lim = c(1, 2),
                           para.steps = 0.1,
                           threshold.lim = c(0, 0.5),
                           threshold.steps = 0.05,
                           K,
                           Sigma,
                           cor,
                           splits = c("multiple", "single"),
                           greedy = TRUE,
                           n.cores,
                           ...) {
  chkDots(...)

  method <- match.arg(method)
  criterion <- match.arg(criterion)
  splits <- match.arg(splits)
  
  N <- nrow(x)
  P <- ncol(x)
  
  if(length(colnames(x)) == 0){colnames(x) <- as.character(1:P)}
  colnames(Sigma) <- colnames(x)

  if(missing(cor)){cor = FALSE}
  
  if(missing(K)){
    K <- min(N,P)
    if(method == "rSVD"){K <- min(N,P)-1}
  }else{
    if(K > min(N,P)){
      if("PMD" == method){
        warning("K must be less than or equal to min(nrow(x), ncol(x)).\nMethod continues with K <- min(nrow(x), ncol(x))")
        K <- min(N,P)
      }
      
      if(K >= min(N,P) || any(c("rSVD", "PP") == method)){
        warning("K must be strictly less than min(nrow(x), ncol(x)) when method = rSVD or method = PP.\nMethod continues with K <- min(nrow(x), ncol(x)) - 1")
        K <- min(N,P) - 1
      }
    }
  }
  
  if(missing(Sigma)){
    if(N > P){  
      Sigma <- cov(x)
    } else{  
      Sigma <- scadEst(x, 0.2)  #or 0.1 ??
      #Sigma <- poetEst(dat = x, k = 2, lambda = 0.1)
    }
  } 
  
  if(missing(n.cores)){n.cores <- length(seq(threshold.lim[1], threshold.lim[2], threshold.steps))}
  
  sub_matrices <- list(colnames(x))
  results <- list()
  
  while (TRUE) {
    for (i in 1:length(sub_matrices)) {
      
      if (length(sub_matrices) == 0) {
        break  # Stop if no more splitting can be done
      }
      
      if(length(sub_matrices[[1]]) == 1){
        results <- c(results, sub_matrices[[1]])
        sub_matrices <- sub_matrices[-1]
        next
      }
            
      sub_results <- split(x = x[, colnames(x) %in% sub_matrices[[1]] ], EC.min = EC.min,
                           method = method, criterion = criterion,
                           para.lim = para.lim, para.steps = para.steps,
                           threshold.lim = threshold.lim, threshold.steps = threshold.steps,
                           Sigma = Sigma[colnames(Sigma) %in% sub_matrices[[1]],
                                         colnames(Sigma) %in% sub_matrices[[1]] ],
                           splits = splits, greedy = greedy, K = K, n.cores = n.cores) 
      
      
      if(length(sub_matrices) == 1){
        sub_matrices <- list()
      } else{
        sub_matrices <- sub_matrices[-1]
      }
      
      if(length(sub_results) == 1){
        results <- c(results, sub_results)
      }else{
        sub_matrices <- c(sub_matrices, sub_results)
      }
      
    }
    if (length(sub_matrices) == 0) {
      break  # Stop if no more splitting can be done
    }
    
  }
  
  
  return(results)
}

