#' @include block.R
#' @include utils.R

#' @title Principle Loading Analysis
#'
#' This function performs principle loading analysis on a data set.
#'
#' @param x data.frame matrix, the raw data.
#' @param manipulator character string, indicating the function which is used on the
#' data set; Options: "cov" for covariance, "cor" for correlation; defaults to "cov".
#' @param scaled_ev boolean, indicating whether the eigen vectors should be
#' scaled; defaults to FALSE.
#' @param threshold numeric, determines "small" values inside
#' the eigen vectors; defaults to 0.33.
#' @param expvar character string, determines the method used for calculating the
#' explained variance; Options: "approx" for an approximation, "exact" for the
#' exact calculation; defaults to "approx".
#' @param type character string, determines how eigen vectors are used; Options: "columns"
#' uses eigen vectors as they are, "rows" uses rows of eigen vector matrix, "rnc" requires
#' equality on columns and rows; defaults to "columns"
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' pla class list of the following attributes:
#' \item{eigen_vectors}{matrix, eigen vector matrix obtained from \code{eigen}}
#' \item{threshold}{numeric, threshold}
#' \item{blocks}{list of blocks, blocks describing the explained variance}
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
                manipulator = "cov",
                scaled_ev = FALSE,
                threshold = 0.33,
                expvar = "approx",
                type = "columns", #or "rows" or "rnc"
                ...) {
  colnames <- get_colnames(x)
  x <- as.matrix(x)
  x <- manipulate_matrix(x, manipulator = manipulator)
  eigen <- eigen(x)
  result <- apply_type(eigen, scaled_ev, threshold, colnames, expvar, type)
  return(result)
}

#' @title Principle Loading Analysis - Thresholds
#'
#' This function allows you to perform principle loading analysis on a data set
#' using multiple different thresholds. This allows you to determine the best 
#' threshold on a given data set.
#' @seealso [pla()] for further information about pla.
#'
#' @param x data.frame matrix, the raw data.
#' @param thresholds vector of numeric, determines "small" values inside
#' the eigen vectors.
#' @param ... further arguments passed to pla.
#'
#' @return
#' List of the following parameters which are provided as a pla class list:
#' \item{eigen_vectors}{matrix, eigen vector matrix obtained from \code{eigen}}
#' \item{threshold}{numeric, threshold}
#' \item{blocks}{list of blocks, blocks describing the explained variance}
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9)
#' )
#' thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' pla.thresholds(data, thresholds = thresholds)
#' 
#' @export
pla.thresholds <- function(x, thresholds, ...) {
  results <- list()
  for (threshold in thresholds) {
    results[[length(results) + 1]] <- pla(x, threshold = threshold, ...)
  }
  return(results)
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
  i = 1
  cat("Explained Variances for each block with threshold",
      x$threshold,
      "\n")
  for (block in x$blocks) {
    cat("Block ", i, ": ", str(block), "\n", sep = "")
    i = i + 1
  }
  invisible(x)
}

#' @title Keep Blocks
#'
#' Used to only keep each variable of the original data set which is part of
#' any of the passed blocks.
#'
#' @param x data.frame matrix, the raw data; should be the same used to obtain the blocks.
#' @param blocks list of blocks; blocks that are obtained by the pla function and should
#' be kept.
#'
#' @return data.frame, data without the removed variables
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9),
#' d = c(10:12)
#' )
#' obj <- pla(data)
#' data <- pla.keep_blocks(data, obj$blocks[1])
#' 
#' @export
pla.keep_blocks <- function(x, blocks) {
  colnames <- get_colnames(x)
  indices <- vector()
  for (block in blocks) {
    indices <- c(indices, block@indices)
  }
  col_idxs <- match(indices, colnames)

  x <- x[,col_idxs, drop = FALSE]
  return(x)
}

#' @title Drop Blocks
#'
#' Used to remove each variable from the original data set which is part of
#' any of the passed blocks.
#'
#' @param x data.frame matrix, the raw data; should
#' be the same used to obtain the blocks.
#' @param blocks list of blocks; blocks obtained by the pla function and should be dropped.
#'
#' @return data.frame, data without the removed variables
#'
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9),
#' d = c(10:12)
#' )
#' obj <- pla(data)
#' data <- pla.drop_blocks(data, obj$blocks[1])
#' 
#' @export
pla.drop_blocks <- function(x, blocks) {
  colnames <- get_colnames(x)
  indices <- vector()
  for (block in blocks) {
    indices <- c(indices, block@indices)
  }
  col_idxs <- match(indices, colnames)

  x <- x[,-col_idxs, drop = FALSE]
  return(x)
}
