#' @include block.R
#' @include utils.R

#' @title Principle Loading Analysis
#'
#' This function allows you to perform principle loading analysis on a data set.
#'
#' @param x data.frame matrix, the raw data.
#' @param manipulator character string, indicating the function which is used on the
#' data set; Options: "cov" for covariance, "cor" correlation; defaults to "cov".
#' @param scaled_ev boolean, indicating whether the eigen vectors should be
#' scaled; defaults to FALSE.
#' @param threshold numeric, determines "small" values inside
#' the eigen vectors; defaults to 0.3.
#' @param expvar character string, determines the method used for calculating the
#' explained variance; Options: "approx" for an approximation, "exact" for the
#' exact calculation; defaults to "approx".
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' \item{eigen_vectors}{matrix, eigen vector matrix obtained from \code{eigen}}
#' \item{threshold}{numeric, threshold}
#' \item{blocks}{list of blocks, blocks describing the explained variance}
#' See \url{https://arxiv.org/pdf/2007.05215.pdf},
#' \url{https://arxiv.org/pdf/2102.09912.pdf} for more information.
#' @examples
#' data <- data.frame(
#' a = c(1:3),
#' b = c(4:6),
#' c = c(7:9)
#' )
#' pla(data)
#' @export
pla <- function(x,
                manipulator = "cov",
                scaled_ev = FALSE,
                threshold = 0.3,
                expvar = "approx",
                ...) {
  x <- as.matrix(x)
  x <- .manipulate_matrix(x, manipulator = manipulator)
  eigen <- eigen(x)

  # scale eigen vectors (column wise) between -1 and 1
  if (scaled_ev)
    eigen$vectors = apply(
      eigen$vectors,
      MARGIN = 2,
      FUN = function(x) {
        x / max(abs(x))
      }
    )

  # create threshold matrix: each element is assigned a 0 or 1 depending
  # if it is underneath or above the threshold
  if (threshold > 1 || threshold < 0)
    warning("Threshold should be between 0 and 1.")
  x <- ifelse(abs(eigen$vectors) < threshold, 1, 0)

  blocks <- .create_blocks(x)

  # get explained variance for each block
  blocks <- lapply(blocks, function(block) {
    block@explained_variance <- .explained_variance(eigen$values,
                                                    eigen$vectors,
                                                    block@columns,
                                                    expvar)
    return(block)
  })


  result <- list(
    "eigen_vectors" = eigen$vectors,
    "threshold" = threshold,
    "blocks" = blocks
  )

  class(result) <- "pla"
  return(result)
}

#' @title Principle Loading Analysis - Thresholds
#'
#' This function allows you to perform principle loading analysis on a data set
#' using multiple different thresholds. This allows you to determine
#' the best threshold on a given data set.
#' @seealso [pla()] for further information about pla.
#'
#' @param x data.frame matrix, the raw data.
#' @param thresholds vector of numeric, determines "small" values inside
#' the eigen vectors.
#' @param ... further arguments passed to pla.
#'
#' @return
#' List of the following parameters which are provided as a list:
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
#' thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' pla.thresholds(data, thresholds = thresholds)
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
#' Prints the explained variances of the blocks inside the object.
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
#' Used to remove each variable from the original data set which is part of
#' any of the blocks that have not been excluded by the block_indizes.
#'
#' @param x data.frame matrix, the raw data; should
#' be the same used to obtain the blocks.
#' @param blocks list of blocks; obtained by the pla function.
#' @param block_indizes vector of numeric; indizes of the blocks from the
#' block list which variables should not be removed from the data set.
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
#' data <- pla.keep_blocks(data, obj$blocks, c(2, 3))
#' @export
pla.keep_blocks <- function(x, blocks, block_indizes) {
  blocks <- blocks[-block_indizes]
  block_indizes <- 1:length(blocks)
  return(pla.drop_blocks(x, blocks, block_indizes))
}
#' @title Drop Blocks
#'
#' Used to remove each variable from the original data set which is part of
#' any of the specified blocks.
#'
#' @param x data.frame matrix, the raw data; should
#' be the same used to obtain the blocks.
#' @param blocks list of blocks; obtained by the pla function.
#' @param block_indizes vector of numeric, indizes of the blocks from the
#' block list which variables should be removed from the data set.
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
#' data <- pla.drop_blocks(data, obj$blocks, c(1))
#' @export
pla.drop_blocks <- function(x, blocks, block_indizes) {
  indizes <- vector()
  for (idx in block_indizes) {
    indizes <- c(indizes, blocks[[idx]]@columns)
  }

  x <- x[,-indizes, drop = FALSE]
  return(x)
}
