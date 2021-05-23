.manipulate_matrix <- function(x, manipulator = "cov") {
  if (manipulator == "cor") return(cor(x, method = c("pearson")))
  if (manipulator == 'cov') return(cov(x, method = c("pearson")))
  stop(paste("'", manipulator, "'", " is not a valid value for manipulator.", sep = ""))
}


.create_blocks <- function(x) {
  blocks <- list()
  violated <- c()
  for (idx in 1:nrow(x)) {
    columns <- .get_row_structure(x[idx, ])
    block_idx <- .get_block_index(blocks, columns)
    violated <- unique(c(violated, block_idx$violates))

    if (block_idx$equals == -1) {
      block <- new("Block", rows = c(idx), columns = columns)
      blocks[[length(blocks)+1]] <- block
    } else {
      blocks[[block_idx$equals]]@rows <- c(blocks[[block_idx$equals]]@rows, idx)
    }
  }

  # remove blocks that violate the n x n structure
  uncontained_rows = c()
  if (length(blocks) > 0) {
    idx <- 1
    for (i in 1:length(blocks)) {
      block <- blocks[[idx]]
      if ((length(block@rows) != length(block@columns)) || (length(intersect(block@columns, violated)) > 0)) {
        uncontained_rows <- unique(c(uncontained_rows, block@rows))
        blocks[[idx]] <- NULL
      } else {
        idx <- idx+1
      }
    }
  }

  # add a new 1x1 block for each row that is not contained within any block
  for (row_nr in uncontained_rows) {
    blocks[[length(blocks)+1]] <- new("Block", rows = c(row_nr), columns = c(row_nr))
  }

  return(blocks)
}

# returns the index of the block containing the exact columns
# returns -1 if no block was found
# returns the column indexes of blocks that can be confirmed to violate the block structure based on the new columns
# (e.g. one block contains columns 2 and 8 --> if this functions retrieves c(8,9), the required structure cannot be given since 9 is not part of c(2,8))
.get_block_index <- function(blocks, columns) {
  violates <- vector() # indexes of blocks that have a confirmed violated structure
  equals <- -1 # index of block that has the same structure

  if (length(blocks) > 0) {
    for (idx in 1:length(blocks)) {
      block <- blocks[[idx]]

      if (setequal(block@columns, columns)) {
        equals <- idx
      } else {
        violates <- c(violates, intersect(block@columns, columns))
      }
    }
  }

  return <- list("equals" = equals, "violates" = violates)
}

# checks if a vector has a continous sequence of 1s
.get_row_structure <- function(vector) {
  columns <- vector()

  if (length(vector) > 0) {
    for (idx in 1:length(vector)) {
      value <- vector[idx]
      if (value == 0) {
          columns <- c(columns, idx)
      }
    }
  }

  return(columns)
}

.explained_variance <- function(eigen_values, eigen_vectors, discard_cols, type = "approx") {
  switch(type,
         "approx" = { return(.explained_variance.approx(eigen_values, discard_cols)) },
         "exact" = { return(.explained_variance.exact(eigen_values, eigen_vectors, discard_cols)) }, {
           stop(paste("'", type, "'", " is not a valid value for explained variance.", sep = ""))
         })
}

.explained_variance.exact <- function(eigen_values, eigen_vectors, discard_cols) {
  sum <- 0

  if (length(eigen_vectors) > 0) {
    for (col in 1:ncol(eigen_vectors)) {
      sum_vec <- 0
      for (row in discard_cols) {
        sum_vec <- sum_vec + eigen_vectors[row, col]
      }
      sum <- sum + (eigen_values[col] * sum_vec ^ 2)
    }
  }

  return(sum / sum(eigen_values))
}

.explained_variance.approx <- function(eigen_values, discard_cols) {
  sum <- 0
  for (col in discard_cols) {
    sum <- sum + eigen_values[col]
  }

  return(sum / sum(eigen_values))
}
