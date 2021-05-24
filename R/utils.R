.manipulate_matrix <- function(x, manipulator = "cov") {
  if (manipulator == "cor") return(cor(x, method = c("pearson")))
  if (manipulator == 'cov') return(cov(x, method = c("pearson")))
  stop(paste("'", manipulator, "'", " is not a valid value for manipulator.", sep = ""))
}


.create_blocks <- function(x) {
  blocks <- list()
  violated <- c()
  for (idx in 1:nrow(x)) {
    rows <- .get_column_structure(x[, idx])
    block_idx <- .get_block_index(blocks, rows)
    violated <- unique(c(violated, block_idx$violates))

    if (block_idx$equals == -1) {
      block <- new("Block", rows = rows, columns = c(idx))
      blocks[[length(blocks)+1]] <- block
    } else {
      blocks[[block_idx$equals]]@columns <- c(blocks[[block_idx$equals]]@columns, idx)
    }
  }

  # remove blocks that violate the n x n structure
  uncontained_columns = c()
  if (length(blocks) > 0) {
    idx <- 1
    for (i in 1:length(blocks)) {
      block <- blocks[[idx]]
      if ((length(block@columns) != length(block@rows)) || (length(intersect(block@rows, violated)) > 0)) {
        uncontained_columns <- unique(c(uncontained_columns, block@columns))
        blocks[[idx]] <- NULL
      } else {
        idx <- idx+1
      }
    }
  }

  # add a new 1x1 block for each column that is not contained within any block
  for (column_nr in uncontained_columns) {
    blocks[[length(blocks)+1]] <- new("Block", rows = vector(), columns = c(column_nr))
  }

  return(blocks)
}

# returns the index of the block containing the exact columns
# returns -1 if no block was found
# returns the row indexes of blocks that can be confirmed to violate the block structure based on the new columns
# (e.g. one block contains rows 2 and 8 --> if this functions retrieves c(8,9), the required structure cannot be given since 9 is not part of c(2,8))
.get_block_index <- function(blocks, rows) {
  violates <- vector() # indexes of blocks that have a confirmed violated structure
  equals <- -1 # index of block that has the same structure

  if (length(blocks) > 0) {
    for (idx in 1:length(blocks)) {
      block <- blocks[[idx]]

      if (setequal(block@rows, rows)) {
        equals <- idx
      } else {
        violates <- c(violates, intersect(block@rows, rows))
      }
    }
  }

  return <- list("equals" = equals, "violates" = violates)
}

# checks where the rows are that contain 0s
.get_column_structure <- function(vector) {
  rows <- vector()

  if (length(vector) > 0) {
    for (idx in 1:length(vector)) {
      value <- vector[idx]
      if (value == 0) {
          rows <- c(rows, idx)
      }
    }
  }

  return(rows)
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
