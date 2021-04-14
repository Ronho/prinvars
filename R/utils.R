.manipulate_matrix <- function(x, manipulator = "cov") {
  if (manipulator == "cor") return(cor(x, method = c("pearson")))
  if (manipulator == 'cov') return(cov(x, method = c("pearson")))
  stop(paste("'", manipulator, "'", " is not a valid value for manipulator.", sep = ""))
}


.create_blocks <- function(x) {
  blocks = list()

  for (idx in 1:nrow(x)) {
    obj <- .get_structure(x[, idx])

    if (obj$isValid == TRUE) {
      block_found <- FALSE
      if (length(blocks) > 0) {
        for (block_index in 1:length(blocks)) {
          block <- blocks[[block_index]]
          if (block@start == obj$start && block@end == obj$end) {
            blocks[[block_index]]@variables <- c(block@variables, idx)
            block_found <- TRUE
          }
        }
      }

      if (!block_found) {
        new_block <- new("Block", start = obj$start, end = obj$end, variables = c(idx))
        blocks <- c(blocks, new_block)
      }
    }
  }

  return(blocks)
}

# checks if a vector has a continous sequence of 1s
.get_structure = function(vector) {
  change <- 0 # number of changes from 0 to 1
  prev_val <- vector[1] # previous value
  start <- -1 # start of the 1s
  end <- -1 # end of the 1s

  for (i in seq_along(vector)) {
    value <- vector[i]

    if (value == 1 && start == -1) start <- i
    else if (value == 0 && prev_val == 1 && end == -1) end <- i - 1
    else if (i == length(vector) && end == -1) end <- i

    if (value != prev_val) {
      change <- change + 1
      prev_val <- value
    }
  }

  data <- list(
    "start" = start,
    "end" = end,
    "isValid" = (change <= 2 && change > 0) # Invalid if 0s are between 1s or only 0s
  )

  return(data)
}

.explained_variance = function(eigen_values, eigen_vectors, discard_cols, type = "approx") {
  switch(type,
         "approx" = { return(.explained_variance.approx(eigen_values, discard_cols)) },
         "exact" = { return(.explained_variance.exact(eigen_values, eigen_vectors, discard_cols)) }, {
           stop(paste("'", type, "'", " is not a valid value for explained variance.", sep = ""))
         })
}

.explained_variance.exact = function(eigen_values, eigen_vectors, discard_cols) {
  sum <- 0

  for (col in 1:ncol(eigen_vectors)) {
    sum_vec <- 0
    for (row in discard_cols) {
      sum_vec <- sum_vec + eigen_vectors[row, col]
    }
    sum <- sum + (eigen_values[col] * sum_vec ^ 2)
  }

  return(sum / sum(eigen_values))
}

.explained_variance.approx = function(eigen_values, discard_cols) {
  sum <- 0
  for (col in discard_cols) {
    sum <- sum + eigen_values[col]
  }

  return(sum / sum(eigen_values))
}
