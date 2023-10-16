find_blocks <- function(matrix, column_idx) {
    if (length(column_idx) == 1) {
      rows <- which(matrix[, column_idx] != 0)
    } else {
      rows <- which(rowSums(matrix[, column_idx] != 0) > 0)
    }
    
    if (length(rows) == 1) {
      cols <- which(matrix[rows, ] != 0)
    } else {
      cols <- which(colSums(matrix[rows, ] != 0) > 0)
    }
    
    return(cols)
  }
  
  
get_blocks <- function(threshold_matrix, feature_names){
  
    P <- nrow(threshold_matrix)
    Q <- ncol(threshold_matrix)
    columns <- 1:Q
    blocks <- list()
    
    if (any(colSums(threshold_matrix == 0) == P)) {
      stop("Zero column. Reduce penalization/threshold value.\n")
    } 
    
    if (any(rowSums(threshold_matrix == 0) == Q)) {
      if(Q < P){
        stop("Zero row. Add more loadings, or reduce penalization (para) and/or threshold value (tau).\n")
      }
      else{
        stop("Zero row. Reduce penalization (para) and/or threshold value (tau).\n")
      }
      
    } 
    
    
    
    while(length(columns) != 0){
      block_columns <- columns[1]
      while(!identical(block_columns, find_blocks(threshold_matrix, block_columns))){
        block_columns <- find_blocks(threshold_matrix, block_columns)
      }
      
      if (length(block_columns) == 1) {
        block_idx <- which(threshold_matrix[, block_columns] != 0)
      } else {
        block_idx <- which(rowSums(threshold_matrix[, block_columns] != 0) > 0)
      }
      
      blocks[[length(blocks) + 1]] <- create_block(
        feature_names=feature_names,
        selected_features=block_idx,
        is_valid=TRUE,
        ev_influenced=block_columns
      )
      
      columns <- columns[!columns %in% block_columns]
    }
    
    return(blocks)
    
  } 
