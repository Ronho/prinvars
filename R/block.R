#' Block
#'
#' Classed used within the package to keep the structure and information about
#' the generated blocks.
#'
#' @slot features vector of numeric, indices which belong to the block
#' @slot explained_varaiance numeric, explained variance of the blocks variables
#' based on the whole data set.
setClass(
  "Block",
  representation(
    features = "vector",
    explained_variance = "numeric"
  ),
  prototype(explained_variance = 0)
)

#' @title Block - Show
#'
#' Prints the blocks structure.
#'
#' @param object block.
#'
#' @examples
#' block <- new("Block", features = c(2, 5), explained_variance = 0.03)
#' print(block)
#' @export
setMethod(
  f = "show",
  signature = "Block",
  definition = function(object) {
    print(str(object))
  }
)

#' @title Block - str
#'
#' Generic function to create a string out of the blocks structure.
#'
#' @param object block.
#'
#' @examples
#' block <- new("Block", features = c(2, 5), explained_variance = 0.03)
#' str(block)
#' @export
setMethod(
  f = "str",
  signature = "Block",
  definition = function(object) {
    features = paste(unlist(object@features), collapse = ", ")
    expvar = round(object@explained_variance * 100, 2)
    str = paste("Features (",
                features,
                ") explain ",
                expvar,
                "% of the overall explained variance",
                sep = "")
    return(str)
  }
)
