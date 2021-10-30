#' Block
#'
#' Class used within the package to keep the structure and information about
#' the generated blocks.
#'
#' @slot features a vector of numeric which contains the indices of the block.
#' @slot explained_varaiance a numeric which contains the variance explained of
#' the blocks variables based on the whole data set.
#' @slot is_valid a logival which indicates if the block structure is valid.
setClass(
  "Block",
  representation(
    features = "vector",
    explained_variance = "numeric",
    is_valid = "logical"
  ),
  prototype(explained_variance = 0, is_valid=TRUE)
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
    if (object@is_valid) {
      str = paste(
        "Features (",
        features,
        ") explain ",
        expvar,
        "% of the overall explained variance",
        sep = ""
      )
    } else {
      str = paste(
        "Features (",
        features,
        ") remain without a block structure row-wise.",
        sep = ""
      )
    }

    return(str)
  }
)
