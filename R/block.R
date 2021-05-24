#' Block
#'
#' Classed used within the package to keep the structure and information about
#' the generated blocks. A Block is a a continuous sequence of 1s
#' within a vector.
#'
#' @slot rows vector of numeric, indizes of the rows which belong to the block
#' @slot columns vector of numeric, indizes of the columns which belong to the block
#' @slot explained_varaiance numeric, explained variance of the blocks variables
#' based on the whole data set.
setClass(
  "Block",
  representation(
    rows = "vector",
    columns = "vector",
    explained_variance = "numeric"
  )
)

#' @title Block - Show
#'
#' Prints the blocks structure.
#'
#' @param object block.
#'
#' @examples
#' block <- new("Block", rows = c(1, 4), columns = c(2, 5),
#' explained_variance = 0.03)
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
#' block <- new("Block", rows = c(1, 4), columns = c(2, 5),
#' explained_variance = 0.03)
#' str(block)
#' @export
setMethod(
  f = "str",
  signature = "Block",
  definition = function(object) {
    rows = paste(unlist(object@rows), collapse = ", ")
    columns = paste(unlist(object@columns), collapse = ", ")
    expvar = round(object@explained_variance * 100, 2)
    str = paste("Columns(",
                columns,
                ")-Rows(",
                rows,
                ") explains ",
                expvar,
                "% of the overall explained variance",
                sep = "")
    return(str)
  }
)
