#' Block
#'
#' Classed used within the package to keep the structure and information about
#' the generated blocks. A Block is a a continuous sequence of 1s
#' within a vector.
#'
#' @slot start numeric, beginning index of the sequence; -1 indicates that no 1 is found.
#' @slot end numeric, ending index of the sequence; -1 indicates that no 1 is found.
#' @slot variables vector of numeric, indizes of the variables having the same
#' start and end.
#' @slot explained_varaiance numeric, explained variance of the blocks variables
#' based on the whole data set.
setClass(
  "Block",
  representation(
    start = "numeric",
    end = "numeric",
    variables = "vector",
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
#' block <- new("Block", start = 1, end = 4, variables = c(2, 5),
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
#' block <- new("Block", start = 1, end = 4, variables = c(2, 5),
#' explained_variance = 0.03)
#' str(block)
#' @export
setMethod(
  f = "str",
  signature = "Block",
  definition = function(object) {
    s1 = paste(unlist(object@variables), collapse = ", ")
    s2 = round(object@explained_variance * 100, 2)
    str = paste("(",
                s1,
                ")",
                " explains ",
                s2,
                "% of the overall explained variance",
                sep = "")
    return(str)
  }
)
