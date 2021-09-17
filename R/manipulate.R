select_manipulator <- function(x, manipulator) {
  x = switch(manipulator,
    "cor"=cor(x, method=c("pearson")),
    "cov"=cov(x, method=c("pearson")),
    "none"=validate_matrix_structure(x=x),
    wrong_manipulator(manipulator=manipulator)
  )

  return(x) 
}

validate_matrix_structure <- function(x) {
  if (!is_quadratic_matrix(x=x)) {
    warning("The given matrix is not quadratic.")
  }

  return(x)
}

is_quadratic_matrix <- function(x) {
  return(nrow(x) == ncol(x))
}

wrong_manipulator <- function(manipulator) {
  stop(
    paste(
      "'",
      manipulator,
      "'",
      " is not a valid value for manipulator.",
      sep=""
    )
  )
}