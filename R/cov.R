select_cov <- function(x, cov) {
  if (cov == TRUE) {
    x <- cov(x)
  } else if (cov == FALSE) {
    x <- cor(x)
  } else {
    wrong_cov(cov=cov)
  }

  return(x) 
}

wrong_cov <- function(cov) {
  stop(
    paste(
      "'",
      cov,
      "'",
      " is not a valid value for cov.",
      sep=""
    )
  )
}