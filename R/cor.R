select_cor <- function(x, cor) {
  if (cor == TRUE) {
    x <- cor(x)
  } else if (cor == FALSE) {
    x <- cov(x)
  } else {
    wrong_cor(cor=cor)
  }

  return(x) 
}

wrong_cor <- function(cor) {
  stop(
    paste(
      "'",
      cor,
      "'",
      " is not a valid value for cor.",
      sep=""
    )
  )
}