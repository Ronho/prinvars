select_eigen_vector_scaling <- function(eigen_vectors, scale=TRUE) {
  if (scale) {
    eigen_vectors <- scale_eigen_vectors(eigen_vectors=eigen_vectors)
  }

  return(eigen_vectors)
}

# scale eigen vectors (column wise) between -1 and 1
scale_eigen_vectors <- function(eigen_vectors) {
  scaled_eigen_vectors <- apply(
    eigen_vectors,
    MARGIN=2,
    FUN=function(x) {
      x / max(abs(x))
    }
  )

  return(scaled_eigen_vectors)
}