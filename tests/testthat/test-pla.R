obj <- pla(
  x=USArrests,
  cor=FALSE,
)

test_that("pla", {  
  blocks <- list()
  blocks[[1]] <- new("Block", features=c("Murder"), explained_variance=0.966)
  blocks[[2]] <- new("Block", features=c("Assault"), explained_variance=0.028)
  blocks[[3]] <- new("Block", features=c("UrbanPop"), explained_variance=0.006)#58)
  blocks[[4]] <- new("Block", features=c("Rape"), explained_variance=0.0008)
  expect_equal(obj$blocks, blocks, tolerance=1e-2)
})

test_that("pla.keep_blocks", {
  result <- pla.keep_blocks(object=obj, blocks=c(1,3))

  expect_equal(result$x, USArrests[c(1,3)])
})

test_that("pla.drop_blocks", {
  result <- pla.drop_blocks(object=obj, blocks=c(1,3))

  expect_equal(result$x, USArrests[c(2,4)])
})