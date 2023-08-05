obj <- spla(
            USArrests,
            method = "spca",
            para = c(0.5, 0.5, 0.5, 0.5),
            cor = TRUE)

test_that("spla blocks", {
  blocks <- list()
  blocks[[1]] <- new(
                     "Block",
                     features=c("UrbanPop"),
                     explained_variance=0.25,
                     ev_influenced = c(1))
  blocks[[2]] <- new(
                     "Block",
                     features=c("Murder", "Assault", "Rape"),
                     explained_variance=0.6849,
                     ev_influenced = c(2, 3, 4))

  expect_equal(obj$blocks, blocks, tolerance=1e-2)
})

test_that("spla EC", {
  ec <- matrix(c(0, 0.92, 0, 0), nrow = 1, ncol = 4)
  dim(ec) <- NULL
  names(ec) <- c("[,2]", "[,1]", "[,3]", "[,4]")

  expect_equal(obj$EC, ec, tolerance=1e-2)
})

test_that("spla threshold", {
  expect_equal(obj$threshold, 0, tolerance=1e-2)
})

test_that("spla threshold_mode", {
  expect_equal(obj$threshold_mode, "cutoff", tolerance=1e-2)
})