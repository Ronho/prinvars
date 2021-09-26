matrix <- matrix(c(1:9), nrow=3, ncol=3)

test_that("select_cor", {
  expect_equal(select_cor(x=matrix, cor=TRUE), cor(matrix))
  expect_equal(select_cor(x=matrix, cor=FALSE), cov(matrix))
  expect_error(select_cor(x=matrix, cor=""))  
})

test_that("wrong_cor", {
  expect_error(wrong_cor(cor="Test"))
})
