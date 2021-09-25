matrix <- matrix(c(1:9), nrow=3, ncol=3)

test_that("select_cov", {
  expect_equal(select_cov(x=matrix, cov=TRUE), cov(matrix))
  expect_equal(select_cov(x=matrix, cov=FALSE), cor(matrix))
  expect_error(select_cov(x=matrix, cov=""))  
})

test_that("wrong_cov", {
  expect_error(wrong_cov(cov="Test"))
})
