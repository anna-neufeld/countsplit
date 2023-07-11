context("basic testing")
library(countsplit)

test_that("Fake Data Poisson", {
  
  set.seed(1)
  n <- 1000
  p <- 200
  X <- matrix(rpois(n*p, lambda=5), nrow=n)
  set.seed(2)
  split <- countsplit(X)
  Xtrain <- split[[1]]
  Xtest <- split[[2]]
  clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
  
  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  }
)


context("basic testing")
library(countsplit)

test_that("Fake Data NB", {
  
  set.seed(1)
  n <- 1000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=5, size=3), nrow=n)
  set.seed(2)
  split <- countsplit(X, overdisps = rep(3,p))
  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 1)
}
)
