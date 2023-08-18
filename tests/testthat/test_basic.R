library(countsplit)

context("Poisson, 2, equal")

test_that("Poisson, 2, equal", {
  n <- 5000
  p <- 200
  X <- matrix(rpois(n*p, lambda=3), nrow=n)
  time <- system.time(split <- countsplit(X))[3]
  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]
  cors <- sapply(1:ncol(X), function(u) cor(Xtrain[,u], Xtest[,u]))
  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(sum(Xtrain)+sum(Xtest)==sum(X))
  expect_true(identical(Xtrain+Xtest, as(X, "sparseMatrix")))

  expect_true(max(abs(cors)) < 0.07)
  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.05)
  }
)

context("Poisson, 2, unequal")

test_that("Poisson, 2, unequal", {
  n <- 5000
  p <- 200
  X <- matrix(rpois(n*p, lambda=3), nrow=n)

  time <- system.time(split <- countsplit(X, epsilon=c(0.2,0.8)))[3]
  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(sum(Xtrain)+sum(Xtest)==sum(X))
  expect_true(identical(Xtrain+Xtest, as(X, "sparseMatrix")))

  cors <- sapply(1:ncol(X), function(u) cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.07)

  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)
  expect_true(abs(mean(Xtrain/0.2)-mean(X)) < 0.05)
}
)

context("Poisson, > 2, equal")

test_that("a", {
  n <- 5000
  p <- 200
  X <- matrix(rpois(n*p, lambda=3), nrow=n)


  ### Fast again
  time <- system.time(split <- countsplit(X, folds=5))[3]
  cat(time)
  Xtrain <- split[[1]]
  Xtest <-  as(X, "sparseMatrix") - Xtrain

  expect_true(identical(Xtest, split[[2]]+split[[3]]+split[[4]]+split[[5]]))

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(sum(Xtrain)+sum(Xtest)==sum(X))
  expect_true(identical(Xtrain+Xtest, as(X, "sparseMatrix")))

  cors <- sapply(1:ncol(X), function(u) cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.07)

  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)
  expect_true(abs(mean(Xtrain/0.2)-mean(X)) < 0.05)

  cat(time)
}
)

context("Poisson, >2, unequal")


test_that("a", {
  n <- 5000
  p <- 200
  X <- matrix(rpois(n*p, lambda=3), nrow=n)


  ### Incredibly slow lol, but this is not a very important case.
  time <- system.time(split <- countsplit(X, folds=6, epsilon=c(0.1,0.1,0.1,0.1,0.1,0.5)))[3]
  cat(time)

  expect_true(length(split)==6)
  expect_true(abs(mean(split[[1]]) - mean(split[[3]])) < 0.05)
  expect_true(abs(mean(split[[4]]) - mean(split[[2]])) < 0.05)

  Xtrain <- split[[1]]
  Xtest <-  as(X, "sparseMatrix") - Xtrain

  expect_true(identical(Xtest, split[[2]]+split[[3]]+split[[4]]+split[[5]] + +split[[6]]))

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(sum(Xtrain)+sum(Xtest)==sum(X))
  expect_true(identical(Xtrain+Xtest, as(X, "sparseMatrix")))

  cors <- sapply(1:ncol(X), function(u) cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.07)
  expect_true(abs(mean(Xtrain/0.1)-mean(Xtest/0.9)) < 0.05)
  expect_true(abs(mean(Xtrain/0.1)-mean(X)) < 0.05)


}
)

context("Overdisp=0, 2, equal")

test_that("overdisp=0, 2, equal", {
  n <- 5000
  p <- 200
  X <- matrix(rpois(n*p, lambda=3), nrow=n)

  #### Implements the fast version but we must make sure it works.
  time <- system.time( split <- countsplit(X, overdisps = rep(0,p)))[3]
  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(identical(Xtrain@x[Xtrain@x !=0], as(X, "sparseMatrix")@x[Xtrain@x !=0]))
  expect_true(identical(Xtest@x[Xtest@x !=0], as(X, "sparseMatrix")@x[Xtest@x !=0]))

  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.05)

  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(mean(cors < 0) > 0.9)
}
)

context("Overdisp=0, 2, unequal")

test_that("overdisp 0 2 unequal folds", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=5, size=3), nrow=n)

  ### Implements the second slowest version
  time <- system.time(split <- countsplit(X, overdisps = rep(0,p), folds=2, epsilon =c(0.2,0.8)))[3]
  expect_that(time, is_less_than(2),  info = cat(time))


  #cat(time)
  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(identical(Xtrain@x[Xtrain@x !=0], as(X, "sparseMatrix")@x[Xtrain@x !=0]))
  expect_true(identical(Xtest@x[Xtest@x !=0], as(X, "sparseMatrix")@x[Xtest@x !=0]))


  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)

  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(mean(cors < 0) > 0.9)
})

context("Overdisp=0, >2, equal")

test_that("a", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=5, size=3), nrow=n)

  ### Implements the slowest version
  time <- system.time(split <- countsplit(X, overdisps = rep(0,p), folds=5))[3]
  cat(time)

  Xtrain <- split[[4]]
  Xtest <- split[[2]]+split[[3]]+split[[5]]+split[[1]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(identical(Xtrain@x[Xtrain@x !=0], as(X, "sparseMatrix")@x[Xtrain@x !=0]))
  expect_true(identical(Xtest@x[Xtest@x !=0], as(X, "sparseMatrix")@x[Xtest@x !=0]))


  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)

  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(mean(cors < 0) > 0.95)
})

context("Overdisp=0, >2, unequal")

test_that("a", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=5, size=3), nrow=n)

  ### Implements the slowest version
  time <- system.time(split <- countsplit(X, overdisps = rep(0,p), folds=3, epsilon =c(0.2,0.3,0.5)))[3]
  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]+split[[3]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(identical(Xtrain@x[Xtrain@x !=0], as(X, "sparseMatrix")@x[Xtrain@x !=0]))
  expect_true(identical(Xtest@x[Xtest@x !=0], as(X, "sparseMatrix")@x[Xtest@x !=0]))


  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)

  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(mean(cors < 0) > 0.9)
})





context("NBright, 2, equal")

test_that("NBright, 2, equal", {
  n <- 5000
  p <- 200

  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)
  time <- system.time( split <- countsplit(X, overdisps = rep(1,p)))[3]

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.1)

  cat(time)
})

context("NBwrong, 2, equal")

test_that("NBwrong, 2, equal", {
  n <- 5000
  p <- 200

  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)

  time <- system.time(split <- countsplit(X, overdisps = exp(seq(-5,5,length.out=p))))[3]


  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(cor(1:p, cors) > 0.8)
}
)


context("NBright, 2, unequal")

test_that("NBright, 2, unequal", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)

  #### Implements the fast version but we must make sure it works.
  time <- system.time( split <- countsplit(X, overdisps = rep(1,p), epsilon=c(0.2,0.8)))[3]

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.1)
  cat(time)
})

context("NBwrong, 2, unequal")
test_that("NBwrong, 2, unequal", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)

  time <- system.time(split <- countsplit(X, overdisps = exp(seq(-5,5,length.out=p)), epsilon=c(0.2,0.8)))[3]

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain/0.2)-mean(Xtest/0.8)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(cor(1:p, cors) > 0.8)
  cat(time)
}
)


context("NBright, >2, equal")

test_that("NBright, >2, equal", {
  n <- 5000
  p <- 200

  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)

  time <- system.time( split <- countsplit(X, overdisps = rep(1,p), folds=4))[3]

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x+split[[3]]@x+split[[4]]@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.02)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(max(abs(cors)) < 0.1)

  cat(time)
})

context("NBwrong, >2, equal")

test_that("NBwrong, >2, equal", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)


  time <- system.time(split <- countsplit(X, overdisps = exp(seq(-5,5,length.out=p)), folds=4))[3]

  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x+split[[3]]@x+split[[4]]@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain)-mean(Xtest)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(cor(1:p, cors) > 0.8)
}
)

context("NBwrong, >2, unequal")

test_that("NBwrong, >2, unequal", {
  n <- 5000
  p <- 200
  X <- matrix(rnbinom(n*p, mu=3, size=1), nrow=n)


  time <- system.time(split <- countsplit(X, overdisps = exp(seq(-5,5,length.out=p)), folds=4, eps=c(0.1,0.4,0.1,0.4)))[3]

  cat(time)

  Xtrain <- split[[1]]
  Xtest <- split[[2]]+split[[3]]

  expect_true(all.equal(dim(Xtrain), dim(Xtest)))
  expect_true(all.equal(rownames(X), rownames(Xtrain)))
  expect_true(identical(Xtrain@x+Xtest@x+split[[4]]@x, as(X, "sparseMatrix")@x))

  expect_true(abs(mean(Xtrain/0.1)-mean(Xtest/0.5)) < 0.05)
  cors <- sapply(1:p, function(u)  cor(Xtrain[,u], Xtest[,u]))
  expect_true(cor(1:p, cors) > 0.8)
}
)

