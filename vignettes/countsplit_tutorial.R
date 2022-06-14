## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("anna-neufeld/countsplit")

## ---- message=FALSE-----------------------------------------------------------
library(countsplit)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
p <- 200
X <- matrix(rpois(n*p, lambda=5), nrow=n)

## -----------------------------------------------------------------------------
clusters.full <- kmeans(log(X+1), centers=2)$cluster
results.naive <- t(apply(X, 2, function(u) summary(glm(u~clusters.full, family="poisson"))$coefficients[2,]))
head(results.naive)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(data=NULL, aes(sample=results.naive[,4]))+geom_qq(distribution="qunif")+geom_abline(col="red")

## -----------------------------------------------------------------------------
set.seed(2)
split <- countsplit(X, epsilon=0.5)
names(split)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
results.countsplit <- t(apply(Xtest, 2, function(u) summary(glm(u~clusters.train, family="poisson"))$coefficients[2,]))
head(results.countsplit)

## -----------------------------------------------------------------------------
ggplot(data=NULL, aes(sample=results.countsplit[,4]))+geom_qq(distribution="qunif")+geom_abline(col="red")

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
p <- 200
clusters.true <- rbinom(n, size=1, prob=0.5)
Lambda <- matrix(5, nrow=n, ncol=p)
Lambda[clusters.true==1, 1:10] <- 10
X <-apply(Lambda,1:2,rpois,n=1)

## -----------------------------------------------------------------------------
set.seed(111)
split <- countsplit(X, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
set.seed(222)
clusters.full <- kmeans(log(X+1), centers=2)$cluster
table(clusters.true, clusters.full)

## -----------------------------------------------------------------------------
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
table(clusters.train, clusters.full)

## -----------------------------------------------------------------------------
intercepts.ideal <- t(apply(X, 2, function(u) summary(glm(u~clusters.true, family="poisson"))$coefficients[1,1]))
slopes.ideal <- t(apply(X, 2, function(u) summary(glm(u~clusters.true, family="poisson"))$coefficients[2,1]))

## -----------------------------------------------------------------------------
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
intercepts.countsplit <- t(apply(Xtest, 2, function(u) summary(glm(u~clusters.train, family="poisson"))$coefficients[1,1]))
slopes.countsplit <- t(apply(Xtest, 2, function(u) summary(glm(u~clusters.train, family="poisson"))$coefficients[2,1]))

## -----------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
p1 <- ggplot(data=NULL, aes(x=intercepts.ideal, y=intercepts.countsplit))+geom_point()+geom_abline(intercept= log(0.5), slope=1, col="red")
p2 <- ggplot(data=NULL, aes(x=slopes.ideal, y=slopes.countsplit))+geom_point()+geom_abline(intercept=0, slope=1, col="red")
p1 + p2

