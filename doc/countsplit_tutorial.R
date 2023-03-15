## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("anna-neufeld/countsplit")

## ---- message=FALSE-----------------------------------------------------------
library(countsplit)
library(ggplot2)
library(patchwork)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1000
p <- 200
X <- matrix(rpois(n*p, lambda=5), nrow=n)

## -----------------------------------------------------------------------------
clusters.full <- kmeans(log(X+1), centers=2)$cluster
results.naive <- t(apply(X, 2, function(u) summary(glm(u~as.factor(clusters.full), family="poisson"))$coefficients[2,]))
head(results.naive)

## -----------------------------------------------------------------------------
ggplot(data=NULL, aes(sample=results.naive[,4]))+geom_qq(distribution="qunif")+geom_abline(col="red")

## -----------------------------------------------------------------------------
set.seed(2)
split <- countsplit(X, epsilon=0.5)
names(split)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
results.countsplit <- t(apply(Xtest, 2, function(u) summary(glm(u~as.factor(clusters.train), family="poisson"))$coefficients[2,]))
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
split <- countsplit(X, epsilon=0.5)
Xtrain <- split$train
Xtest <- split$test

## -----------------------------------------------------------------------------
clusters.full <- kmeans(log(X+1), centers=2)$cluster
table(clusters.true, clusters.full)

## -----------------------------------------------------------------------------
clusters.train <- kmeans(log(Xtrain+1), centers=2)$cluster
table(clusters.true, clusters.train)

## -----------------------------------------------------------------------------
coeffs.X <- t(apply(X, 2, function(u) summary(glm(u~as.factor(clusters.true), family="poisson"))$coefficients[,1]))
coeffs.Xtest <- t(apply(Xtest, 2, function(u) summary(glm(u~as.factor(clusters.true), family="poisson"))$coefficients[,1]))

## -----------------------------------------------------------------------------
differentially_expressed = as.factor(c(rep(1,10), rep(0,190)))
p1 <- ggplot(data=NULL, aes(x=coeffs.X[,1], y=coeffs.Xtest[,1], col=differentially_expressed))+
  geom_point()+
  geom_abline(intercept= log(0.5), slope=1, col="red")+
  geom_abline(intercept= 0, slope=1, col="red", lty=2)+
  coord_fixed()+xlim(0,2)+ylim(0,2)+
  xlab("Intercepts from X")+ ylab("Intercepts from Xtest")+
  ggtitle("Intercepts")+theme_bw()
p2 <- ggplot(data=NULL, aes(x=coeffs.X[,2], y=coeffs.Xtest[,2], col=differentially_expressed))+
  geom_point()+
  geom_abline(intercept=0, slope=1, col="red")+
  coord_fixed()+xlim(-0.15,1)+ylim(-0.15,1)+
  xlab("Slopes from X")+ ylab("Slopes from Xtest")+
  ggtitle("Slopes")+theme_bw()
p1 + p2 + plot_layout(guides="collect")

## -----------------------------------------------------------------------------
coeffs.ideal <- t(apply(X, 2, function(u) summary(glm(u~as.factor(clusters.true), family="poisson"))$coefficients[,1]))
coeffs.countsplit <- t(apply(Xtest, 2, function(u) summary(glm(u~as.factor(clusters.train==1), family="poisson"))$coefficients[,1]))
p1 <- ggplot(data=NULL, aes(x=coeffs.ideal[,1], y=coeffs.countsplit[,1], col=differentially_expressed))+
  geom_point()+
  geom_abline(intercept= log(0.5), slope=1, col="red")+
  geom_abline(intercept= 0, slope=1, col="red", lty=2)+
  coord_fixed()+xlim(0,2)+ylim(0,2)+
  xlab("Intercepts from ideal method")+ ylab("Intercepts from count splitting")+
  ggtitle("Intercepts")
p2 <- ggplot(data=NULL, aes(x=coeffs.ideal[,2], y=coeffs.countsplit[,2], col=differentially_expressed))+
  geom_point()+
  geom_abline(intercept=0, slope=1, col="red")+
  coord_fixed()+xlim(-0.15,1)+ylim(-0.15,1)+
  xlab("Slopes from ideal method")+ ylab("Slopes from count splitting")+
  ggtitle("Slopes")
p1 + p2 + plot_layout(guides="collect") & theme_bw()

