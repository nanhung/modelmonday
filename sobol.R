library(mrgsolve)
library(tidyverse)
library(PKPDmisc)
library(sensitivity)
library(randtoolbox)

mod <- mread("sunit", "model") %>% 
  update(end = 24, delta = 1) %>% zero_re

see(mod)

sunev <- function(amt = 50,...) ev(amt = amt, ...)

gen_samples_1 <- function(n, l, which = names(l), factor = c(0.01,100)) {
  vars <- select_vars(names(l), !!(enquo(which)))
  l <- as.list(l)[vars]
  l <- lapply(l, function(x) {
    x*factor  
  })
  n <- length(l)*n*2
  df <- as.data.frame(l)
  len <- length(df)
  X <- matrix(ncol=len, nrow=n)
  colnames(X) <- names(df)
  Y <- X
  for(i in seq(len)){
    r <- runif(n, df[1,i], df[2,i])
    X[,i] <- r
    r <- runif(n, df[1,i], df[2,i])
    Y[,i] <- r
  }
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}

xx=4
X <- sobol(100, dim = 5)
plot(X[,1], X[,2])
X[,1] <- 10^(X[,1] * xx) 
plot(X[,1])

gen_samples_2 <- function(n, l, which = names(l), factor = c(0.01,100)) {
  vars <- select_vars(names(l), !!(enquo(which)))
  l <- as.list(l)[vars]
  l <- lapply(l, function(x) {x*factor})
  xx <- log(factor, 10)[2] - log(factor, 10)[1]
  len <- length(vars)

  X <- sobol(n = length(l)*n*2, dim = 5)
  Y <- sobol(n = length(l)*n*2, dim = 5, seed = 2345, scrambling = 3)
  
  for(i in seq(len)){
    X[,i] <- l[[i]][[1]] * 10^(X[,i] * xx)
    Y[,i] <- l[[i]][[1]] * 10^(Y[,i] * xx)
    colnames(X) <- colnames(Y) <- vars 
  }
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}

sim_chunk <- function(mod, x) {
  mrgsim_ei(x = mod, ev = sunev(), idata = x, obsonly = TRUE) %>% 
    as_data_frame
}

batch_run <- function(x) {
  out <- sim_chunk(mod,x)
  out <- 
    group_by(out,ID) %>% 
    summarise(AUC = auc_partial(time,CP))
  return(out$AUC)
}

samp1 <- gen_samples_1(1000, param(mod), TVCL:TVVP)
samp2 <- gen_samples_2(1000, param(mod), TVCL:TVVP)
head(samp1$x1)
head(samp2$x1)

i=2
plot(density(samp1$x1[,i]))
plot(density(samp1$x2[,i]))

heatscatter(samp1$x1[,1], samp1$x1[,2], add.contour=T, nlevels=3)
heatscatter(samp2$x1[,1], samp2$x1[,2], add.contour=T, nlevels=3)

heatscatter(log(samp1$x1[,1]), log(samp1$x1[,2]), add.contour=T, nlevels=3)
heatscatter(log(samp2$x1[,1]), log(samp2$x1[,2]), add.contour=T, nlevels=3)

x1 <- sobol2007(batch_run, X1=samp1$x1, X2=samp1$x2, nboot=100)
x2 <- sobol2007(batch_run, X1=samp2$x1, X2=samp2$x2, nboot=100)

x1
x2





par(mfrow = c(2,1))
plot(x1)
plot(x2)

devtools::session_info()
