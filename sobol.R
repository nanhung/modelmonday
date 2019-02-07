library(mrgsolve)
library(tidyverse)
library(PKPDmisc)
library(sensitivity)
library(randtoolbox)

mod <- mread("sunit", "model") %>% 
  update(end = 24, delta = 1) %>% zero_re

see(mod)

sunev <- function(amt = 50,...) ev(amt = amt, ...)

gen_samples <- function(n, l, which = names(l), factor = c(0.01,100)) {
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

gen_samples <- function(n, l, which = names(l), factor = c(0.01,100)) {
  vars <- select_vars(names(l), !!(enquo(which)))
  l <- as.list(l)[vars]
  l <- lapply(l, function(x) {
    x*factor  
  })
  n <- length(l)*n*2
  df <- as.data.frame(l)
  len <- length(df)
  #X <- matrix(ncol=len, nrow=n)
  #colnames(X) <- names(df)
  #Y <- X
  
  X <- sobol(n = 1000, dim = 5, seed = 1234, scrambling = 3)
  Y <- sobol(n = 1000, dim = 5, seed = 2345, scrambling = 3)
  
  for(i in seq(len)){
    r <- runif(n, df[1,i], df[2,i])
    X[,i] <- r
    r <- runif(n, df[1,i], df[2,i])
    Y[,i] <- r
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

samp <- gen_samples(1000, param(mod), TVCL:TVVP)
head(samp$x1)

x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=100)

x

plot(x)

devtools::session_info()
