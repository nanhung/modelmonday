# This nblog post is inspired by 
# https://github.com/metrumresearchgroup/pbpk-qsp-mrgsolve/blob/master/docs/global_sensitivity_analysis.md
# I want to use some knowledge that I learn from 
# https://ec.europa.eu/jrc/en/event/training-course/samo-2018
# log scale; qusi Monte Carlo; Always plot the data 

library(mrgsolve) # mrgsim_ei
library(tidyverse) # select_vars
library(PKPDmisc) # auc_partial	
library(sensitivity) # sobol2007
library(randtoolbox)
library(LSD)

mod <- mread("sunit", "models") %>% 
  update(end = 24, delta = 1) %>% zero_re

see(mod)

sunev <- function(amt = 50,...) ev(amt = amt, ...)

gen_samples<- function(n, l, which = names(l), factor = c(0.01,100)) {
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

set.seed(88771)
samp <- gen_samples(6000, param(mod), TVCL:TVVP) 
head(samp$x1)
dim(samp$x1)

system.time(x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=100))
x
plot(x)

####

x$T[,"max. c.i."] - x$T[,"min. c.i."]

gen_samples_1 <- function(n, l, which = names(l), factor = c(0.01,100)) {
  vars <- select_vars(names(l), !!(enquo(which)))
  l <- as.list(l)[vars]
  l <- lapply(l, function(x) {x*factor})
  xx <- log(factor, 10)[2] - log(factor, 10)[1]
  len <- length(vars)
  X <- matrix(runif(len * length(l)*n*2), ncol = len)
  Y <- matrix(runif(len * length(l)*n*2), ncol = len)
  for(i in seq(len)){
    X[,i] <- l[[i]][[1]] * 10^(X[,i] * xx)
    Y[,i] <- l[[i]][[1]] * 10^(Y[,i] * xx)
    colnames(X) <- colnames(Y) <- vars 
  }
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}

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

set.seed(88771)
samp <- gen_samples(1000, param(mod), TVCL:TVVP)
samp1 <- gen_samples_1(1000, param(mod), TVCL:TVVP)
samp2 <- gen_samples_2(1000, param(mod), TVCL:TVVP)

head(samp$x1)
head(samp1$x1)
head(samp2$x1)

i=1
range(samp$x1[,i])
range(samp1$x1[,i])
range(samp2$x1[,i])


# The devil is in the details
samp$x1[,i] %>% density() %>% plot()
samp1$x1[,i] %>% density() %>% plot()
samp1$x1[,i] %>% log(10) %>% density() %>% plot()
samp2$x1[,i] %>% log(10) %>% density() %>% plot()
samp$x1[,i] %>% log(10) %>% density() %>% plot()

heatscatter(samp$x1[,1], samp$x1[,2], add.contour=T, nlevels=3, xlab = "TVCL", ylab = "TVKA")
heatscatter(samp1$x1[,1], samp1$x1[,2], add.contour=T, nlevels=3, xlab = "TVCL", ylab = "TVKA")
heatscatter(log(samp1$x1[,1]), log(samp1$x1[,2]), add.contour=T, nlevels=3, xlab = "TVCL", ylab = "TVKA")
heatscatter(log(samp2$x1[,1]), log(samp2$x1[,2]), add.contour=T, nlevels=3, xlab = "TVCL", ylab = "TVKA")
heatscatter(log(samp$x1[,1]), log(samp$x1[,2]), add.contour=T, nlevels=3, xlab = "TVCL", ylab = "TVKA")

system.time(x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=100))
system.time(x1 <- sobol2007(batch_run, X1=samp1$x1, X2=samp1$x2, nboot=100))
system.time(x2 <- sobol2007(batch_run, X1=samp2$x1, X2=samp2$x2, nboot=100))

x
x1
x2

plot(x)
plot(x1)
plot(x2)

x$T[,"max. c.i."] - x$T[,"min. c.i."]
x1$T[,"max. c.i."] - x1$T[,"min. c.i."]
x2$T[,"max. c.i."] - x2$T[,"min. c.i."]

# Convergence



devtools::session_info()
