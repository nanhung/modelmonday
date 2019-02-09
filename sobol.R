# This nblog post is inspired by 
# https://github.com/metrumresearchgroup/pbpk-qsp-mrgsolve/blob/master/docs/global_sensitivity_analysis.md
# https://github.com/mrgsolve/gallery
# I want to use some knowledge that I learn from 
# https://ec.europa.eu/jrc/en/event/training-course/samo-2018
# log scale; qusi Monte Carlo; Always plot the data 

library(tidyverse) 
library(mrgsolve) # mrgsim_ei
library(PKPDmisc) # auc_partial	
library(sensitivity) # sobol2007
library(randtoolbox) # sobol
library(reshape2) # melt
library(LSD) # heatscatter

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

par(mfrow = c(3,2), mar = c(2,2,3,1))
for (i in 1:5){
  heatscatter(log(x$X[,i]), log(x$y), xlab = "", ylab = "", main = names(x$X)[i])
}

par(mfrow = c(3,2), mar = c(2,2,3,1))
for (i in 1:5){
  heatscatter(log(x1$X[,i]), log(x1$y), xlab = "", ylab = "", main = names(x1$X)[i])
}

x$T[,"max. c.i."] - x$T[,"min. c.i."]
x1$T[,"max. c.i."] - x1$T[,"min. c.i."]
x2$T[,"max. c.i."] - x2$T[,"min. c.i."]

sample_converge <- function(n, l, which = names(l)){
  vars <- select_vars(names(l), !!(enquo(which)))
  m <- matrix(NA, length(n), length(vars))
  colnames(m) <- vars
  rownames(m) <- n
  m2 <- m1 <- m
  for (i in seq(length(n))){
    samp <- gen_samples(n[i], l, names(vars))
    samp1 <- gen_samples_1(n[i], l, names(vars))
    samp2 <- gen_samples_2(n[i], l, names(vars))
    x <- sobol2007(batch_run, X1=samp$x1, X2=samp$x2, nboot=100)
    x1 <- sobol2007(batch_run, X1=samp1$x1, X2=samp1$x2, nboot=100)
    x2 <- sobol2007(batch_run, X1=samp2$x1, X2=samp2$x2, nboot=100)
    m[i,] <- x$T[,"max. c.i."] - x$T[,"min. c.i."]
    m1[i,] <- x1$T[,"max. c.i."] - x1$T[,"min. c.i."]
    m2[i,] <- x2$T[,"max. c.i."] - x2$T[,"min. c.i."]
  } 
  X <- list(MC = m, log_MC = m1, log_QMC = m2)
  m %>% melt()
  
  return(X)
}

sample <- c(500, 1000, 2000, 4000, 8000)

set.seed(88771)
system.time(converge_list <- sample_converge(sample, param(mod), TVCL:TVVP))


df <- do.call(rbind, list(converge_list[[1]] %>% melt() %>% cbind(type = "MC"),
                          converge_list[[2]] %>% melt() %>% cbind(type = "log_MC"),
                          converge_list[[3]] %>% melt() %>% cbind(type = "log_QMC")))

theme_set(theme_light())


df %>% `colnames<-`(c("sample.no", "parameter", "index", "type")) %>%
  ggplot(aes(sample.no, index, group = parameter)) + geom_line(aes(color = parameter)) + 
  facet_wrap(~type) + 
  expand_limits(y= c(0, 0.5)) + geom_hline(yintercept = 0.05, linetype="dashed", size = 0.2) +
  labs(y = "Convergence index", x = "Sample number")

devtools::session_info()
