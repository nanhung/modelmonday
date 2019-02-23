# https://statmodeling.stat.columbia.edu/2014/03/10/stan-model-pk-iv-oral-dosing/
library(rstan)

nIV <- 12
nOral <- 11
doseIV <- 500
doseOral <- 500
timeIV <- c(0.33, 0.5, 0.67, 1.5, 2, 4, 6, 10, 16, 24, 32, 48)
concIV <- c(14.7, 12.6, 11, 9, 8.2, 7.9, 6.6, 6.2, 4.6, 3.2, 2.3, 1.2)
timeOral <- c(0.5, 1, 1.5, 2, 4, 6, 10, 16, 24, 32, 48)
concOral <- c(2.4, 3.8, 4.2, 4.6, 8.1, 5.8, 5.1, 4.1, 3, 2.3, 1.3)


fit <- stan('models/pk_iv_oral.stan', 
            data=c("nIV","nOral","doseIV","doseOral", 
                   "timeIV","concIV","timeOral","concOral"),
            chains=4, warmup=5000, iter=20000, control = list(adapt_delta = 0.85))

fit

pairs(fit)
