# https://github.com/stan-dev/stancon_talks
# https://mc-stan.org/events/stancon2017-notebooks/stancon2017-margossian-gillespie-ode.html
# rm(list=ls())

library(tidyverse)
library(rstan)
library(bayesplot)

modelName <- "twoCpt"
data <- read_rdump(file.path(paste0("data/", modelName,".data.R")))

init <- function(){
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       ke0 = exp(rnorm(1,log(1),0.2)),
       EC50 = exp(rnorm(1,log(100),0.2)),
       sigma = 0.5,
       sigmaResp = 20)
}

## Specify the variables for which you want history plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

nChains <- 4
nPost <- 1000 ## Number of post-warm-up samples per chain after thinning
nBurn <- 1000 ## Number of warm-up samples per chain after thinning
nThin <- 1
nIter <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(paste("models/", modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))

#save(fit, file = file.path(modelDir, paste(modelName, "Fit.Rsave", sep = "")))
#mcmcplots::mcmcDensity(fit, parametersToPlot, byChain = TRUE)

xdata <- data.frame(data$cObs, data$time[data$evid != 1])
xdata <- plyr::rename(xdata, c("data.cObs" = "cObs", "data.time.data.evid....1." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

plot(fit)
fit
