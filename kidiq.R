## 3.1 one_predictor ####
library(rstan)
library(grid)
library(ggplot2)

### Data
source("ARM/kidiq.data.R", echo = TRUE)

### First model: kid_score ~ mom_hs 
data.list.1 <- c("N", "kid_score", "mom_hs")
kidscore_momhs <- stan(file = "ARM/kidscore_momhs.stan", data = data.list.1, iter = 500)
kidscore_momhs
pairs(kidscore_momhs)

### Second model: lm(kid_score ~ mom_iq) 
data.list.2 <- c("N", "kid_score", "mom_iq")
kidscore_momiq <- stan(file = 'ARM/kidscore_momiq.stan', data = data.list.2, iter = 500)
kidscore_momiq
pairs(kidscore_momiq)

# Figure 3.1
kidiq.data <- data.frame(kid_score, mom_hs, mom_iq)
beta.post.1 <- extract(kidscore_momhs, "beta")$beta
beta.mean.1 <- colMeans(beta.post.1)
p1 <- ggplot(kidiq.data, aes(x = mom_hs, y = kid_score)) +
  geom_jitter(position = position_jitter(width = 0.05, height = 0)) +
  geom_abline(aes(intercept = beta.mean.1[1], slope = beta.mean.1[2])) +
  scale_x_continuous("Mother completed high school", breaks = c(0, 1)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  theme_bw()
print(p1)

# Figure 3.2
beta.post.2 <- extract(kidscore_momiq, "beta")$beta
beta.mean.2 <- colMeans(beta.post.2)
p2 <- ggplot(kidiq.data, aes(x = mom_iq, y = kid_score)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean.2[1], slope = beta.mean.2[2])) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  theme_bw()
print(p2)

## 3.2_MultiplePredictors ####
### Model: kid_score ~ mom_hs + mom_iq 
data.list <- c("N", "kid_score", "mom_hs", "mom_iq")
kidiq_multi_preds <- stan(file = 'ARM/kidiq_multi_preds.stan',
                          data = data.list,
                          iter = 500, chains = 4)
kidiq_multi_preds
pairs(kidiq_multi_preds)

# Figure 3.3
beta.post <- extract(kidiq_multi_preds, "beta")$beta
beta.mean <- colMeans(beta.post)
kidiq.data <- data.frame(kid_score, mom_hs = as.factor(mom_hs), mom_iq)
levels(kidiq.data$mom_hs) <- c("No", "Yes")

p <- ggplot(kidiq.data, aes(x = mom_iq, y = kid_score, color = mom_hs)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean[1] + beta.mean[2] * (mom_hs == "Yes"),
                  slope = beta.mean[3], color = mom_hs)) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  scale_color_manual("Mother\ncompleted\nhigh\nschool",
                     values = c("No" = "black", "Yes" = "gray")) +
  theme_bw()
print(p)

## 3.3_Interactions ####
### Model: lm(kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq) 
data.list <- c("N", "kid_score", "mom_hs", "mom_iq")
kidiq_interaction <- stan(file = 'ARM/kidiq_interaction.stan',
                          data = data.list,
                          iter = 1000, chains = 4)
kidiq_interaction
pairs(kidiq_interaction)

### Figures
beta.post <- extract(kidiq_interaction, "beta")$beta
beta.mean <- colMeans(beta.post)
kidiq.data <- data.frame(kid_score, mom_hs = as.factor(mom_hs), mom_iq)
levels(kidiq.data$mom_hs) <- c("No", "Yes")

# Figure 3.4 (a)
p <- ggplot(kidiq.data, aes(x = mom_iq, y = kid_score, color = mom_hs)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean[1] + beta.mean[2] * (mom_hs == "Yes"),
                  slope = beta.mean[3] + beta.mean[4] * (mom_hs == "Yes"),
                  color = mom_hs)) +
  scale_color_manual("Mother\ncompleted\nhigh\nschool",
                     values = c("No" = "black", "Yes" = "gray")) +
  theme_bw()
print(p +
        scale_x_continuous("Mother IQ score", breaks = seq(80, 140, 20)) +
        scale_y_continuous("Child test score", breaks = seq(20, 140, 40)))

# Figure 3.4 (b)
print(p +
        scale_x_continuous("Mother IQ score", limits = c(0, 150),
                           breaks = seq(0, 150, 50)) +
        scale_y_continuous("Child test score", limits = c(-15, 150),
                           breaks = c(0, 50, 100)))

## 3.4_StatInference ####
## Model (kidiq_multi_preds.stan): kid_score ~ mom_hs + mom_iq 

data.list <- c("N", "kid_score", "mom_hs", "mom_iq")
kidiq_multi_preds <- stan(file='ARM/kidiq_multi_preds.stan', data=data.list,
                          iter=500, chains=4)
kidiq_multi_preds
pairs(kidiq_multi_preds)

## 3.5_GraphDisplays ####
### Regression line as a function of one input variable 
data.list.2 <- c("N", "kid_score", "mom_iq")
stanfit.2 <- stan(file='ARM/kidscore_momiq.stan', data=data.list.2,
                  iter=500, chains=4)
print(stanfit.2, pars = c("beta", "sigma", "lp__"))

# Figure 3.2
beta.post.2 <- extract(stanfit.2, "beta")$beta
beta.mean.2 <- colMeans(beta.post.2)
kidiq.data.2 <- data.frame(kid_score, mom_iq)
p2 <- ggplot(kidiq.data.2, aes(x = mom_iq, y = kid_score)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean.2[1], slope = beta.mean.2[2])) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  theme_bw()
print(p2)

### Two fitted regression lines 

## Model with no interaction: kid_score ~ mom_hs + mom_iq
data.list.3 <- c("N", "kid_score", "mom_hs", "mom_iq")
stanfit.3 <- stan(file='ARM/kidiq_multi_preds.stan', data=data.list.3,
                  iter=500, chains=4)
print(stanfit.3, pars = c("beta", "sigma", "lp__"))

# Figure 3.3
dev.new()
beta.post.3 <- extract(stanfit.3, "beta")$beta
beta.mean.3 <- colMeans(beta.post.3)
kidiq.data.3 <- data.frame(kid_score, mom_hs = as.factor(mom_hs), mom_iq)
levels(kidiq.data.3$mom_hs) <- c("No", "Yes")
p3 <- ggplot(kidiq.data.3, aes(x = mom_iq, y = kid_score, color = mom_hs)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean.3[1] + beta.mean.3[2] * (mom_hs == "Yes"),
                  slope = beta.mean.3[3], color = mom_hs)) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  scale_color_manual("Mother\ncompleted\nhigh\nschool",
                     values = c("No" = "black", "Yes" = "gray")) +
  theme_bw()
print(p3)

## Model with interaction: kid_score ~ mom_hs + mom_iq + mom_hs:mom_iq
stanfit.4 <- stan(file='ARM/kidiq_interaction.stan', data=data.list.3,
                  iter=500, chains=4)
print(stanfit.4, pars = c("beta", "sigma", "lp__"))

# Figure 3.4 (a)
dev.new()
beta.post.4 <- extract(stanfit.4, "beta")$beta
beta.mean.4 <- colMeans(beta.post.4)
p4 <- ggplot(kidiq.data.3, aes(x = mom_iq, y = kid_score, color = mom_hs)) +
  geom_point() +
  geom_abline(aes(intercept = beta.mean.4[1] + beta.mean.4[2] * (mom_hs == "Yes"),
                  slope = beta.mean.4[3] + beta.mean.4[4] * (mom_hs == "Yes"),
                  color = mom_hs)) +
  scale_color_manual("Mother\ncompleted\nhigh\nschool",
                     values = c("No" = "black", "Yes" = "gray")) +
  theme_bw()
print(p4 + scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
        scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)))

# Figure 3.4 (b)
dev.new()
print(p4 +
        scale_x_continuous("Mother IQ score", limits = c(0, 150), breaks = seq(0, 150, 50)) +
        scale_y_continuous("Child test score", limits = c(-15, 150), breaks = c(0, 50, 100)))

### Displaying uncertainty in the fitted regression (Figure 3.10)

dev.new()
n <- 100
ndx <- sample(nrow(beta.post.2), n)
kidiq.post.2 <- data.frame(sampled.int = beta.post.2[ndx,1],
                           sampled.slope = beta.post.2[ndx,2],
                           id = ndx)
p5 <- ggplot(kidiq.data.2, aes(x = mom_iq, y = kid_score)) +
  geom_point() +
  geom_abline(aes(intercept = sampled.int, slope = sampled.slope, group = id),
              data = kidiq.post.2, alpha = 0.05) +
  geom_abline(aes(intercept = beta.mean.2[1], slope = beta.mean.2[2])) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  theme_bw()
print(p5)

### Displaying using one plot for each input variable (Figure 3.11)

dev.new()
n <- 100
ndx <- sample(nrow(beta.post.3), n)
mom_hs.mean <- mean(mom_hs)
mom_iq.mean <- mean(mom_iq)
kidiq.post.3 <- data.frame(
  sampled.hs.int = beta.post.3[ndx, 1] + beta.post.3[ndx, 3] * mom_iq.mean,
  sampled.iq.int = beta.post.3[ndx, 1] + beta.post.3[ndx, 2] * mom_hs.mean,
  sampled.hs.slope = beta.post.3[ndx,2],
  sampled.iq.slope = beta.post.3[ndx,3],
  id = ndx
)

p6 <- ggplot(kidiq.data.3, aes(x = mom_iq, y = kid_score)) +
  geom_point() +
  geom_abline(aes(intercept = sampled.iq.int, slope = sampled.iq.slope, group = id),
              data = kidiq.post.3, alpha = 0.05) +
  geom_abline(aes(intercept = beta.mean.3[1] + beta.mean.3[2] * mom_hs.mean,
                  slope = beta.mean.3[3])) +
  scale_x_continuous("Mother IQ score", breaks = c(80, 100, 120, 140)) +
  scale_y_continuous("Child test score", breaks = c(20, 60, 100, 140)) +
  theme_bw()


p7 <- ggplot(data.frame(kid_score, mom_hs, mom_iq), aes(x = mom_hs, y = kid_score)) +
  geom_jitter(position = position_jitter(width = 0.05, height = 0.05)) +
  geom_abline(aes(intercept = sampled.hs.int, slope = sampled.hs.slope, group = id),
              data = kidiq.post.3, alpha = 0.05) +
  geom_abline(aes(intercept = beta.mean.3[1] + beta.mean.3[3] * mom_iq.mean,
                  slope = beta.mean.3[2])) +
  scale_x_continuous("Mother completed high school", breaks = c(0, 1)) +
  scale_y_continuous("Child test score", breaks = c(0, 50, 100, 150)) +
  theme_bw()

pushViewport(viewport(layout = grid.layout(1, 2)))
print(p6, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p7, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

## 3.6_Diagnostics####
### Model: kid_score ~ mom_iq
data.list.2 <- c("N", "kid_score", "mom_iq")
kidscore_momiq.sf <- stan(file='ARM/kidscore_momiq.stan', data=data.list.2,
                          iter=500, chains=4)
print(kidscore_momiq.sf, pars = c("beta", "sigma", "lp__"))

### Figure 3.12

beta.post <- extract(kidscore_momiq.sf, "beta")$beta
beta.mean <- colMeans(beta.post)
resid <- kid_score - (beta.mean[1] + beta.mean[2] * mom_iq)
resid.sd <- sd(resid)

p <- ggplot(data.frame(mom_iq, resid), aes(x = mom_iq, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-resid.sd, resid.sd), linetype = "dashed") +
  scale_x_continuous("Mother IQ score", breaks = seq(80, 140, 20)) +
  scale_y_continuous("Residuals", breaks = seq(-60, 40, 20))
print(p)

## 3.7_Prediction####
### Model: kid_score ~ mom_hs + mom_iq
data.list <- c("N", "kid_score", "mom_hs", "mom_iq", "mom_hs_new", "mom_iq_new")
kidiq_prediction.sf <- stan(file = 'ARM/kidiq_prediction.stan',
                            data = data.list,
                            iter = 500, chains = 4)
print(kidiq_prediction.sf, pars = c("kid_score_pred"), prob = c(0.025, 0.975))
pairs(kidiq_prediction.sf)

### Model: ppvt ~ hs + afqt
data.list <- c("N", "ppvt", "hs", "afqt")
kidiq_pre1987.sf <- stan(file = 'ARM/kidiq_validation.stan',
                         data = data.list,
                         iter = 500, chains = 4)
print(kidiq_pre1987.sf, pars = c("beta", "sigma", "lp__"))
pairs(kidiq_prediction.sf)

### Data

source("ARM/kids_before1987.data.R", echo = TRUE)

### Model: ppvt ~ hs + afqt
data.list <- c("N", "ppvt", "hs", "afqt")
kidiq_pre1987.sf <- stan(file = 'ARM/kidiq_validation.stan',
                         data = data.list,
                         iter = 500, chains = 4)
print(kidiq_pre1987.sf, pars = c("beta", "sigma", "lp__"))
pairs(kidiq_prediction.sf)

### External validation

## Data

source("ARM/kids_after1987.data.R", echo = TRUE)

## Predicted scores

beta.post <- extract(kidiq_pre1987.sf, "beta")$beta
beta.mean <- colMeans(beta.post)
cscores.new <- beta.mean[1] + beta.mean[2] * hs_ev + beta.mean[3] * afqt_ev
resid <- ppvt_ev - cscores.new
resid.sd <- sd(resid)

## Figure 3.13
# left
p1 <- ggplot(data.frame(cscores.new, ppvt_ev), aes(x = cscores.new, y = ppvt_ev)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous("Predicted score", limits = c(20, 140), breaks = seq(20, 140, 20)) +
  scale_y_continuous("Actual score", limits = c(20, 140), breaks = seq(20, 140, 20)) +
  theme_bw()
print(p1)
# right
dev.new()
p2 <- ggplot(data.frame(cscores.new, resid), aes(x = cscores.new, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-resid.sd, resid.sd), linetype = "dashed") +
  scale_x_continuous("Predicted score", breaks = seq(70, 100, 10)) +
  scale_y_continuous("Prediction error", breaks = seq(-60, 40, 20)) +
  theme_bw()
print(p2)


