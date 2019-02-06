#https://mc-stan.org/workshops/Rmedicine2018/

library(rstanarm)
post <- stan_nlmer(conc ~ SSfol(Dose, Time, lKe, lKa, lCl) ~ 
                      (0 + lKe + lKa + lCl | Subject), data = Theoph,
                    prior = normal(location = c(-2, 0.5, -3), scale = 1, autoscale = FALSE),
                    seed = 982018, refresh = 0)
post