# http://wiekvoet.blogspot.com/2015/04/hierarchical-two-compartimental-pk-model.html
rm(list=ls())

Theoph.1 <- Theoph[ Theoph$Subject == 1, ]
with(Theoph.1, SSfol(Dose, Time, -2.5, 0.5, -3)) # response only
with(Theoph.1, local({  lKe <- -2.5; lKa <- 0.5; lCl <- -3
SSfol(Dose, Time, lKe, lKa, lCl) # response _and_ gradient
}))
getInitial(conc ~ SSfol(Dose, Time, lKe, lKa, lCl), data = Theoph.1)
## Initial values are in fact the converged values
fm1 <- nls(conc ~ SSfol(Dose, Time, lKe, lKa, lCl), data = Theoph.1)
summary(fm1)

coplot(conc ~ Time | Subject, data = Theoph, show.given = FALSE)
Theoph.4 <- subset(Theoph, Subject == 4)
fm1 <- nls(conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
           data = Theoph.4)
summary(fm1)
plot(conc ~ Time, data = Theoph.4,
     xlab = "Time since drug administration (hr)",
     ylab = "Theophylline concentration (mg/L)",
     main = "Observed concentrations and fitted model",
     sub  = "Theophylline data - Subject 4 only",
     las = 1, col = 4)
xvals <- seq(0, par("usr")[2], length.out = 55)
lines(xvals, predict(fm1, newdata = list(Time = xvals)),
      col = 4)



library(lattice)
library(latticeExtra)
library(tidyverse)

Theoph.MW <- 180.164 #g/mol
Theoph$conc.mol<-Theoph$conc * 1000 / Theoph.MW
Theoph$dose.mol <- Theoph$Dose * Theoph$Wt * 10^6 / Theoph.MW # uM

par(mfrow=c(3,4))
for (i in 0:11){
  plot(Theoph[1:11,4], Theoph[i*11+c(1:11), 6], ylim=c(1,100), 
       main = paste("Ind.", LETTERS[1+i], ":", " ", " ", mean(Theoph[i*11+c(1:11),3]), "mg/kg", sep = ""),
       xlab="Time (hr)", ylab=expression(paste("Plasma concentration, ",mu, "M")), type="b")
}

xyplot(conc ~ Time| Subject, data=Theoph, xlab="Time (hr)",
       auto.key = T, ylab=expression(paste("Plasma concentration, ",mu, "M")))
theme_set(theme_light())
ggplot(Theoph, aes(Time, conc)) + geom_point() + facet_wrap(~Subject)

la <- lapply(levels(Theoph$Subject), function(iSubject)
{Theoph.i <- subset(Theoph, Subject == iSubject) 
  fm1 <- nls(conc ~ SSfol(Dose, Time, lKe, lKa, lCl), data = Theoph.i)
  ee <- exp(coef(fm1))
  dd <- data.frame(Ke=ee[1],Ka=ee[2],CL=ee[3])
  row.names(dd) <- iSubject
  dd
}
)

nlsres <- do.call(rbind,la)
nlsres

for ( i in 1:12){
  Theoph.n <- Theoph[ Theoph$Subject == i, ]
  Ke <- nlsres[paste(i), "Ke"]
  Ka <- nlsres[paste(i), "Ka"]
  CL <- nlsres[paste(i), "CL"] 
  Time <- Theoph.n$Time
  x <- Theoph.n$Dose*((Ke*Ka)/CL)/(Ka-Ke) * (exp(-Ke*Time)-exp(-Ka*Time))
  if (i == 1) Prediction <- x else Prediction <- c(Prediction, x)
}

Theoph$Prediction <- Prediction

xyplot(Prediction ~ Time | Subject, data=Theoph, type="l") + 
  as.layer(xyplot(conc~Time | Subject, pch=16, col="black", data=Theoph))
ggplot(Theoph, aes(Time, conc)) + geom_point() + geom_line(aes(Time, Prediction)) + 
  facet_wrap(~Subject)

# http://wiekvoet.blogspot.com/2015/04/pk-analysis-theoph-again.html
library(tidyverse)
library(rstan) 

Theoph$dose <- Theoph$Dose*Theoph$Wt
datain <-  list(
  Time=Theoph$Time,
  conc=Theoph$conc,
  n=nrow(Theoph),
  subject =c(1:nlevels(Theoph$Subject))[Theoph$Subject],
  nsubject = nlevels(Theoph$Subject),
  dose=unique(subset(Theoph,, c(Subject,dose)))$dose)

smodel <- '
data {
int <lower=1> n;
int <lower=1> nsubject;
int <lower=1, upper=nsubject>  subject[n];
real  dose[nsubject];
real  Time[n];
real  conc[n];
}
parameters {
real <lower=0> sdr;
real <lower=0> kas[nsubject];
real <lower=0> kes[nsubject];
real <lower=0> CLs[nsubject];
real lke;
real lka;
real <lower=0> CL;
real kesd;
real kasd;
real CLsd;
}

transformed parameters {
real pred[n];
real <lower=1> c0star[nsubject];

for (i in 1:nsubject) {
if (kas[i]>kes[i]){
  c0star[i] <- dose[i]*((kes[i]*kas[i])/CLs[i])/(kas[i]-kes[i]);
} else c0star[i] <- 1;
}
for (i in 1:n) {
pred[i] <- c0star[subject[i]]*(exp(-kes[subject[i]]*Time[i])-exp(-kas[subject[i]]*Time[i]));
}
}

model {
    conc ~ normal(pred,sdr);
    kes ~ lognormal(lke,kesd);
    kas ~ lognormal(lka,kasd);
    lke ~ uniform(-3,3);
    lka ~ uniform(-3,3);
    CL ~ uniform(0.01,300);
    CLs ~ normal(CL,CLsd);
    kesd ~ uniform(0.01,2);
    kasd ~ uniform(0.01,2);
    CLsd ~ uniform(0.01,10);
}
generated quantities {
real ka;
real ke;
ka <- exp(lka);
ke <- exp(lke);
}' 

init <- function(){
  list(lka = rnorm(1,log(2),.3),
       kas=   rnorm(12,2,.2),
       lke=   rnorm(1,log(.08),.01),
       kes=   rnorm(12,.08,.01),
       CLs = rnorm(12,1,.1),
       CL = rnorm(1,1,.1),
       kesd=runif(1,0.04,0.1),
       kasd=runif(1,0.2,2),
       sdr= runif(1,1,3),
       CLsd=runif(1,.1,1))
}

#(kes[i]<kas[i])* 
fstan <- stan(model_code = smodel, 
              data = datain, 
              pars=c('CL','kes','kas','CLs',
                     'c0star','ka','ke'),
              chains = 4,
              init = init)
fstan

stansampls <- rstan::extract(fstan, c('CLs[7]','kes[7]','kas[7]'))
Time <- seq(1,max(datain$Time))

la2 <- sapply(1:nrow(stansampls$CL),
              function(i) {
                CL <- stansampls$CL[i]
                Ka <- stansampls$ka[i]
                Ke <- stansampls$ke[i]
                c0star <- mean(datain$dose)*((Ke*Ka)/CL)/(Ka-Ke)
                y<- c0star* (exp(-Ke*Time)-exp(-Ka*Time))
              })

res2 <- data.frame(
  Time=Time,
  Conc=apply(la2,1,mean),
  C05=apply(la2,1,FUN=function(x) quantile(x,0.05)),
  C95=apply(la2,1,FUN=function(x) quantile(x,0.95)))

ggplot(res2, aes(x=Time)) +
  geom_line(aes(y=Conc))+
  scale_y_log10('theophylline (mg/L)')+
  ylab("theophylline (mg/L)")+
  geom_point(aes(y=conc,x=Time,col=Subject),alpha=1,data=Theoph.1)+
  geom_ribbon(aes(ymin=C05, ymax=C95),alpha=.2)+
  theme(legend.position='bottom')

#

flist <- rstan::extract(fstan,c("CLs",'kes','kas','c0star'))

CLs<-flist[[1]]
kes<-flist[[2]]
kas<-flist[[3]]
c0star<-flist[[4]]

subject<-function(par1,par2,par3,par4){
  stansampls.3 <- rstan::extract(fstan,c(par1,par2,par3,par4))
  la3 <- as.data.frame(c(stansampls.3[1], stansampls.3[2], stansampls.3[3], stansampls.3[4]))
  
  y.01 <- la3[,4]* (exp(-la3[,2]*0.1)-exp(-la3[,3]*0.1)) * 1000 / Theoph.MW
  y.02 <- la3[,4]* (exp(-la3[,2]*0.25)-exp(-la3[,3]*0.25)) * 1000 / Theoph.MW
  y.03 <- la3[,4]* (exp(-la3[,2]*0.5)-exp(-la3[,3]*0.5)) * 1000 / Theoph.MW
  y.04 <- la3[,4]* (exp(-la3[,2]*1)-exp(-la3[,3]*1)) * 1000 / Theoph.MW
  y.05 <- la3[,4]* (exp(-la3[,2]*1.5)-exp(-la3[,3]*1.5)) * 1000 / Theoph.MW
  y.06 <- la3[,4]* (exp(-la3[,2]*2)-exp(-la3[,3]*2)) * 1000 / Theoph.MW
  y.07 <- la3[,4]* (exp(-la3[,2]*3)-exp(-la3[,3]*3)) * 1000 / Theoph.MW
  y.08 <- la3[,4]* (exp(-la3[,2]*4)-exp(-la3[,3]*4)) * 1000 / Theoph.MW
  y.09 <- la3[,4]* (exp(-la3[,2]*5)-exp(-la3[,3]*5)) * 1000 / Theoph.MW
  y.10 <- la3[,4]* (exp(-la3[,2]*6)-exp(-la3[,3]*6)) * 1000 / Theoph.MW
  y.11 <- la3[,4]* (exp(-la3[,2]*8)-exp(-la3[,3]*8)) * 1000 / Theoph.MW
  y.12 <- la3[,4]* (exp(-la3[,2]*10)-exp(-la3[,3]*10)) * 1000 / Theoph.MW
  y.13 <- la3[,4]* (exp(-la3[,2]*12)-exp(-la3[,3]*12)) * 1000 / Theoph.MW
  y.14 <- la3[,4]* (exp(-la3[,2]*16)-exp(-la3[,3]*16)) * 1000 / Theoph.MW
  y.15 <- la3[,4]* (exp(-la3[,2]*20)-exp(-la3[,3]*20)) * 1000 / Theoph.MW
  y.16 <- la3[,4]* (exp(-la3[,2]*25)-exp(-la3[,3]*25)) * 1000 / Theoph.MW
  
  do.call(cbind, list(la3,y.01,y.02,y.03,y.04,y.05,y.06,y.07,y.08,y.09,y.10,
                      y.11,y.12,y.13,y.14,y.15,y.16))
}

sub.1<-subject('CLs[1]','kes[1]','kas[1]','c0star[1]')
sub.2<-subject('CLs[2]','kes[2]','kas[2]','c0star[2]')
sub.3<-subject('CLs[3]','kes[3]','kas[3]','c0star[3]')
sub.4<-subject('CLs[4]','kes[4]','kas[4]','c0star[4]')
sub.5<-subject('CLs[5]','kes[5]','kas[5]','c0star[5]')
sub.6<-subject('CLs[6]','kes[6]','kas[6]','c0star[6]')
sub.7<-subject('CLs[7]','kes[7]','kas[7]','c0star[7]')
sub.8<-subject('CLs[8]','kes[8]','kas[8]','c0star[8]')
sub.9<-subject('CLs[9]','kes[9]','kas[9]','c0star[9]')
sub.10<-subject('CLs[10]','kes[10]','kas[10]','c0star[10]')
sub.11<-subject('CLs[11]','kes[11]','kas[11]','c0star[11]')
sub.12<-subject('CLs[12]','kes[12]','kas[12]','c0star[12]')

t <- c(0.1,0.25,0.5,1,1.5,2,3,4,5,6,8,10,12,16,20,25)

#png(file="0505_1.7.png",width=3000,height=2000,res=300)
par(mfrow=c(3,4), mar=c(3,2,3,1))
plot(t, apply(sub.1[,5:20], 2, quantile, probs= c(.95)), 
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. A")
lines(t, apply(sub.1[,5:20], 2,  median),lty=1, lwd=2)
lines(t, apply(sub.1[,5:20], 2,  quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[1:11,4],Theoph[1:11,7], col = "red", pch=19)
plot(t, apply(sub.2[,5:20], 2, quantile, probs= c(.95)), 
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. B")
lines(t, apply(sub.2[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.2[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[1*11+c(1:11),4],Theoph[1*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.3[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. C")
lines(t, apply(sub.3[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.3[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[2*11+c(1:11),4],Theoph[2*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.4[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. D")
lines(t, apply(sub.4[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.4[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[3*11+c(1:11),4],Theoph[3*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.5[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. E")
lines(t, apply(sub.5[,5:20], 2,  median),lty=1, lwd=2)
lines(t, apply(sub.5[,5:20], 2,  quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[4*11+c(1:11),4],Theoph[4*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.6[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. F")
lines(t, apply(sub.6[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.6[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[5*11+c(1:11),4],Theoph[5*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.7[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. G")
lines(t, apply(sub.7[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.7[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[6*11+c(1:11),4],Theoph[6*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.8[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. H")
lines(t, apply(sub.8[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.8[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[7*11+c(1:11),4],Theoph[7*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.9[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. I")
lines(t, apply(sub.9[,5:20], 2,  median),lty=1, lwd=2)
lines(t, apply(sub.9[,5:20], 2,  quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[8*11+c(1:11),4],Theoph[8*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.10[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. J")
lines(t, apply(sub.10[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.10[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[9*11+c(1:11),4],Theoph[9*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.11[,5:20], 2, quantile, probs= c(.95))
     , type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. K")
lines(t, apply(sub.11[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.11[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[10*11+c(1:11),4],Theoph[10*11+c(1:11),7], col = "red", pch=19)
plot(t, apply(sub.12[,5:20], 2, quantile, probs= c(.95)),
     type="l", lwd=1, lty=2, log="y", ylim=c(1,100),
     main="Ind. L")
lines(t, apply(sub.12[,5:20], 2, median),lty=1, lwd=2)
lines(t, apply(sub.12[,5:20], 2, quantile, probs= c(.05)),lty=2, lwd=1 )
points(Theoph[11*11+c(1:11),4],Theoph[11*11+c(1:11),7], col = "red", pch=19)
#dev.off()


posterior <- as.array(fstan)
mcmc_trace(posterior, pars = c("ka"))

library("bayesplot")
library("rstanarm")

#plot(fstan,pars="kas")
#plot(fstan,pars="kes")
posterior <- as.matrix(fstan)

#png(file="0505_1.8.png",width=3000,height=2000,res=300)
mcmc_areas(posterior, 
           pars = c("kas[1]","kas[2]","kas[3]","kas[4]","kas[5]","kas[6]",
                    "kas[7]","kas[8]","kas[9]","kas[10]","kas[11]","kas[12]"), 
           prob = 0.95) + ggtitle("Absorption rate (/hr) posterior distributions",
                                  "with medians and 95% intervals")
#dev.off()

#png(file="0505_1.9.png",width=3000,height=2000,res=300)
mcmc_areas(posterior, 
           pars = c("kes[1]","kes[2]","kes[3]","kes[4]","kes[5]","kes[6]",
                    "kes[7]","kes[8]","kes[9]","kes[10]","kes[11]","kes[12]"), 
           prob = 0.95) + ggtitle("Elimination rate (/hr) posterior distributions",
                                  "with medians and 95% intervals")
#dev.off()

#png(file="0505_1.91.png",width=3000,height=2000,res=300)
pairs(fstan, pars=c("ka","ke", "CL"))
#dev.off()

##########################

library("bayesplot")
x <- example_mcmc_draws(chains = 4, params = 6)

dim(x)
color_scheme_set("green")
dimnames(x)

parameter<-dimnames(x)[[3]]
mcmc_acf(x, pars = parameter)





  