library(pksensi)

mName <- "pbtk1cpt"
compile_model(mName, application = "R")
#> * Created file 'pbtk1cpt.so'.
source(paste0(mName, "_inits.R"))

parms <- initParms()
parms["vdist"] <- 0.5
parms["ke"] <- 0.2
parms["km"] <- 0.5
parms["kgutabs"] <- 2.0
initState <- initStates(parms=parms)
initState["Agutlument"] <- 10

parms

initState

Outputs

times <- seq(from = 0.01, to = 24.01, by = 1)


y <- deSolve::ode(initState, times, func = "derivs", parms = parms, 
                  dllname = mName, initfunc = "initmod", nout = 1, outnames = Outputs)

#
par(mar = c(3,3,2,2))
png(file="png/sim1.png",width=2400,height=1500,res=300)
plot(y)
dev.off()

LL <- 0.5 
UL <- 1.5
q <- "qunif"
q.arg <- list(list(min = parms["vdist"] * LL, max = parms["vdist"] * UL),
              list(min = parms["ke"] * LL, max = parms["ke"] * UL),
              list(min = parms["km"] * LL, max = parms["km"] * UL),
              list(min = parms["kgutabs"] * LL, max = parms["kgutabs"] * UL)) 
set.seed(1234)
params <- c("vdist", "ke", "km", "kgutabs")
x <- rfast99(params, n = 800, q = q, q.arg = q.arg, replicate = 10)

#
par(mfrow = c(2,5), mar = c(0,3,2,2))


for(i in 1:10){
  plot(x$a[,i,"kgutabs"], type = "b")
}



Outputs <- c("Ccompartment", "Agutlument", "Ametabolized", "Aelimination")
y <- solve_fun(x, times, initState = initState, outnames = Outputs, dllname = mName)
tell2(x,y)


par(mfrow = c(2,2), mar = c(2,2,2,1), oma = c(2,2,0,0))
pksim(y, vars = "Agutlument", main = "Agutlument")
pksim(y, vars = "Aelimination", legend = F, main = "Aelimination")
pksim(y, vars = "Ccompartment", legend = F, main = "Ccompartment")
pksim(y, vars = "Ametabolized", legend = F, main = "Ametabolized")
mtext("Time", SOUTH<-1, line=0.4, outer=TRUE)
mtext("Quantity", WEST<-2, line=0.4, outer=TRUE)


#
plot(x, vars = "Ccompartment")
plot(x, vars = "Ametabolized")
check(x, SI.cutoff = 0.5, vars = "Ccompartment")
check(x, SI.cutoff = 0.5, vars = "Ametabolized")
heat_check(x)
heat_check(x, index = "CI")

conditions <- c("Agutlument = 10") # Set the initial state of Agutlument = 10 
y<-solve_mcsim(x, mName = mName, 
               params = params,
               vars = Outputs,
               time = times,
               condition = conditions)




#
par(mfrow = c(3,3), mar = c(2,2,2,2), oma = c(2,2,1,1))
plot(x$a[,1,"vdist"], y[,1,"0.01","Ccompartment"], main = "vdist")
text(0.7, .7, "t=0.01",cex = 1.2)
plot(x$a[,1,"ke"], y[,1,"0.01","Ccompartment"], main = "ke")
plot(x$a[,1,"kgutabs"], y[,1,"0.01","Ccompartment"], main = "kgutabs")
plot(x$a[,1,"vdist"], y[,1,"2.01","Ccompartment"])
text(1, 18, "t=2.01",cex = 1.2)
plot(x$a[,1,"ke"], y[,1,"2.01","Ccompartment"])
plot(x$a[,1,"kgutabs"], y[,1,"2.01","Ccompartment"])
plot(x$a[,1,"vdist"], y[,1,"24.01","Ccompartment"])
text(1, .7, "t=24.01",cex = 1.2)
plot(x$a[,1,"ke"], y[,1,"24.01","Ccompartment"])
plot(x$a[,1,"kgutabs"], y[,1,"24.01","Ccompartment"])
mtext("parameter", SOUTH<-1, line=0.4, outer=TRUE)
mtext("Ccompartment", WEST<-2, line=0.4, outer=TRUE)

