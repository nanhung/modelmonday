if(!require(pksensi)) install.packages("pksensi"); library(pksensi)

mName <- "pbtk1cpt_v2"
compile_model(mName, application = "mcsim",　version = "6.0.1")
parms <- c(vdist = 0.5, ke = 0.2, km = 0.5, kgutabs = 2.0)
times <- seq(from = 0.01, to = 24.01, by = 1)

#
LL <- 0.5 
UL <- 1.5
q <- "qunif"
q.arg <- list(list(min = parms["vdist"] * LL, max = parms["vdist"] * UL),
              list(min = parms["ke"] * LL, max = parms["ke"] * UL),
              list(min = parms["km"] * LL, max = parms["km"] * UL),
              list(min = parms["kgutabs"] * LL, max = parms["kgutabs"] * UL)) 
set.seed(1234)
params <- names(parms)
x <- rfast99(params, n = 800, q = q, q.arg = q.arg, replicate = 20)
cex <- 0.2

png(file="png/sim2.png",width=3300,height=2100,res=300)
par(mfrow=c(4,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,1,"vdist"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"vdist"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"ke"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,2,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,3,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,4,"km"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], ylab = "", cex = cex)
plot(x$a[,2,"kgutabs"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,3,"kgutabs"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,4,"kgutabs"], ylab = "", yaxt="n", cex = cex)
dev.off()

#
Outputs <- c("Agutlument", "Aelimination", "Acompartment", "Ccompartment", "AUC", "Ametabolized")
conditions <- c("Agutlument = 10") # Set the initial state of Agutlument = 10 
y<-solve_mcsim(x, mName = mName, 
               params = params,
               vars = Outputs,
               time = times,
               condition = conditions)
tell2(x,y)

png(file="png/sim3.png",width=3300,height=1800,res=300)
par(mfrow = c(2,3), mar = c(2,2,2,1), oma = c(2,2,0,0))
pksim(y, vars = "Agutlument", main = "Agutlument")
pksim(y, vars = "Aelimination", legend = F, main = "Aelimination")
pksim(y, vars = "Acompartment", legend = F, main = "Acompartment")
pksim(y, vars = "Ccompartment", legend = F, main = "Ccompartment")
pksim(y, vars = "Ametabolized", legend = F, main = "Ametabolized")
pksim(y, vars = "AUC", legend = F, main = "AUC")
mtext("Time", SOUTH<-1, line=0.4, outer=TRUE)
mtext("Quantity", WEST<-2, line=0.4, outer=TRUE)
dev.off()

#
png(file="png/sim4.png",width=3300,height=1800,res=300)
plot(x, var = 4)
dev.off()

png(file="png/sim5.png",width=3300,height=1800,res=300)
plot(x, var = 6)
dev.off()


png(file="png/sim6.png",width=3300,height=2100,res=300)
heat_check(x)
dev.off()

png(file="png/sim7.png",width=3300,height=2100,res=300)
heat_check(x, index = "CI")
dev.off()


check(x, SI.cutoff = 0.05, vars = "Ccompartment")
check(x, SI.cutoff = 0.05, vars = "Ametabolized")
check(x, SI.cutoff = 0.05)


#
png(file="png/sim8.png",width=3300,height=2100,res=300)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,1,"vdist"], y[,1,"0.01","Ccompartment"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"0.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"0.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"0.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"1.01","Ccompartment"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"1.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"1.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"1.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"2.01","Ccompartment"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"2.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"2.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"2.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"4.01","Ccompartment"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"4.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"4.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"4.01","Ccompartment"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"12.01","Ccompartment"], ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"12.01","Ccompartment"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,1,"km"], y[,1,"12.01","Ccompartment"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"12.01","Ccompartment"], ylab = "", yaxt="n", cex = cex)
mtext("Ccompartment", NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)
dev.off()

png(file="png/sim9.png",width=3300,height=2100,res=300)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(3,3,2,1));
plot(x$a[,1,"vdist"], y[,1,"0.01","Ametabolized"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"0.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"0.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"0.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"1.01","Ametabolized"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"1.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"1.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"1.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"2.01","Ametabolized"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"2.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"2.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"2.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"4.01","Ametabolized"], xaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"4.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"km"], y[,1,"4.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"4.01","Ametabolized"], xaxt="n", yaxt="n", ylab = "", cex = cex)
plot(x$a[,1,"vdist"], y[,1,"12.01","Ametabolized"], ylab = "", cex = cex)
plot(x$a[,1,"ke"], y[,1,"12.01","Ametabolized"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,1,"km"], y[,1,"12.01","Ametabolized"], ylab = "", yaxt="n", cex = cex)
plot(x$a[,1,"kgutabs"], y[,1,"12.01","Ametabolized"], ylab = "", yaxt="n", cex = cex)
mtext("Ametabolized", NORTH<-3, line=0.4, adj=0, cex=1.5, outer=TRUE)
dev.off()
