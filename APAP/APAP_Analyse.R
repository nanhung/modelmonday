rm(list = ls())
system("./mcsim.APAP Forward_APAP1.in")
system("./mcsim.APAP Forward_APAP2.in")
system("./mcsim.APAP Forward_APAP3.in")

APAP_forward1 <- read.delim("Forward_APAP1.out", skip = 2)
APAP_forward2 <- read.delim("Forward_APAP2.out", skip = 2)
APAP_forward3 <- read.delim("Forward_APAP3.out", skip = 2)

par(mfrow = c(3,3))
for (i in 2:4) {
  plot(APAP_forward1$Time, APAP_forward1[,i], xlab = "Time (hr)", ylab = "", 
       main = names(APAP_forward1)[i], las = 1, col = "red", lwd = 2,
       type = "l")
}
for (i in 2:4) {
  plot(APAP_forward2$Time, APAP_forward2[,i], xlab = "Time (hr)", ylab = "",
       main = names(APAP_forward2)[i], las = 1, col = "red", lwd = 2,
       type = "l")
}
for (i in 2:4) {
  plot(APAP_forward3$Time, APAP_forward3[,i], xlab = "Time (hr)", ylab = "",
       main = names(APAP_forward3)[i], las = 1, col = "red", lwd = 2,
       type = "l")
}
