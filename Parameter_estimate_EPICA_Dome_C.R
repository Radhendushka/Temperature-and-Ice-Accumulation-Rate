rm(list=ls())
library(Iso)
library(Rfast)
library(matrixStats)
library(gplm)
library(LSTS)
library(forecast)
print(date())
source("All_functions.R")
set.seed(1)
##Reading data file
data=read.csv("EDC_AICC_model2.csv")
## Data with depth larger than 40 meter
Ind=which(data$d>40)
y=data$y[Ind]
x=data$x[Ind]
z=data$z[Ind]

##Estimation of parameters
n=length(z)
h=14
S_Estimates=Smooth_Estimate(h,y,x,z)
tilde_gamma=S_Estimates$tilde_gamma
tilde_log_g=S_Estimates$tilde_log_g
tilde_res=S_Estimates$tilde_res

##ARMA order specification and AAR error fitting
p=4
q=3
err_fit=Arima(tilde_res, order = c(p,0,q),include.mean = F)
innov_est=err_fit$residuals


## Bootstrap se and ci estimation
BS=BS_Estimate(x,z,tilde_gamma,tilde_log_g,tilde_res,p,q)



##Figure 6 (top panel)
par(mar = c(5, 5, 2, 2),cex=1)
plot(z, tilde_log_g,
     type = "l", 
     xaxt = "n", 
     yaxt = "n",
     ylab = "", 
     xlab = "", 
     col = 1, lty = 1,lwd=2,
     xlim = c(0,820),
     ylim = c(-0.5,4)
)
lines(z,BS$BS_log_g[2,],lty=3,lwd=2,col = 1,)
lines(z,BS$BS_log_g[3,],lty=3,lwd=2,col = 1,)
axis(side=2,at = seq(0,4,1),cex.axis=2)
mtext("Smooth log g estimate", side = 2, line = 3,cex=2)
axis(side=1,at = seq(0,820,100),cex.axis=2)
mtext("Age (KyrBP)", side = 1, line = 3,cex=2)
legend("topright",
       c("Smooth log g Estimate ", "95% Confidence limits "),
       col=c(1,1),
       lty = c(1,2),
       bty="n",
       cex=2
)

##Figure 6 bottom panel
par(mar = c(5, 5, 2, 2),cex=1)
plot(z,y,
     type = "l", 
     xaxt = "n", 
     yaxt = "n",
     ylab = "", 
     xlab = "", 
     col = "#00008B", lty = 1,lwd=2,
     xlim = c(0,820),
     ylim = c(-1,4.5)
)
lines(z, log(1+tilde_gamma*x)+tilde_log_g,
      type = "l", 
      xaxt = "n", 
      yaxt = "n",
      ylab = "", 
      xlab = "", 
      col = "#ff8c00", 
      lty = 2,
      lwd=2,
      xlim = c(0,820),
      ylim = c(-1,4.5)
)
axis(side=2,at = seq(-1,4.5,1),cex.axis=2)
mtext("log AAR (m/KYr)", side = 2, line = 3,cex=2)
axis(side=1,at = seq(0,820,100),cex.axis=2)
mtext("Age (KyrBP)", side = 1, line = 3,cex=2)
legend("topright",
       c("log AAR ", "Fitted log AAR"),
       col=c("#00008B","#ff8c00"),
       lty = c(1,2),
       lwd=c(2,2),
       bty="n",
       cex=2
)


##Figure 7 (top panel)
resid=y-tilde_log_g
par(mar = c(5, 5, 2, 2),cex=1)
plot(x,exp(resid),
     xlab=" ",
     ylab=" ",
     xaxt = "n", 
     yaxt = "n",
     xlim=c(-10,5),
     ylim=c(0,2.5),
     pch=19,
     cex=0.8,
     col = gray.colors(5),
)
axis(side=1,at = seq(-10,5,3),cex.axis=2)
mtext("Temperature anomaly", side = 1, line = 3,cex=2)
mtext("AAR / Smooth estimated g", side = 2, line = 3,cex=2)
axis(side=2,at = seq(0,2.5,0.5),cex.axis=2)
f1=lm(exp(resid)~x)
lines(x,f1$fitted.values,col=1,lty=1,lwd=4)
legend("topright",
       c("AAR corrected for thinning", "Fitted line"),
       col=c(gray.colors(5),"black"),
       lty=c(NA,1),
       lwd=c(NA,4),
       pch=c(19,NA),
       pt.cex=c(0.8,NA),
       merge=FALSE,
       bty="n",
       cex=2
)

##Figure 7 (bottom panel)
resid1=y-tilde_log_g-log(1+tilde_gamma*x)
par(mar = c(5, 5, 2, 2),cex=1)
plot(resid1[1:(n-1)],resid1[2:n],
     xlab=" ",
     ylab=" ",
     xaxt = "n", 
     yaxt = "n",
     xlim=c(-0.5,0.5),
     ylim=c(-0.5,0.5),
     pch=19,
     cex=0.8,
     col = gray.colors(5)
)
axis(side=1,at = seq(-0.5,0.5,0.2),cex.axis=2)
mtext("Residual at lag 1", side = 1, line = 3,cex=2)
mtext("Residual", side = 2, line = 3,cex=2)
axis(side=2,at = seq(-0.5,0.5,0.2),cex.axis=2)
points(innov_est[1:(n-1)],innov_est[2:n],
       xlab=" ",
       ylab=" ",
       xaxt = "n", 
       yaxt = "n",
       xlim=c(-0.5,0.5),
       ylim=c(-0.5,0.5),
       pch=21,
       cex=0.8,
       col="darkolivegreen",
       bg = "chocolate1"
)
abline(a=0,b=1,lwd=2,col="#800000")
legend("topleft",
       c("Residuals of AAR", "Predicted innovations of AAR errors"),
       col=c(gray.colors(5),"darkolivegreen"),
       pch=c(19,21),
       pt.bg=c(NA,"chocolate1"),
       pt.cex=c(0.8,0.8),
       merge=FALSE,
       bty="n",
       cex=2
)



##Figure 8 (top panel) 
innov_acf=acf(innov_est,lag=30,ci.type="white",plot=F)
plot(innov_acf$acf[2:31], 
     type="h", 
     xlab=" ",     
     ylab=" ", 
     ylim=c(-0.2,0.2), 
     xaxt="n",
     yaxt="n",lwd=2)
abline(h=0,lwd=2)
abline(h=1.96/sqrt(length(innov_est)),col="blue",lty=2,lwd=2)
abline(h=-1.96/sqrt(length(innov_est)),col="blue",lty=2,lwd=2)
# Add labels to the x-axis
xticks= c(1,seq(5,30, 5))
yticks= seq(-0.2,0.2,0.05)
axis(1, at=xticks, labels=xticks,cex.axis=2)
axis(2, at=yticks, labels=yticks,cex.axis=2)
mtext("Lag", side = 1, line = 3,cex=2)
mtext("Auto-Correlation Function", side = 2, line = 3,cex=2)


#Figure 8 (middel panel)
innov_pacf=pacf(innov_est,lag=30,ci.type="white",plot=F)
plot(innov_pacf$acf, 
     type="h", 
     xlab=" ",     
     ylab=" ", 
     ylim=c(-0.2,0.2), 
     xaxt="n",
     yaxt="n",lwd=2)
abline(h=0,lwd=2)
abline(h=1.96/sqrt(length(innov_est)),col="blue",lty=2,lwd=2)
abline(h=-1.96/sqrt(length(innov_est)),col="blue",lty=2,lwd=2)
# Add labels to the x-axis
xticks= c(1,seq(5,30, 5))
yticks= seq(-0.2,0.2,0.05)
axis(1, at=xticks, labels=xticks,cex.axis=2)
axis(2, at=yticks, labels=yticks,cex.axis=2)
mtext("Lag", side = 1, line = 3,cex=2)
mtext("Partial Auto-Correlation Function", side = 2, line = 3,cex=2)

#Figure8 (bottom panel)

hist(innov_est, breaks = 30,
     probability = T,
     main = " ",
     xlim = c(-0.2,0.2), ylim = c(0,30),
     xlab = "", ylab="",xaxt="n",yaxt="n"
)
axis(1, at=seq(-0.2,0.2,0.1), labels=seq(-0.2,0.2,0.1),cex.axis=2)
axis(2, at=seq(0,30,5), labels=seq(0,30,5),cex.axis=2)
mtext("Residual", side = 1, line = 3,cex=2)
mtext("Density", side = 2, line = 3,cex=2)

print(date())