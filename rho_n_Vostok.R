rm(list=ls())
library(Iso)
library(Rfast)
library(matrixStats)
library(gplm)
source("All_functions.R")

data=read.csv("Vostok_AICC_model2.csv")

Ind=which(data$d>40)
y=data$y[Ind]
x=data$x[Ind]
z=data$z[Ind]

n=length(y)
h=14
beta_rho_n=c()
n1=c()
l=0
for(i in 1:8)
{
  for (j in 1:i)
  {
    sub=seq(j,n,i)
    y1=y[sub]
    x1=x[sub]
    z1=z[sub]
    n1[l+1]=length(y1)
    LS=LS_Estimate(y1,x1,z1)
    tilde_log_g1=kreg(z1, LS$log_g_hat, kernel = c("epanechnikov"), bandwidth = h,grid=z1)$y  
    model=nls(y1-tilde_log_g1~log(1+v*x1),start = list(v=LS$gamma_hat) )
    tilde_gamma1=as.numeric(coef(model))
    tilde_res=y1-log(1+tilde_gamma1*x1)-tilde_log_g1
    beta_rho_n[l+1]=-log(sum(tilde_res[2:n1[l+1]]*tilde_res[1:(n1[l+1]-1)])/sum(tilde_res[1:n1[l+1]]^2))
    l=l+1
  }
}
M=lm(log(beta_rho_n)~log(n1))
summary(M)

##Figure 6 
par(mar = c(5, 5, 2, 2),cex=1)
plot(log(n1),log(beta_rho_n),xlab="",ylab="",pch=21,
     cex=2,
     col = 1,bg="#AAAAAA",xlim=c(6,9),ylim=c(-4.5,-1.5),xaxt = "n", 
     yaxt = "n")
lines(log(n1),predict(M),lty=1,lwd=2)
mtext("Estimated log(βρ)",side=2,line=3,cex=2)
mtext("log(n)", side=1,line=3,cex=2 )
axis(side=1,at = seq(6,9,0.5),cex.axis=2)
axis(side=2,at = seq(-4.5,-1.5,0.5),cex.axis=2)
legend("bottomright",
       c( "Fitted line","Estimated log(βρ)"),
       col=c(1,"black"),
       pt.bg=c(NA,"#AAAAAA"),
       lty=c(1,NA),
       lwd=c(2,NA),
       pch=c(NA,21),
       pt.cex=c(NA,2),
       merge=FALSE,
       bty="n",
       cex=2
)
