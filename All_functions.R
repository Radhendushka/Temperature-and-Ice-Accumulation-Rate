LS_Estimate=function(y,x,z)
{
  log_g_hat=pava(y,decreasing =  T)
  model=nls(y-log_g_hat~log(1+v*x),start = list(v=sum((exp(y-log_g_hat)-1)*x)/sum(x*x)) )
  gamma_hat=as.numeric(coef(model))
  Q=mean((y-log(1+gamma_hat*x)-log_g_hat)^2)
  Criterion=1
  while(Criterion>10^(-10))
  {
    log_g_hat=pava(y-log(1+gamma_hat*x),decreasing =  T)
    model=nls(y-log_g_hat~log(1+v*x),start = list(v=gamma_hat))
    gamma_hat=as.numeric(coef(model))  
    Q_update=mean((y-log(1+gamma_hat*x)-log_g_hat)^2)
    Criterion=abs(Q_update-Q)
    Q=Q_update
  }
  Estimates=list(gamma_hat=gamma_hat,
                 log_g_hat=log_g_hat,
                 Avg_LS=Q)
  return(Estimates)
  
}


Smooth_Estimate=function(h,y,x,z,rho_n)
{
  LS=LS_Estimate(y,x,z)
  n=length(y)
  tilde_log_g=kreg(z, LS$log_g_hat, kernel = c("epanechnikov"), bandwidth = h,grid=z)$y  
  model=nls(y-tilde_log_g~log(1+v*x),start = list(v=LS$gamma_hat) )
  tilde_gamma=as.numeric(coef(model))
  tilde_res=y-log(1+tilde_gamma*x)-tilde_log_g
  tilde_beta=-(1/rho_n)*log(sum(tilde_res[2:n]*tilde_res[1:(n-1)])/sum(tilde_res[1:n]^2))
  tilde_sig= sqrt(mean(tilde_res[1:n]^2)*2*tilde_beta)
  Estimates=list(tilde_gamma=tilde_gamma,
                 tilde_log_g=tilde_log_g,
                 tilde_beta=tilde_beta,
                 tilde_sig=tilde_sig)
  return(Estimates)
  
}  



BS_Estimate=function(y,x,z,tilde_gamma,tilde_log_g,tilde_beta,rho_n)
{
  
  n=length(z)
  M=1000
  ## Let theta=(tilde_gamma, tilde_beta,tilde_sig2)
  tilde_theta_star=c()
  tilde_log_g_star=c()
  tilde_res=y-log(1+tilde_gamma*x)-tilde_log_g
  tilde_U=tilde_res[2:n]-exp(-tilde_beta*rho_n)*tilde_res[1:(n-1)]
  s=0
  while (s <M)
  {
  U_star=sample(tilde_U,n,replace = T)
  if(var(U_star)>=var(tilde_U))
  {  
  eps_star=arima.sim(list(order=c(1,0,0),ar=exp(-tilde_beta*rho_n)),
                     n,innov =  U_star,
                     n.start=1000,
                     start.innov=sample(tilde_U,1000,replace = T))
  y_star=log(1+tilde_gamma*x)+tilde_log_g+eps_star
  E_star=Smooth_Estimate(h,y_star,x,z,rho_n)
  tilde_theta_star=rbind(c(tilde_theta_star, 
                                     c(E_star$tilde_gamma,
                                       E_star$tilde_beta,
                                       E_star$tilde_sig)))
  tilde_log_g_star=rbind(c(tilde_log_g_star,
                           c(E_star$tilde_log_g)))
  s=s+1
  #print(s)
  }
  
  
  }
  
  Mat_theta_BS=matrix(c(tilde_theta_star),ncol=3,byrow=T)
  Mat_log_g_BS=matrix(c(tilde_log_g_star),ncol=n,byrow=T)
  
  
  
  theta_quantile_025=colQuantiles(Mat_theta_BS,probs=0.025)
  theta_quantile_975=colQuantiles(Mat_theta_BS,probs=0.975)
  theta_sd=(colVars(Mat_theta_BS))^(0.5)
  BS_theta=matrix(rbind(c(theta_sd,
                  theta_quantile_025,
                  theta_quantile_975)),ncol=3,byrow=T)
  
  log_g_quantile_025=colQuantiles(Mat_log_g_BS,probs=0.025)
  log_g_quantile_975=colQuantiles(Mat_log_g_BS,probs=0.975)
  log_g_sd=(colVars(Mat_log_g_BS))^(0.5)
  BS_log_g=matrix(rbind(c(log_g_sd,
                         log_g_quantile_025,
                         log_g_quantile_975)),ncol=n,byrow=T)
  
  
  
  Estimates=list(BS_theta=BS_theta,
                 BS_log_g=BS_log_g)
  
  return(Estimates)
  
}
