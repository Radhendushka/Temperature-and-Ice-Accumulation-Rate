## The function LS_Estimate computes the estimates for
## $\hat{gamma}$ and $\log \hat{g}$ given in section 2.
## Input: log(AAR):= y; temperature:= x; and age:= z
## Output: gamma_hat:= $\hat{gamma}$; log_g_hat:=$\log \hat{g}$, Q_value:=$Q(\hat{gamma},\hat{g})$

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
                 Q_value=Q)
  return(Estimates)
  
}


## The function Smooth_Estimate computes the estimates
## $\tilde{gamma}$ and $\log {\tilde{g}}$ given in section 2.
## It also computes log AAR residual $\tilde{\varepsilon}$ given in section 3.
## Input: Smoothing parameter:=h; log(AAR):=y; temperature:=x; and age:=z
## Output: tilde_gamma:=$\tilde{gamma}$; tilde_log_g:=$\log {\tilde{g}}$; tilde_res:=$\tilde{\varepsilon}$

Smooth_Estimate=function(h,y,x,z)
{
  LS=LS_Estimate(y,x,z)
  n=length(y)
  tilde_log_g=kreg(z, LS$log_g_hat, kernel = c("epanechnikov"), bandwidth = h,grid=z)$y  
  model=nls(y-tilde_log_g~log(1+v*x),start = list(v=LS$gamma_hat) )
  tilde_gamma=as.numeric(coef(model))
  tilde_res=y-log(1+tilde_gamma*x)-tilde_log_g
  
  Estimates=list(tilde_gamma=tilde_gamma,
                 tilde_log_g=tilde_log_g,
                 tilde_res=tilde_res
                 )
  return(Estimates)
  
}  




## The function BS_Estimate computes the bootstrap standard error, bootstarp quantiles of 
## $\tilde{gamma}$ and $\log {\tilde{g}}$ given in section 4.
## Input:  temperature:=x; age:=z; estimated gamma:=tilde_gamma; estimated log g:=tilde_log_g; log AAR residual:=tilde_res AR order:=p; MA order:=q
## Output: BS_gamma:= A matrix of order 3 times 1: first row gives bootstrap se, 2.5% bootstrap quantile and 97.5% bootstrap quantile of $\tilde{gamma}$ are given in second and third row, 
##         BS_log:= A matrix of size 3 times n: The first row gives the bootstrap se of $\log \tilde{g}$ 
##                  whereas the second and third row gives 2.5% bootstrap quantile and 97.5% bootstrap quantile of $\log \tilde{g}$. The columns correspond to age values.
BS_Estimate=function(x,z,tilde_gamma,tilde_log_g,tilde_res,p,q)
{
  err_fit=Arima(tilde_res, order = c(p,0,q),include.mean = F)
  ar_coef=err_fit$coef[1:p]
  ma_coef=err_fit$coef[(p+1):(p+q)]
  innov_est=err_fit$residuals
  n=length(z)
  M=1000
  tilde_gamma_star=c()
  tilde_log_g_star=c()
  s=0
  while (s <M)
  {
  innov_est_star=sample(innov_est,n,replace = T)  
  if(var(innov_est_star)>=var(innov_est))
  {  
  eps_star=arima.sim(list(ar=ar_coef,ma=ma_coef),
                         n,innov =  innov_est_star,
                         n.start=1000,
                         start.innov=sample(innov_est_star,1000,replace = T))
  y_star=log(1+tilde_gamma*x)+tilde_log_g+eps_star
  E_star=Smooth_Estimate(h,y_star,x,z)
  tilde_gamma_star=rbind(c(tilde_gamma_star, c(E_star$tilde_gamma)))
  tilde_log_g_star=rbind(c(tilde_log_g_star,c(E_star$tilde_log_g)))
  s=s+1
  }
  }
  Mat_gamma_BS=matrix(c(tilde_gamma_star),ncol=1,byrow=T)
  Mat_log_g_BS=matrix(c(tilde_log_g_star),ncol=n,byrow=T)
  gamma_quantile_025=colQuantiles(Mat_gamma_BS,probs=0.025)
  gamma_quantile_975=colQuantiles(Mat_gamma_BS,probs=0.975)
  gamma_sd=(var(Mat_gamma_BS))^(0.5)
  BS_gamma=matrix(rbind(c(gamma_sd,
                          gamma_quantile_025,
                          gamma_quantile_975)),ncol=1,byrow=T)
  
  log_g_quantile_025=colQuantiles(Mat_log_g_BS,probs=0.025)
  log_g_quantile_975=colQuantiles(Mat_log_g_BS,probs=0.975)
  log_g_sd=(colVars(Mat_log_g_BS))^(0.5)
  BS_log_g=matrix(rbind(c(log_g_sd,
                          log_g_quantile_025,
                          log_g_quantile_975)),ncol=n,byrow=T)
  
  Estimates=list(BS_gamma=BS_gamma,
                 BS_log_g=BS_log_g)
  
  return(Estimates)
  
}
