

##  ***************************   ash Full Bayes HMC + NUTS  approach  ******************************************** ##

autoselect.mixsd = function(betahat,sebetahat,mult){
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}



ash.HMC <- function(beta_data, sebeta_data, model_filename, multiplier=sqrt(2), nchains=1, niter=1000)
{
  
  source("http://mc-stan.org/rstan/stan.R") ## nchains will run in parallel cores
  
  se_mix=c(0,autoselect.mixsd(beta_data,sebeta_data,multiplier));
  dir_param=c(1,rep(1/(length(se_mix)-1),length(se_mix)-1));
  
  
  test_data.ash <- list(N=length(beta_data),
                        K=length(se_mix),
                        betahat=beta_data,
                        se=sebeta_data,
                        se_mix=se_mix,
                        alpha=dir_param);
  
  ### Fitting STAN
  
  fit1 <- suppressMessages(suppressWarnings(stan(file=model_filename, data=test_data.ash,iter=niter,chains=nchains,verbose=FALSE)));
  
  temp=summary(fit1)$summary[,"mean"]
  beta_temp=as.numeric(temp[grep("beta",names(temp))]);
  out_list = list("model"=fit1, "beta_est"=beta_temp);
  return(out_list)
  
}