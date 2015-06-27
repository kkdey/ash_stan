
###   ****************************************************************************************************************************  ####
###   ********************************** A Stochastic MAP framework for inference in ash  ****************************************  ####
###   ****************************************************************************************************************************  ####

library(gtools)
library(rstan)
library(ashr)
library(snow)
library(parallel)
library(doParallel)
library(foreach)

## ****************** Select the mixture components and corresponding standard errors of mixing components ********************** ##

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



## ***************  reverse transformation of a simplex vector   *****************##

reverse_transform_simplex <- function(x)
{
  K=length(x);
  out=numeric(K-1);
  for(k in 2:K)
  {
    out[k-1]=log((x[k]+1)/(x[1]+1));
  }
  return(out)
}

##  **************************   The main function for Stochastic MAP estimation in ash  ****************************** ##

ash.MAP  <- function(beta_data, sebeta_data, model_filename, batch_size, niter, nstep_descent, stepsize, 
                     multiplier=sqrt(2), model=NULL, exec.method = c("snow", "foreach", "serial"))
{
  
  ## Creating the Stan fit model for a dummy variable set of batch_size length
  
  if (is.null(model))
  {
    dummy_subsample <- sample(1:length(beta_data),batch_size);
    se_mix_dummy=c(0,autoselect.mixsd(beta_data[dummy_subsample],sebeta_data[dummy_subsample],multiplier));
    dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
    dummy_data.ash <- list(N=batch_size,
                           K=length(se_mix_dummy),
                           betahat=beta_data[dummy_subsample],
                           se=sebeta_data[dummy_subsample],
                           se_mix=se_mix_dummy,
                           alpha=dir_param);
    model <- stan(model_filename, data = dummy_data.ash, chains=0)  ## compiling the STAN model without running any chain
  }
  
  n=length(beta_data);
  params <- rep(0,get_num_upars(model));  ## get unconstrained parameters in the STAN model
  grad_scale <- n / batch_size;
  
  score <- numeric (niter); 
  
  doichunk <- function(ichunk)
  {
    
    reverse_transform_simplex <- function(x)
    {
      K=length(x);
      out=numeric(K-1);
      for(k in 2:K)
      {
        out[k-1]=log((x[k]+1)/(x[1]+1));
      }
      return(out)
    }
    
    subsample_iter <- sample(1:length(beta_data),batch_size);
    dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
    iter_data.ash <- list(N=batch_size,
                          K=length(se_mix_dummy),
                          betahat=beta_data[subsample_iter],
                          se=sebeta_data[subsample_iter],
                          se_mix=se_mix_dummy,
                          alpha=dir_param);
    library(rstan)
    library(gtools)
    iter_model <- stan(fit = model, data = iter_data.ash, chains = 0);
    prop <- rdirichlet(1,dir_param);
    prop_rev <- reverse_transform_simplex(prop);
    params <- c(beta_data[subsample_iter],prop_rev);
    
    for (step in 1:nstep_descent)
    {
      grad <- grad_scale * grad_log_prob(iter_model,params) ;
      params <- params + stepsize * grad
    }
    
    out_list <- list("model"=model, "index"=subsample_iter, "params" = params, "score"=log_prob(iter_model, params) / batch_size);
    return(out_list);
  }
  
  
  
  if(exec.method=="snow")
  {
      cls <- makeSOCKcluster(mc <- getOption("cl.cores", 8));
      clusterExport(cls,"grad_scale",envir=environment());
      clusterExport(cls,"model",envir=environment());
      clusterExport(cls,"sebeta_data",envir=environment());
      clusterExport(cls,"beta_data",envir=environment());
      clusterExport(cls,"se_mix_dummy",envir=environment());
      clusterExport(cls,"stepsize",envir=environment());
      clusterExport(cls,"batch_size",envir=environment());
      clusterExport(cls,"nstep_descent",envir=environment());
        
      ichunks <- 1:niter;
      iter_out <- clusterApply(cls,ichunks,doichunk);
      stopCluster(cls)
        
  }
      
  if(exec.method=="foreach")
  {
    cl<-makeCluster(mc <- getOption("cl.cores", 8))
    registerDoParallel(cl)
    
    iter_out <- foreach(i=1:niter, .packages=c('gtools','rstan'))%dopar%
    {
       doichunk(i);
    }
    
  }
  
  if(exec.method=="serial")
  {
    
    iter_out <- foreach(i=1:niter, .packages=c('gtools','rstan'))%do%
    {
      doichunk(i);
    }
    
  }
  
  
  
}