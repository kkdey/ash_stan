
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
    out[k-1]=log((x[k])/(x[1]));
  }
  return(out)
}

##  *****************  simplex transformation of a vector  ********************************  ##

transform_simplex <- function(x)
{
  exp_x=exp(x);
  y= c( 1/(1+sum(exp_x)),exp_x/(1+sum(exp_x)));
  return(y);
  
}

## ******************************  random partition of a vector  ********************************************************** ##

chunk2 <- function(x,n) split(sample(x, length(x), replace=FALSE), ceiling(seq_along(x)/n)) 

##  ************************************  Fixing the Robbins Monro stepsize function **************************************  ##

stepsize <- function(t,kappa, tau) return (t+tau)^{-kappa};


##  ***************************  Stochastic MAP on each batch- worker function run parallel ********************************** ##

doichunk <- function(ichunk)
{
  stepsize <- function(t,kappa, tau) return ((t+tau)^{-kappa});
  
  subsample_iter = subsample_chunks[[ichunk]];
  #dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
  library(rstan)
  library(gtools)
  params_beta=params[1:length(beta_data)];
  params_prop_rev = params[-(1:length(beta_data))];
  
  params_sub_iter=c(params_beta[subsample_iter],params_prop_rev);
  
  # prop <- rdirichlet(1,dir_param);
  # prop_rev <- reverse_transform_simplex(prop);
  # params <- c(beta_data[subsample_iter],prop_rev);
  
  
  iter_data.ash <- list(N=batch_size,
                        K=length(se_mix_dummy),
                        betahat=beta_data[subsample_iter],
                        se=sebeta_data[subsample_iter],
                        se_mix=se_mix_dummy,
                        alpha=dir_param);
  
  iter_model <- stan(fit = model, data = iter_data.ash, chains = 0);
 ## iter_model <- stan(model_filename_MAP, data = iter_data.ash, chains=0)
  for (step in 1:nstep_descent)
  {
    grad <- grad_log_prob(iter_model,params_sub_iter) ;
    params_sub_iter <- params_sub_iter + 0.01*stepsize(iter, kappa, tau) * grad;
  }
  
  out_list <- list("model"=iter_model, "index"=subsample_iter, "params" = params_sub_iter, "score"=log_prob(iter_model, params_sub_iter) / batch_size);
  return(out_list);
}




##  **************************   The main function for Stochastic MAP estimation in ash  ****************************** ##

ash.MAP  <- function(beta_data, sebeta_data, model_filename_MAP, batch_size, niter=50, nstep_descent=1,  
                     multiplier=sqrt(2), model=NULL, exec.method = c("snow", "foreach", "serial"),kappa=0.5, tau=1)
{
  
  ## Creating the Stan fit model for a dummy variable set of batch_size length
  
  if (is.null(model))
  {
    dummy_subsample <- sample(1:length(beta_data),batch_size);
    se_mix_dummy=c(0,autoselect.mixsd(beta_data,sebeta_data,multiplier));
    dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
    dummy_data.ash <- list(N=batch_size,
                           K=length(se_mix_dummy),
                           betahat=beta_data[dummy_subsample],
                           se=sebeta_data[dummy_subsample],
                           se_mix=se_mix_dummy,
                           alpha=dir_param);
    model <- stan(model_filename_MAP, data = dummy_data.ash, chains=0)  ## compiling the STAN model without running any chain
  }
  
  n=length(beta_data);
  prop <- rdirichlet(1,dir_param);
  prop_rev <- reverse_transform_simplex(prop);
  params <- c(beta_data,prop_rev); ## prior specification of the parameters before iteration
  
 # num_params_batch <- get_num_upars(model);  ## get unconstrained parameters in the STAN model
  grad_scale <- floor(length(beta_data) / batch_size); ## grad scale is same as number of chunks
  
  if(length(beta_data) %% batch_size !=0) warning ( "Data overflow : batch size does not divide the original beta vector size")
  
  score <- numeric (niter); 
  
  for(iter in 1:niter)
  {
    subsample_chunks =chunk2(1:length(beta_data),batch_size); ## breaking the vector into chunks (each of batch size)
    subsample_indices=as.numeric(unlist(subsample_chunks));
  
    
    if(exec.method=="snow")
    {
      ##  send the chunks with additional data to the parallel clusters (cores on localhost) 
      
        cls <- makeSOCKcluster(mc <- getOption("cl.cores", 4)); ## make clusters using socket connection (default=4)
        
        ## Exporting all the variables required for the parallel worker function
        
        clusterExport(cls,"grad_scale",envir=environment());
        clusterExport(cls,"model",envir=environment());
        clusterExport(cls,"sebeta_data",envir=environment());
        clusterExport(cls,"beta_data",envir=environment());
        clusterExport(cls,"se_mix_dummy",envir=environment());
        clusterExport(cls,"stepsize",envir=environment());
        clusterExport(cls,"batch_size",envir=environment());
        clusterExport(cls,"nstep_descent",envir=environment());
        clusterExport(cls,"subsample_chunks",envir=environment());
        clusterExport(cls,"params",envir=environment());
        clusterExport(cls, "iter", envir=environment());
        clusterExport(cls, "tau", envir=environment());
        clusterExport(cls, "kappa", envir=environment());
        clusterExport(cls, "dir_param", envir=environment());
        
        
        ichunks <- 1:grad_scale; ## define the number of parallel working chunks you need 
      
        suppressWarnings(iter_out <- clusterApply(cls,ichunks,doichunk)); ## applying the worker function "doichunk" on each chunk in parallel
        stopCluster(cls);  ## stop the cluster 
      
        ## updating the parameters (beta and the reverse transformed mixture proportions)
        
        new_params_beta = unlist(lapply(iter_out, function(x) x$params[1:batch_size]))[order(subsample_indices)];
        prop_rev_list=lapply(iter_out, function(x) x$params[-(1:batch_size)])
        new_prop_rev = Reduce("+", prop_rev_list) / length(prop_rev_list);
        new_params = c(new_params_beta, new_prop_rev);
        total_score = sum(unlist(lapply(iter_out,function(x) x$score)));
        
    }
      
    if(exec.method=="foreach")
    {
        ##  send the chunks with additional data to the parallel clusters (cores on localhost) 
      
        cl<-makeCluster(mc <- getOption("cl.cores", 4))
        registerDoParallel(cl)
        
        ## automatic export of environmental variables to the clusters 
    
        suppressWarnings(iter_out <- foreach(i=1:grad_scale, .packages=c('gtools','rstan'))%dopar% doichunk(i)); ## parallel application of doichunk worker
        
        ## updating the parameters (beta and the reverse transformed mixture proportions)
        
        new_params_beta = unlist(lapply(iter_out, function(x) x$params[1:batch_size]))[order(subsample_indices)];
        prop_rev_list=lapply(iter_out, function(x) x$params[-(1:batch_size)])
        new_prop_rev = Reduce("+", prop_rev_list) / length(prop_rev_list);
        total_score=sum(unlist(lapply(iter_out,function(x) x$score)));
    }
  
   if(exec.method=="serial")
   {
    
     suppressWarnings(iter_out <- foreach(i=1:grad_scale, .packages=c('gtools','rstan'))%do% doichunk(i)); ## serial appliction of doichunk worker 
     
     ## updating the parameters (beta and the reverse transformed mixture proportions)
     
     new_params_beta = unlist(lapply(iter_out, function(x) x$params[1:batch_size]))[order(subsample_indices)];
     prop_rev_list=lapply(iter_out, function(x) x$params[-(1:batch_size)])
     new_prop_rev = Reduce("+", prop_rev_list) / length(prop_rev_list);
     total_score=sum(unlist(lapply(iter_out,function(x) x$score)))
    
   }
  
   params=c(new_params_beta,new_prop_rev);
   print(paste("end of iteration",iter))
  
  }
  
  params_beta_out = new_params_beta;
  params_prop_rev_out = new_prop_rev;
  
  out_main = list("beta_est"=params_beta_out, "mix_prop_est"= transform_simplex(params_prop_rev_out),  "score"=total_score);
  return(out_main);
  
}