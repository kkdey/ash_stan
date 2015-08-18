
## ash full bayes STAN approach
library(rstan)
library(ashr)
library(foreach)
library(parallel)

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



ash.stan <- function(betahat, sebetahat, mult=sqrt(2), stan_iter, cores=detectCores(), shrink_factor=10)
{
  betahat <- as.matrix(betahat);
  sebetahat <- as.matrix(sebetahat);
  
  num_cases <- dim(as.matrix(betahat))[2];
  
  ind=1;
  
  se_mix=c(0,autoselect.mixsd(betahat[,ind],sebetahat[,ind],mult));
  dir_param=c(1,rep(1/length(se_mix),length(se_mix)-1));
  test.ash <- list(N=length(betahat[,ind]),
                   K=length(se_mix),
                   betahat=betahat[,ind],
                   sebetahat=sebetahat[,ind],
                   se_mix=se_mix,
                   alpha=dir_param);
  
  ### Fitting STAN
  fit1 <- stan(file="ash_fullbayes.stan",data=test.ash,chains=0);
  
  local_chunk <- function(fit, betahat, sebetahat, ind)
  {
    se_mix=c(0,autoselect.mixsd(betahat[,ind],sebetahat[,ind],mult));
    dir_param=c(1,rep(1/(shrink_factor*(length(se_mix)-1)),length(se_mix)-1));
    test.ash <- list(N=length(betahat[,ind]),
                     K=length(se_mix),
                     betahat=betahat[,ind],
                     sebetahat=sebetahat[,ind],
                     se_mix=se_mix,
                     alpha=dir_param);
    out <- stan(fit = fit, data = test.ash, chains = 1, iter=stan_iter);
    temp=summary(out)$summary[,"mean"];
    prop_mix_vec <- as.numeric(temp[1:length(se_mix)]); 
    
    #  beta_update_mat <- matrix(0, ncol=0, nrow=length(betahat[,ind]));
    
    beta_update_vec <- array(0,length(betahat[,ind]));
    scaling <- do.call(cbind,lapply(1:length(se_mix), function(k) se_mix[k]^2 /(se_mix[k]^2+sebetahat[,ind]^2)));
    omega_mix <- do.call(rbind,lapply(1:length(betahat[,ind]), function(num) prop_mix_vec*dnorm(betahat[num,ind],0,sd=sqrt(sebetahat[num,ind]^2 + se_mix^2))));
    omega_mix <- sweep(omega_mix,1,rowSums(omega_mix),'/');
    beta_update_vec <- rowSums(scaling*omega_mix)*betahat[,ind];
    
#    Serial execution of the above codes (beta_update)
    
#    for(num in 1:length(betahat[,ind]))
#    {
     # scaling <- se_mix^2 /(se_mix^2+sebetahat[num,ind]^2);
#      omega_mix <- prop_mix_vec*dnorm(betahat[num,ind],0,sd=sqrt(sebetahat[num,ind]^2 + se_mix^2));
#      omega_mix <- omega_mix/sum(omega_mix);
#      beta_update_vec[num] <- sum(scaling*omega_mix)*betahat[num,ind];
#    }
    
    out <- list("pi_local"=prop_mix_vec,"beta_local"=beta_update_vec)
    return(out)
  }

  out <- mclapply(1:num_cases, function(ind) local_chunk(fit1, betahat, sebetahat, ind), mc.cores=cores);
  beta_mat <- do.call(cbind, lapply(1:num_cases, function(x) out[[x]]$beta_local))
  prop_mix <- lapply(1:num_cases, function(x) out[[x]]$pi_local)
    
  out <- list("beta_update"=beta_mat, "pi_update"=prop_mix)
  return(out)
}