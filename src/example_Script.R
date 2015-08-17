
# ash experiments - full bayes 

betahat1=c(0,1.03,0.05,0.67,2.06,3.11,2.22,1.19,1.43,1.2,0.03,10.54,8.66);
sebetahat1=c(0.01,1,0.2,0.15,0.44,0.62,0.31,0.29,0.01,0.11,0.50,10.11,10.43);

betahat <- cbind(0*betahat1, rnorm(length(betahat1),1,5), rnorm(length(betahat1),0,5));
sebetahat <- cbind(sebetahat1, rep(5,length(betahat1)), rep(5,length(betahat1)));

source('ash_fullbayes.R')

res_fullbayes_ash <- ash.stan(betahat, sebetahat, stan_iter = 1000)

res_standard_ash <- do.call(cbind, lapply(1:dim(betahat)[2],function(ind) ash(betahat[,ind],sebetahat[,ind])$PosteriorMean))
