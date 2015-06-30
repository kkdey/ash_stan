
setwd("Documents/ash_stan/src/");

# niter=300;
# nstep_descent=1;
# stepsize=0.01;
exec.maethod <- 'foreach';

source('ash_stochastic_map.R');
beta_data=c(rnorm(10,0,1),rnorm(10,1,1),rnorm(10,0,2),rnorm(10,0,10),rnorm(10,4,1));
sebeta_data=c(rep(1,10),rep(1,10),rep(2,10),rep(10,10),rep(1,10));

model_filename_MAP <- 'ash_reg_stan_propdir.stan';
batch_size = 50;
res1 <- ash.MAP (beta_data, sebeta_data, model_filename_MAP, batch_size, exec.method="snow")
res1$beta_est;


res2 <- ash(beta_data,sebeta_data);
res2$PosteriorMean;

source('ash_hmc_nuts.R')
res3 <- ash.HMC(beta_data, sebeta_data, model_filename_MAP)
res3$beta_est
