
###  Stochastic MAP code for ash data

stochastic.map.ash.stan <- function(model_filename, beta_data, sebeta_data, subsample, n, batch_size,
                           niter, stepsize, model = NULL)
{
   # compile the model (if needed) and initial parameters.
   # here, parameters are initialized to zero.
  
  
   if (is.null(model))
   {
     dummy_subsample <- sample(1:length(beta_data),batch_size);
     se_mix_dummy=c(0,autoselect.mixsd(beta_data[dummy_subsample],sebeta_data[dummy_subsample],sqrt(2)));
     dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
     dummy_data.ash <- list(N=batch_size,
                            K=length(se_mix_dummy),
                            betahat=beta_data[dummy_subsample],
                            se=sebeta_data[dummy_subsample],
                            se_mix=se_mix_dummy,
                            alpha=dir_param);
     model <- stan(model_filename, data = dummy_data.ash, chains=0)
   }
  
   n=length(beta_data);
   params <- rep(0,get_num_upars(model));
   grad_scale <- n / batch_size;
  
  
    # run stochastic gradient ascent by subsampling a batch, computing a
    # noisy gradient, and taking a step.  at each point, compute the
    # the per-data log probability as a score to track.
  
    score <- numeric(niter);
    
    
    for (t in 1:niter)
    {
      subsample_iter <- sample(1:length(beta_data),batch_size);
      dir_param=c(1,rep(1/(length(se_mix_dummy)-1),length(se_mix_dummy)-1));
      iter_data.ash <- list(N=batch_size,
                             K=length(se_mix_dummy),
                             betahat=beta_data[subsample_iter],
                             se=sebeta_data[subsample_iter],
                             se_mix=se_mix_dummy,
                             alpha=dir_param);
      iter_model <- stan(fit = model, data = iter_data.ash, chains = 0);
      params <- c(beta_data[subsample_iter],prop_rev);
      prop <- rdirichlet(1,dir_param);
      prop_rev <- reverse_transform_simplex(prop);
      
      for (step in 1:nstep_descent)
      {
        grad <- grad_scale * grad_log_prob(iter_model,params) ;
        params <- params + stepsize * grad
      }
      
      score[t] <- log_prob(iter_model, params) / batch_size
    }
  
    # return the model object (in case we want to run again)
    #        the params vector
    #        the score trajectory
  
   list(model = model, params = params, score = score)
 }
