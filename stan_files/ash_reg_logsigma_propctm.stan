/* Stan Model Tutorial : Adaptive Shrinkage Model (log sigma uniform prior)
 * Author: Kushal Kumar Dey
 * -------------------------------------------------------
 * Date: 24th May, 2015
 */
  
  
  /* Model: Adaptive Shrinkage Method (ASH) - log sigma uniform prior, CTM prior for pi
   * Model description: b[j] ~ N(beta[j],s[j])
   * beta[j] ~ \sum_{k=1}^{K} p_k N(0, sigma[k])    for each j
   * sigma[k] ~ 1/sigma[k]
   * p_k ~ CTM_N(mu, Sigma)
   */
  
  
  data {
          int<lower=1> N;
          int<lower=1> K;
          vector[N] betahat;
          vector<lower=0>[N] se;

       }

parameters {
              vector<lower=0>[K] sigma; //... scales of individual mixing components
              vector[K] eta;  //... mixing proportions un-transformed
              vector[N] beta ; //... the true mean
              vector[K] mu_prop; //... the means for CTM prop
              vector<lower=0>[K] sigma_prop_vec; //... the scales for CTM prop
              corr_matrix[K] Omega_prop;  //... correlation matrix
           }

transformed parameters {
                            vector[K] log_sigma;
                            cov_matrix[K] Sigma_prop;
                            simplex[K] prop;
                            
                            log_sigma <- log(sigma);
                            Sigma_prop <- diag_matrix(sigma_prop_vec) * Omega_prop * diag_matrix(sigma_prop_vec);
                            prop <- softmax(eta);

                       }

model {
        real ps[K];
        mu_prop ~ normal(0,100); // vectorized, diffuse
        Omega_prop ~ lkj_corr(3.0); // regularize to unit correlation
        sigma_prop_vec ~ cauchy(0,100); // half-Cauchy due to constraint
        
        eta ~ multi_normal(mu_prop,Sigma_prop);
        
        for(n in 1:N)
        {
            //....  Mixture normal prior for beta 
            for(k in 1:K)
            {
                ps[k] <- log(prop[k])+ normal_log(beta[n],0,exp(log_sigma[k])); 
            }
            increment_log_prob(log_sum_exp(ps));
            betahat[n] ~ normal(beta[n],se[n]); //...  the data likelihood
         }
      }
