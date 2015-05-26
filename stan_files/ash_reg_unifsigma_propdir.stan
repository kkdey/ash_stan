
/* Stan Model Tutorial : Emp Bayes approach for Adaptive Shrinkage Model  
 * Author: Kushal Kumar Dey
 * -------------------------------------------------------
 * Date: 24th May, 2015
 */
  
  
  /* Model: Adaptive Shrinkage Method (ASH) - Emp Bayes Approach
   * Model description: b[j] ~ N(beta[j],s[j])
   * beta[j] ~ \sum_{k=1}^{K} p_k N(0, sigma[k])    for each j
   * sigma[k] ~ Uniform(0, infty)
   */
  
  
  data {
          int<lower=1> N;
          int<lower=1> K;
          vector[N] betahat;
          vector<lower=0>[N] se;
       }

parameters {
              vector<lower=0>[K] sigma;
              simplex[K] prop;
              vector[N] beta ;
           }

model {
          real ps[K];
          sigma ~ uniform(0,1000);  //.... prior for the variation

          for(n in 1:N)
          {
              //....  Mixture normal prior for beta 
              for(k in 1:K)
              {
                  ps[k] <- log(prop[k])+ normal_log(beta[n],0,sigma[k]); 
              }
              increment_log_prob(log_sum_exp(ps));
              betahat[n] ~ normal(beta[n],se[n]); //...  the data likelihood
          }
      }

