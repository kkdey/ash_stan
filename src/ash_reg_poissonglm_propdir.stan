/* Stan Model Tutorial : Adaptive Shrinkage for Poisson GLM Model (grid prior for beta)
 * Author: Kushal Kumar Dey
 * -------------------------------------------------------
 * Date: 26th May, 2015
 */
  
  
  /* Model: Adaptive Shrinkage Method (ASH) - Poisson GLM model
   * Model description: y[j] ~ Poi(mu+beta[t(j)],s[j])
   * mu ~ N(0,1000)
   * beta[t] ~ p_0 \delta_{0} + \sum_{k=1}^{K} p_k N(0, sigma[k])    for each j 
   * sigma[k] ~ m^{k} s_{min}  
   * m <- sqrt(2)
   * k chosen  to be the last k such that m^{k} s_{min} < s_{max}
   */
  
  
  data {
          int<lower=1> N; //...  the sample size 
          int<lower=1> K; //... the number of mixing components
          int<lower=0> B; //...  the number of batches 
          int<lower=0> counts[N]; //.. the count response  
          int<lower=1> batch[N] ;   //... the vector of tissue or batch labels
          vector<lower=0>[K] se_mix; //... standard error of mixing comp, obtained from grid method
          vector<lower=0>[K] alpha; //... the known prior parameter for mixing proportions
       }

parameters {
              simplex[K] prop; //... proportion vector for different mix comp.
              vector[B] beta ; //... the true parameter of interest
              real mu; //... the intercept term in Poisson GLM
           }



model {
          real ps[K-1];
          //.. int batch_lab;
          prop ~ dirichlet(alpha);
          mu ~ normal(0,1000);
  
          for(b in 1:B)
          {
              //....  Mixture normal prior for beta 
              //... ps[1] <- log(prop[1]);
    
              if (beta[b]==0) 
      
                  increment_log_prob(bernoulli_log(1,prop[1]));
    
              else
              {
      
                  for(k in 2:K)
                  {
                      ps[k-1] <- log(prop[k]) + normal_log(beta[b],0,se_mix[k]); 
                  }
      
              increment_log_prob(log_sum_exp(ps));
              
              }
      
          }
          
         
         
    
        for(n in 1:N)
       {
            //... batch_lab <- batch[n];
           counts[n] ~ poisson_log(mu+beta[batch[n]]) ; //...  the data likelihood
        }
    
        
      }
