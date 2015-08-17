/* Stan Model Tutorial : (ash) Adaptive Shrinkage Model (Full Bayes approach)
 * Author: Kushal Kumar Dey
 * -------------------------------------------------------
 * Date: 16th August, 2015
 */
  
  
/* Model: Adaptive Shrinkage Method (ASH) - Full Bayes Approach
 * Model description: b[j] ~ N(beta[j],s[j])
 * beta[j] ~ p_0 \delta_{0} + \sum_{k=1}^{K} p_k N(0, sigma[k])    for each j
 * p_k ~ Dir(alpha);
 * sigma[k] ~ m^{k} s_{min}  
 * m <- sqrt(2)
 * k chosen  to be the last k such that m^{k} s_{min} < s_{max}
 */
 
 data {
         int<lower=1> N; //...  the sample size 
         int<lower=1> K; //... the number of mixing components 
         vector[N] betahat;   //... the betahat (effects) values vector
         vector<lower=0>[N] sebetahat; //... the standard error of the betahat (effects) values
         vector<lower=0>[K] se_mix; //... standard error of mixing comp, obtained from grid 
         vector<lower=0>[K] alpha; //... the known prior parameter for mixing proportions
      }
      
 parameters {
              simplex[K] prop; //... proportion vector for different mix comp.
            }
            
transformed parameters {
                          vector[K-1] prop_reverse;  // we reverse the simplex 
                          for(k in 2:K)
                          {
                            prop_reverse[k-1] <- log((prop[k])/(prop[1]));
                          }
                         
                       }
                       
model {
        real ps[K];
        //.. int batch_lab;
        prop ~ dirichlet(alpha);
        
        for(n in 1:N)
        {
           //...mixture normal for betahat given pi (mixture prop)
           
           for(k in 1:K)
           {
             ps[k] <- log(prop[k]) + normal_log(betahat[n],0,sqrt(pow(se_mix[k],2)+pow(sebetahat[n],2)));
           }
           
           increment_log_prob(log_sum_exp(ps));
        }
}


  

  