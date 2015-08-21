
## contour plots for LASSO and ash (adaptive shrinkage)

betahat1=c(0,1.03,0.05,0.67,2.06,3.11,2.22,1.19,1.43,1.2,0.03,10.54,8.66);
sebetahat1=c(0.01,1,0.2,0.15,0.44,0.62,0.31,0.29,0.01,0.11,0.50,10.11,10.43);

source('ash_fullbayes.R')

se_mix <- autoselect.mixsd(betahat1, sebetahat1, sqrt(2))
pi <- rep(1/length(se_mix), length(se_mix));

#pi <- c(1, rep(0, length(se_mix)-1))

#pi <- rdirichlet(1,c(1,rep(1/(10*(length(se_mix)-1)), length(se_mix)-1)));

contour_function_ash <- function(beta1, beta2, se_mix, pi)
{
  out1 <- log(sum(pi*dnorm(beta1,0,se_mix)));
  out2 <- log(sum(pi*dnorm(beta2,0,se_mix)));
  out <- out1+out2;
  return(out)
}

beta1.vec <- seq(-1, 1, length.out = 1000);
beta2.vec <- seq(-1, 1, length.out = 1000);
mat <- matrix(0, nrow=length(beta1.vec), ncol=length(beta2.vec));
for(m in 1:length(beta1.vec))
{
  for(n in 1:length(beta2.vec))
  {
    mat[m,n] = contour_function_ash(beta1.vec[m], beta2.vec[n], se_mix, pi)
  }
}
contour(x=beta1.vec,y=beta2.vec,mat)

lambda=2
contour_function_lasso <- function(beta1, beta2, lambda)
{
  out1 <- log(ddoublex(beta1,0,1));
  out2 <- log(ddoublex(beta2,0,1));
  out <- out1+out2;
  return(out)
}

beta1.vec <- seq(-1, 1, length.out = 1000);
beta2.vec <- seq(-1, 1, length.out = 1000);
mat <- matrix(0, nrow=length(beta1.vec), ncol=length(beta2.vec));
for(m in 1:length(beta1.vec))
{
  for(n in 1:length(beta2.vec))
  {
    mat[m,n] = contour_function_lasso(beta1.vec[m], beta2.vec[n], lambda)
  }
}
contour(x=beta1.vec,y=beta2.vec,mat)
