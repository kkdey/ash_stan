### Example script ash_cor

R <- matrix(runif(16), ncol=4) 
R <- (R * lower.tri(R)) + t(R * lower.tri(R)) 
diag(R) <- 1 
eigen(R)$val 
Q <- nearPD(R, posd.tol=1.e-04)$mat 
Q <- sweep(Q,diag(Q), MARGIN=1,"/")
det(Q)
cormat <- as.matrix(Q)
out <- ash_cor(cormat,nsamples=10)
det(out)
