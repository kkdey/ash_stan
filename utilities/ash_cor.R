
library(ashr)
library(reshape)
library(reshape2)
library(Matrix)

## The function for correlation matrix ash framework

ash_cor <- function(cormat, nsamples)
{
  cor_table <- melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];
  
  
  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag[,3])/(1-cor_table_non_diag[,3]))
  cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
  options(warn=-1)
  fit=ash(cor_transform_mean_vec,cor_transform_sd_vec,mixcompdist="normal");
  ash_cor_vec=(exp(2*fit$PosteriorMean)-1)/(exp(2*fit$PosteriorMean)+1);
  
  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  new_mat <- dcast(newdata.table, X1~X2, value.var = "value")[,-1];
  new_mat[is.na(new_mat)]=1;
  pd_completion <- nearPD(as.matrix(new_mat), conv.tol=1e-06);
  new_mat <- sweep(pd_completion$mat,diag(pd_completion$mat), MARGIN=1,"/")
  return(new_mat)
}







