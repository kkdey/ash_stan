
## reverse transformation of the simplex

reverse_transform_simplex <- function(x)
{
  K=length(x);
  out=numeric(K-1);
  for(k in 2:K)
  {
    out[k-1]=log((x[k]+1)/(x[1]+1));
  }
  return(out)
}

