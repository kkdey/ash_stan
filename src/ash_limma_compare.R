
# compare between ash, ash.stan and limma 

# In this script, we compare between limma and ash/ ash.stan and show why we feel one should
# use ash prior instead of limma to identify the differentially expressed genes.

library(limma)
#################################################################################

e=read.csv("http://rafalab.jhsph.edu/688/expression.csv",row.names=1)
e=log2(e) # log2 transform

m1=rowMeans(e[,1:3]) # group 1 
m2=rowMeans(e[,4:6]) # group 2
a=rowMeans(e)

which.min(m2-m1)
grp = rep(1:2, each=3)

o=order(m2-m1)
rownames(e)[o[1:10]] 
#source("http://bioconductor.org/biocLite.R")
#biocLite("genefilter")
library(genefilter)
library(limma)
library(ashr)


s1=rowSds(e[,1:3])
s2=rowSds(e[,4:6])
ttest=(m2-m1)/sqrt(s1^2/3 + s2^2/3)

pval=2*(1-pt(abs(ttest),4))

smoothScatter(m2-m1, -log10(pval), main="volcano",xlab="Effect size",ylab="-log10(pvalue)",nbin=128)
plot(m2-m1, -log10(pval), main="volcano",xlab="Effect size",ylab="-log10(pvalue)",pch=19,cex=0.2,lwd=1)

tt = rowttests(as.matrix(e),factor(grp))

design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
fit <- lmFit(as.matrix(e),design);
ebayes_fit <- ebayes(fit);
ebayes_fit.pvalue <- ebayes_fit$p.value[,2];

smoothScatter(m2-m1, -log10(ebayes_fit.pvalue), main="volcano",xlab="Effect size",ylab="-log10(pvalue)",nbin=128)
plot(m2-m1, -log10(ebayes_fit.pvalue), main="volcano",xlab="Effect size",ylab="-log10(pvalue)",pch=19,cex=0.2,lwd=1)

e <- as.matrix(e);
betahat <- array(0,dim(e)[1]);
sebetahat <- array(0, dim(e)[1]);
for(g in 1:dim(e)[1])
{
  out <- summary(lm(e[g,]~as.factor(grp)));
  betahat[g] <- as.numeric(out$coefficients[2,1]);
  sebetahat[g] <- as.numeric(out$coefficients[2,2])
}

res <- ash(betahat,sebetahat,mixcompdist="normal",control=list(maxiter=2000))
smoothScatter(m2-m1, -log10(res$qvalue), main="volcano",xlab="Effect size",ylab="-log10(qvalue)",nbin=128)
plot(m2-m1, -log10(res$qvalue), main="volcano",xlab="Effect size",ylab="-log10(qvalue)",pch=19,cex=0.2,lwd=1)


