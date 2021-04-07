library(rstan)
library(tidyverse)

small <- TRUE
source('lca.r')

sink('temp.out')
mod <- stan_model('lcaPoisCovUnConst.stan')
sink()


for(nclass in 3:10){
    print(nclass)
    print(Sys.time())
    sdat$nclass <- nclass
    fit <- optimizing(mod,data=sdat,iter=1e8,init=list(sigTime=rep(.5,nclass),sigProb=rep(.5,3),sigStud=rep(.5,nclass-1)),verbose=TRUE,draws=1,importance_resampling=TRUE)
    save(fit,sdat,mod,file=paste0('fit',nclass,'.RData'))
}
