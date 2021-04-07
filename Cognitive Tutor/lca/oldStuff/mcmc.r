library(rstan)

nclass <- 3
load(paste0('fit',nclass,'.RData')

ppp <- fit$par

nn <- names(ppp)

getPar <- function(pname) ppp[startsWith(nn,pname)]


init <- lapply(
    c('meanTime','sigTime','effHint','effErr','probEff',
      'OmegaProb','sigProb','alpha','studEff','OmegaStud',
      'sigStud','beta'),
    getPar)

## make things into matrices
init$probEff <- matrix(init$probEff,ncol=3)
init$studEff <- matrix(init$studEff,ncol=nclass-1)
