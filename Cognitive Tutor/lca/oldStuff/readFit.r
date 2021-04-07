library(tidyverse)
#source('logLik.r')

getNu <- function(ppp){
    nclass <- length(grep('meanTime',names(ppp)))
    if(nclass==2){
        nu <- ppp[grep('nu[',names(ppp),fixed=TRUE)]
        nu <- cbind(nu,1-nu)
    } else nu <- matrix(ppp[grep('nu[',names(ppp),fixed=TRUE)],ncol=nclass)
    nu
}

toMat <- function(vec){
    nnn <- names(vec)[length(vec)]
    dims <- strsplit(nnn,'\\[|\\,|\\]')
    dims <- na.omit(as.numeric(dims[[1]]))
    matrix(vec,nrow=dims[1])
}

parList <- function(ppp){
    nnn <- names(ppp)
    parNames <- strsplit(nnn,'[',fixed=TRUE)%>%map_chr(~.[1])
    lst <- split(ppp,parNames)

    for(pn in names(lst))
        if(grepl('\\,',names(lst[[pn]])[1]))
            lst[[pn]] <- toMat(lst[[pn]])
    lst
}

fits <- list()
#LL <- numeric(10)
for(i in 3:10){
    load(paste0('fit',i,'.RData'))
    fits[[as.character(i)]] <- fit
 #   print(system.time(LL[i] <- llik(sdat,fit$par)))
}

conv <- sapply(fits,function(x) x$return_code)

lp <- sapply(fits[conv==0],function(x) x$value)

load('LL.RData')

npar <- sapply((3:10)[conv==0],
               #function(i) i*4+(i-1)*(2+sdat$nstud+sdat$ncov+(i-2)/2)+sdat$nprob)
               function(i) (i-1)*sdat$nstud+i*4+sdat$nprob*3)

bic <- npar*log(sdat$nworked)-2*LL[3:7]
aic <- 2*(npar-LL[3:7])


### get hessian for 3-class model
ppp <- fits[['3']]$par

nn <- names(ppp)

getPar <- function(pname) ppp[startsWith(nn,pname)]


init <- sapply(
    c('meanTime','sigTime','effHint','effErr','probEff',
      'OmegaProb','sigProb','alpha','studEff','OmegaStud',
      'sigStud','beta'),
    getPar,
    simplify=FALSE)

nclass <- 3
## make things into matrices
init$probEff <- matrix(init$probEff,ncol=3)
init$studEff <- matrix(init$studEff,ncol=nclass-1)
init$OmegaStud <- matrix(init$OmegaStud,ncol=nclass-1)
init$OmegaProb <- matrix(init$OmegaProb,ncol=3)
init$beta <- matrix(init$beta,ncol=nclass-1)


mod <- rstan::stan_model('lca.stan')
sdat$nclass <- 3
newFit <- rstan::optimizing(mod,data=sdat,init=init,hessian=TRUE)
save(newFit,file='newFit3.RData')

                                        #invHes <- solve(-newFit$hessian)
## pick out fixed parameters from hessian (is that legal?)

nn <- colnames(newFit$hessian)
pars <- strsplit(nn,'.',fixed=TRUE)%>%map_chr(~.[1])
par.dims <- table(pars)

smallHess <- newFit$hessian[!pars%in%c('probEff','studEff'),!pars%in%c('probEff','studEff')]

invHes <- solve(-smallHess)

