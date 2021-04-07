#### started writing function to compute log likelihood from my LCA model
#### cuz LP (what stan returns) isn't quite it, apparently
### (it's "up to a constant" whatever that means)
### but then realized that to get the log likelihood we have to integrate
### out the random effects. which I don't know how to do
### (I assume there's no analytic solution)

library(tidyverse)
library(rstan)
library(mvtnorm)

makeParams <- function(ppp){
 genPars <- map(
  c('meanTime','sigTime','effHint','effErr')%>%setNames(.,.),
  ~ppp[startsWith(names(ppp),.)]
 )

 nclass <- length(genPars$meanTime)

 genPars$SigmaProb <- matrix(ppp[startsWith(names(ppp),'SigmaProb')],3,3)
 genPars$SigmaStud <-
  matrix(ppp[startsWith(names(ppp),'SigmaStud')],nclass-1,nclass-1)

 probEff <- matrix(ppp[startsWith(names(ppp),'probEff')],ncol=3)
 nu <- matrix(ppp[startsWith(names(ppp),'nu[')],ncol=nclass)
 studEff <- matrix(ppp[startsWith(names(ppp),'studEff')],ncol=nclass-1)

 list(pars=genPars,probEff=probEff,nu=nu,studEff=studEff)
}

oneWPoneClass <- function(pars,cc,hint,err,time,nui,probEffp)
    log(nui[cc])+
        dpois(hint,exp(pars$effHint[cc]+probEffp[1]),log=TRUE)+
        dpois(err,exp(pars$effErr[cc]+probEffp[2]),log=TRUE)+
        dnorm(time,mean=pars$meanTime[cc]+probEffp[3],sd=pars$sigTime[cc],log=TRUE)

oneWP <- function(pars,nclass,hint,err,time,nui,probEffp)
    log(
        sum(
            exp(map_dbl(1:nclass,~oneWPoneClass(pars,.,hint,err,time,nui,probEffp)))))


logLik <- function(fit,dat){
    params <- makeParams(fit$par)

    nclass <- length(params$pars$meanTime)

    wp <- #sum(
        map_dbl(1:dat$nworked,
                function(i)
                    oneWP(pars=params$pars,
                          nclass=nclass,
                          hint=dat$hint[i],
                          err=dat$err[i],
                          time=dat$ltime[i],
                          nui=params$nu[dat$stud[i],],
                          probEffp=params$probEff[dat$prob[i],]
                          )
            )
        #)

    eta <- apply(params$studEff,1,dmvnorm,
                 sigma=params$pars$SigmaStud,log=TRUE)
    delta <- apply(params$probEff,1,dmvnorm,
                   sigma=params$pars$SigmaProb,log=TRUE)

    sum(wp[is.finite(wp)])+sum(eta)+sum(delta)
}



oneObs <- function(hint,err,ltime,probEff,effHint,effErr,meanTime,sigTime,nu){
### log(nu)+log(sum(dnorm*dpois*dpois))
    log(
        sum(
            exp(
                log(nu)+
                dpois(hint,lambda=exp(probEff[1]+effHint),log=TRUE)+
                dpois(err,lambda=exp(probEff[2]+effErr),log=TRUE)+
                dnorm(ltime,mean=probEff[3]+meanTime,sd=sigTime,log=TRUE)
            )
        )
    )
}

library(parallel)


arguments <- function(i,sdat,par,nu,probEff,nclass){
                                        #attach(sdat)
    if(i%%1000==0 & i<18000) print(round(i/18000*100))
    sdat$nclass <- nclass
    with(sdat,
         list(
             hint=hint[i],
             err=err[i],
             ltime=ltime[i],
             probEff=probEff[i,],
             effHint=par[paste0('effHint[',1:nclass,']')],
             effErr=par[paste0('effErr[',1:nclass,']')],
             meanTime=par[paste0('meanTime[',1:nclass,']')],
             sigTime=par[paste0('sigTime[',1:nclass,']')],
             nu=nu[i,]
         )
         )
    #detach(sdat)
}


llik <- function(sdat,ppp){
    cl <- makePSOCKcluster(7,outfile='out.txt')
    on.exit(    stopCluster(cl))
    nclass <- length(grep('effHint',names(ppp)))
    nu <- matrix(ppp[grep('nu[',names(ppp),fixed=TRUE)],ncol=nclass)
    nu2 <- nu[sdat$stud,]
    probEff <- matrix(ppp[grep('probEff[',names(ppp),fixed=TRUE)],ncol=3)
    pe2 <- probEff[sdat$prob,]
    eee <- environment()
    clusterExport(cl,c('oneObs','arguments','sdat','ppp','nu2','pe2','nclass'), envir=eee)
    system.time(args <- parLapply(cl,#1:10000,
                                  1:sdat$nworked,
                      arguments,sdat=sdat,par=ppp,nu=nu2,probEff=pe2,nclass=nclass))
    clusterExport(cl,'args', envir=eee)
    llik <- #sum(
        parSapply(cl,
                  1:sdat$nworked,
                  function(i){
                      do.call("oneObs",args[[i]])

                  }
                  )
        #)

    llik
}

