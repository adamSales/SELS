library(parallel)

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




arguments <- function(i,sdat,par,nu,probEff,nclass){
                                        #attach(sdat)
    if(i%%1000==0 & i<18000) print(round(i/18000*100))
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
    nclass <- length(grep('effHint',names(ppp)))
    nu <- matrix(ppp[grep('nu[',names(ppp),fixed=TRUE)],ncol=nclass)
    nu2 <- nu[sdat$stud,]
    probEff <- matrix(ppp[grep('probEff[',names(ppp),fixed=TRUE)],ncol=3)
    pe2 <- probEff[sdat$prob,]
    eee <- environment()
    clusterExport(cl,c('oneObs','arguments','sdat','ppp','nu2','pe2','nclass'), envir=eee)
    system.time(args <- #mclapply(#1:10000,
                    parLapply(cl,
                                  1:sdat$nworked,
                      arguments,sdat=sdat,par=ppp,nu=nu2,probEff=pe2,nclass=nclass))#,mc.cores=50)
    clusterExport(cl,'args', envir=eee)
    llik <- #sum(
        do.call("c",
    #    mclapply(
                parLapply(cl,
                1:sdat$nworked,
                  function(i){
                      do.call("oneObs",args[[i]])

                  })#,
#                  mc.cores=50)
        )#)
#    stopCluster(cl)
    llik
}

LL <- numeric(10)
for(i in 3:10){
    load(paste0('fit',i,'.RData'))
    LL[i] <- llik(sdat,fit$par)
}
save(LL,file='LL.RData')
