lologit <- function(x,eta,cc)
    log(
        ifelse(x==1,
               1-plogis(eta-cc[1]),
        ifelse(x==2, plogis(eta-cc[1])-plogis(eta-cc[2]),
               plogis(eta-cc[2])))
    )

llik1 <- function(hint,err, ltime,probEff,effHint,cHint,effErr,cErr,meanTime,sigTime,nu)
    log(
        exp(
            log(nu)+
            lologit(hint,probEff[1]+effHint[1],cHint)+
            lologit(err,probEff[2]+effErr[1],cErr)+
            dnorm(ltime,probEff[3]+meanTime[1],sigTime[1],log=TRUE)
        )+
        exp(
            log(1-nu)+
            lologit(hint,probEff[1]+effHint[2],cHint)+
            lologit(err,probEff[2]+effErr[2],cErr)+
            dnorm(ltime,probEff[3]+meanTime[2],sigTime[2],log=TRUE)
        )
    )

llik1oneClass <- function(hint,err, ltime,probEff,effHint,cHint,effErr,cErr,meanTime,sigTime,nu)
    lologit(hint,probEff[1]+effHint,cHint)+
        lologit(err,probEff[2]+effErr,cErr)+
        dnorm(ltime,probEff[3]+meanTime,sigTime,log=TRUE)



llikLoo <- function(data_i,draws){
    if(!is.matrix(draws)) draws <- rbind(draws)

    nnn <- colnames(draws)
    nu <- grep(paste0('nu[',data_i$stud),nnn,fixed=TRUE)
    if(!length(nu)) nu <- grep('^nu',nnn)

    lologit <- function(x,eta,cc)
        log(
            ifelse(x==1,
                   1-plogis(eta-cc[1]),
            ifelse(x==2, plogis(eta-cc[1])-plogis(eta-cc[2]),
                   plogis(eta-cc[2])))
        )

    llik1 <- function(hint,err, ltime,probEff,effHint,cHint,effErr,cErr,meanTime,sigTime,nu)
        log(
            exp(
                log(nu)+
                lologit(hint,probEff[1]+effHint[1],cHint)+
                lologit(err,probEff[2]+effErr[1],cErr)+
                dnorm(ltime,probEff[3]+meanTime[1],sigTime[1],log=TRUE)
            )+
            exp(
                log(1-nu)+
                lologit(hint,probEff[1]+effHint[2],cHint)+
                lologit(err,probEff[2]+effErr[2],cErr)+
                dnorm(ltime,probEff[3]+meanTime[2],sigTime[2],log=TRUE)
            )
        )


    apply(draws,1,function(x)
        llik1(
            hint=data_i$hint,
            err=data_i$err,
            ltime=data_i$ltime,
            probEff=x[grep(paste0('probEff[',data_i$prob),nnn,fixed=TRUE)],
            effHint=x[grep('effHint',nnn)],
            cHint=x[grep('cHint',nnn)],
            effErr=x[grep('effErr',nnn)],
            cErr=x[grep('cErr',nnn)],
            meanTime=x[grep('meanTime',nnn)],
            sigTime=x[grep('sigTime',nnn)],
            nu=x[nu]
        ))

}

llik <- function(ppp,sdat){
    pars <-
        cbind(
            hint=sdat$hint,
            err=sdat$err,
            ltime=sdat$ltime,
            probEff=ppp1$probEff[sdat$prob,]
        )

    if(length(ppp$effErr)==2){
        pars <- cbind(
            pars,
            nu=if(length(ppp$nu)==1){
                   rep(ppp$nu,sdat$nworked)
               } else ppp$nu[sdat$stud]
        )
        func <- llik1
    } else{
        pars <- cbind(
            pars,
            nu=rep(1,nrow(pars))
        )
        func <- llik1oneClass
    }



    lliks <- apply(pars,1,
                   function(x)
                       func(hint=x['hint'],
                            err=x['err'],
                            ltime=x['ltime'],
                            probEff=x[4:6],
                            effHint=ppp$effHint,
                            cHint=ppp$cHint,
                            effErr=ppp$effErr,
                            cErr=ppp$cErr,
                            meanTime=ppp$meanTime,
                            sigTime=ppp$sigTime,
                            nu=x['nu']
                            )
                   )
    sum(lliks)
}

