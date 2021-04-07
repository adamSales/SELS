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

seList <- function(invHes){
    nnn <- rownames(invHes)
    parNames <- strsplit(nnn,'.',fixed=TRUE)
    split(sqrt(-diag(invHes)),map_chr(parNames,~.[1]))
}

### differences between classes
diff <- function(pname,ppp,invHes){
    diff <- ppp[[pname]][1]-ppp[[pname]][2]
    vdiff <- -sum(diag(invHes)[paste0(pname,'.',1:2)])+
        2*invHes[paste0(pname,'.1'),paste0(pname,'.2')]
    c(diff=diff,se=sqrt(vdiff))
}



ologit <- function(x,eta,cc)
    ifelse(x==1,
           1-plogis(eta-cc[1]),
    ifelse(x==2, plogis(eta-cc[1])-plogis(eta-cc[2]),
           plogis(eta-cc[2])))

posterior <- function(i,ppp,sdat){
    ## pr(state 1|X=x,Y=y)=pr(state 1|X=x)pr(Y=y|state 1)/Pr(Y=y)
    ## =exp(log(nu[stud[w]])+
    ## ordered_logistic_lpmf(hint[w]|probEff[prob[w]][1]+effHint[1],cHint)+
    ## ordered_logistic_lpmf(err[w]|probEff[prob[w]][2]+effErr[1],cErr)+
    ## normal_lpdf(ltime[w]| probEff[prob[w]][3]+meanTime[1],sigTime[1]))
    ## normalizing constant: Pr(Y=y|X=x)=Pr(Y=y|state 1)Pr(State 1|X)+Pr(Y=y|state 2)Pr(state 2|X)
    ppp$nu[sdat$stud[i]]*lik(i,1,ppp,sdat)/(ppp$nu[sdat$stud[i]]*lik(i,1,ppp,sdat)+(1-ppp$nu[sdat$stud[i]])*lik(i,2,ppp,sdat))
}


lik <- function(i,cls,ppp,sdat)
    with(sdat,
         ologit(hint[i],ppp$probEff[prob[i],1]+ppp$effHint[cls],ppp$cHint)*
         ologit(err[i],ppp$probEff[prob[i],2]+ppp$effErr[cls],ppp$cErr)*
         dnorm(ltime[i],ppp$probEff[prob[i],3]+ppp$meanTime[cls],ppp$sigTime[cls])
         )


postStud <- function(sss,sdat,ppp)
    sapply(which(sdat$stud==sss),posterior,ppp=ppp,sdat=sdat)


studPlots <- function(sss,sdat,ppp,probDat){
    opar <- par()
    par(mfrow=c(4,1))
    par(mar=c(3.1,4.1,4.1,2.1))
    on.exit(par(mfrow=opar$mfrow,mar=opar$mar))

    dotted <- c(c(1:3)*floor(sum(sdat$stud==sss)/4),sum(sdat$stud==sss))

    plot(sdat$ltime[sdat$stud==sss],type='l',ylab='log Time',
         main= paste0('Subject: ',probDat$field_id[sdat$stud==sss][1],'; Nobs: ',sum(sdat$stud==sss),'\nInteraction Time'),
         xaxt='n')
    axis(side=1, at=seq(0,80,20))
    abline(v=dotted,lty=2)

    eee <- as.table(probDat$nerrs1[sdat$stud==sss])
    names(eee) <- 1:length(eee)
    plot(eee,ylab='Count',xaxt='n',main='Number of Errors')
    axis(side=1, at=seq(0,80,20))
    abline(v=dotted,lty=2)

    hhh <- as.table(probDat$nhints1[sdat$stud==sss])
    names(hhh) <- 1:length(hhh)
    plot(hhh,ylab='Count',xaxt='n',main='Number of Hints')
    axis(side=1,at=seq(0,80,20))
    abline(v=dotted,lty=2)

    plot(postStud(sss,sdat,ppp),type='l',xaxt='n',main='Pr(State=2)',xlab='',ylab='')
    axis(side=1,at=seq(0,80,20))
    abline(v=dotted,lty=2)
}


diff <- function(pname,ppp,invHes){
    diff <- ppp[[pname]][1]-ppp[[pname]][2]
    vdiff <- -sum(diag(invHes)[paste0(pname,'.',1:2)])+
        2*invHes[paste0(pname,'.1'),paste0(pname,'.2')]
    c(diff=diff,se=sqrt(vdiff))
}

classParFun <- function(parname,ppp,SEs,invHes,digits=2){
    rnd <- function(x) sprintf(paste0("%.",digits,"f"), round(x,digits))

    row <- paste0(
        rnd(ppp[[parname]]),
        ' (',
        rnd(SEs[[parname]]),
        ')')
    pardiff <- diff(parname,ppp,invHes)
    p <- 2*pnorm(-abs(pardiff[1]/pardiff[2]))
    row <- c(row,
             paste0(
                 rnd(pardiff[1]),
                 ' (',
                 rnd(pardiff[2]),
                 ')',
                 ifelse(p<0.001,'***',
                 ifelse(p<0.01,'**',
                 ifelse(p<0.05,'*','')))
             )
             )
    tibble(
        Parameter=parname,
        `Class 1`=row[1],
        `Class 2`=row[2],
        `Difference`=row[3]
    )
}


altClassTab1 <- function(parname,ppp,SEs,invHes,sig){
    DIFF <- diff(parname,ppp,invHes)
    tibble(
        Parameter=parname,
        sigStud=sig,
        `Class 1`=ppp[[parname]][1],
        SE1=SEs[[parname]][1],
        `Class 2`=ppp[[parname]][2],
        SE2=SEs[[parname]][2],
        Difference=DIFF[1],
        SEdiff=DIFF[2]
    )

}

altCoefTab1 <- function(ppp,SEs,sig)
    tibble(
        covariate=names(ppp$beta),
        sigStud=sig,
        Est=ppp$beta,
        SE=SEs$beta
    )

altClassTab <- function(parList,ppp,SEs,invHes,sig){

    out <- bind_rows(
        map_dfr(parList,altClassTab1,ppp=ppp,SEs=SEs,invHes=invHes,sig=sig),
        tibble(
            Parameter='Probability',
            sigStud=sig,
            `Class 1`=mean(ppp$nu),
            `Class 2`=1-mean(ppp$nu)
        )
    )
    if(out$`Class 1`[out$Parameter=='meanTime']>out$`Class 2`[out$Parameter=='meanTime']){
        nnn <- names(out)
        nnn[grep('1',names(out))] <- names(out)[grep('2',names(out))]
        nnn[grep('2',nnn)] <- names(out)[grep('1',names(out))]
        names(out) <- nnn
        out$Difference <- -out$Difference
    }
    out
}
