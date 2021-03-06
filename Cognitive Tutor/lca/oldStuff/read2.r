library(tidyverse)
library(flextable)
library(officer)

load('fit2.RData')


nu <- getNu(fit$par)

nu <- nu[sdat$stud,]

colMeans(nu)

ppp <- parList(fit$par)

names(ppp$beta) <- colnames(sdat$X)

invHes <- solve(fit$hessian)

SEs <- seList(invHes)

### higher -> higher prob(class 1)
betaMat <- as.data.frame(cbind(Est=ppp$beta,SE=SEs$beta,T=ppp$beta/SEs$beta)%>%round(2))
betaMat$T <- paste0(betaMat$T,
                    ifelse(2*pnorm(-abs(betaMat$T))<0.001,'***',
                    ifelse(2*pnorm(-abs(betaMat$T))<0.01,'**',
                    ifelse(2*pnorm(-abs(betaMat$T))<0.05,'*',''))))
betaMat <- rownames_to_column(betaMat,'Covariate')
betaMat <- flextable(betaMat)
save_as_docx(betaMat,path='LCAbetas.docx')

classPars <- matrix(nrow=8,ncol=2)
classParNames <- c('effErr','effHint','meanTime','sigTime')
classPars[seq(1,7,2),]<- do.call("rbind",ppp[classParNames])
classPars[seq(2,8,2),]<- do.call("rbind",SEs[classParNames])
rownames(classPars) <- paste(rep(classParNames,each=2),rep(c('','SE'),4))
colnames(classPars) <- paste('Class',1:2)

classPars <- rbind(classPars,prior=colMeans(nu))

classPars <- sapply(c('effErr','effHint','meanTime','sigTime'),
                    function(x) paste0(round(ppp$effHint,2),' (',round(SEs$effHint,2),')')

### differences between classes
diff <- function(pname,ppp,invHes){
    diff <- ppp[[pname]][1]-ppp[[pname]][2]
    vdiff <- -sum(diag(invHes)[paste0(pname,'.',1:2)])+
        2*invHes[paste0(pname,'.1'),paste0(pname,'.2')]
    c(diff=diff,se=sqrt(vdiff))
}

diffs <- t(sapply(classParNames,diff,ppp=ppp,invHes=invHes))


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
