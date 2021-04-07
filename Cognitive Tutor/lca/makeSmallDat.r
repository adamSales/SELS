load('smallSdat.RData')

### only first 500 students, for experimenting
sdat <- within(sdat,
               {
                   prob <- prob[stud<501]
                   err <- err[stud<501]
                   hint <- hint[stud<501]
                   ltime <- ltime[stud<501]
                   stud <- stud[stud<501]
                   nstud <- 500
                   nworked <- length(prob)
                   X <- X[1:500,]
               }
               )

mod <- stan_model('lca2classStudEff2.stan')

sdat$studEffsd <- -1
fit0 <- optimizing(mod,data=sdat)
ppp0 <- parList(fit0$par)


sdat$studEffsd <- 0.7
fit7 <- optimizing(mod,data=sdat)

modOrig <- stan_model('lca2classUnConst.stan')

fitOrig <- optimizing(modOrig,data=sdat,init=list(sigStud=2,sigTime=runif(2,0.5,2)))

pppOrig <- parList(fitOrig$par)


posterior <- function(i,ppp,sdat){
    if(is.null(ppp$probEff)) lik <- lik2
    ## pr(state 1|X=x,Y=y)=pr(state 1|X=x)pr(Y=y|state 1)/Pr(Y=y)
    ## =exp(log(nu[stud[w]])+
    ## ordered_logistic_lpmf(hint[w]|probEff[prob[w]][1]+effHint[1],cHint)+
    ## ordered_logistic_lpmf(err[w]|probEff[prob[w]][2]+effErr[1],cErr)+
    ## normal_lpdf(ltime[w]| probEff[prob[w]][3]+meanTime[1],sigTime[1]))
    ## normalizing constant: Pr(Y=y|X=x)=Pr(Y=y|state 1)Pr(State 1|X)+Pr(Y=y|state 2)Pr(state 2|X)
    ppp$nu[sdat$stud[i]]*lik(i,1,ppp,sdat)/(ppp$nu[sdat$stud[i]]*lik(i,1,ppp,sdat)+(1-ppp$nu[sdat$stud[i]])*lik(i,2,ppp,sdat))
}


lik2 <- function(i,cls,ppp,sdat)
    with(sdat,
         ologit(hint[i],ppp$effHint[cls],ppp$cHint)*
         ologit(err[i],ppp$effErr[cls],ppp$cErr)*
         dnorm(ltime[i],ppp$meanTime[cls],ppp$sigTime[cls])
         )
