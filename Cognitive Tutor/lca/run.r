library(loo)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(flextable)
library(officer)
library(tidyverse)

 source('llik.r')
source('extractFunctions.r')

load('smallSdat.RData')

#mod <- stan('lca2class2.stan',data=sdat,iter=1000,sample_file='samp.txt',refresh=10)

mod <- stan_model('lca2class2.stan')

sdat$sigStud <- 0.001
fit0 <- optimizing(mod,data=sdat,hessian=TRUE)#,draws=200,importance_resampling=TRUE)

sdat$sigStud <- 0.5
fit.5 <- optimizing(mod,data=sdat,hessian=TRUE)#,draws=200,importance_resampling=TRUE)

sdat$sigStud <- 1
fit1 <- optimizing(mod,data=sdat,hessian=TRUE)#,draws=200,importance_resampling=TRUE)

sdat$sigStud <- 2
fit2 <- optimizing(mod,data=sdat,hessian=TRUE)#,draws=2000,importance_resampling=TRUE)

sdat$sigStud <- 5
fit5 <- optimizing(mod,data=sdat,hessian=TRUE)#,draws=2000,importance_resampling=TRUE)

save(fit0,fit.5,fit1,fit2,fit5,sdat,file='mod.RData')

for(sig in c('0','.5','2','5')){
    fit <- get(paste0('fit',sig))
    assign(paste0('ppp',sig),parList(fit$par))
    assign(paste0('invHes',sig),solve(fit$hessian))
}

save(list=outer(c('ppp','invHess'),c('0','.5','2','5'),paste0),file='parSE.RData')

oneClassMod <- stan_model('lca1class.stan')
fit1class <- optimizing(oneClassMod,data=sdat,hessian=TRUE,draws=2000,importance_resampling=TRUE)
save(fit1class,sdat,oneClassMod,file='fit1class.RData')

#### no student effects or covariates
mod0 <- stan_model('lca2class0.stan')
fit00 <- optimizing(mod0,data=sdat,hessian=TRUE,draws=2000,importance_resampling=TRUE)

#### yes student effects no covariates
mod1 <- stan_model('lca2class1.stan')
sdat$sigStud <- 1
fit01 <- optimizing(mod1,data=sdat,hessian=TRUE,draws=2000,importance_resampling=TRUE)

save(fit00,fit01,mod0,mod1,sdat,file='prelimMods.RData')


logLiks <- c(
    oneClass=llik(parList(fit1class$par),sdat),
    twoClass=llik(parList(fit00$par),sdat),
    twoClassStud=llik(parList(fit01$par),sdat),
    twoClassCovs=llik(ppp1,sdat)
    )

npar <- c(
    oneClass=
        4+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4, #ordered logit intercepts for hints & errors
    twoClass=
          4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        1, # alpha = class prob
    twoClassStud=
        4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        sdat$nstud, #student class probs
    twoClassCovs=
        4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        sdat$nstud+ #student class probs
        sdat$ncov # beta
)

npar2 <- c(
    oneClass=
        4+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4, #ordered logit intercepts for hints & errors
    twoClass=
          4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        1, # alpha = class prob
    twoClassStud=
        4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        2, # mean and var of student effects
    twoClassCovs=
        4*2+ # parameters that vary by class (x1)
        3*sdat$nprob+ #problem intercepts for each indicator
        4+ #ordered logit intercepts for hints & errors
        2+ # mean and var of student effects
        sdat$ncov # beta
)


bics <- npar*log(sdat$nworked)-2*logLiks
aics <- 2*(npar-logLiks)
abics <- npar*log((sdat$nworked+2)/24)-2*logLiks
caics <- npar*log(sdat$nworked+1)-2*logLiks

fitTab <- tibble(
    M=c(1,2,'',''),
    Condition=c('','HomStud','StudEff','StudCovs'),
    p=npar,
    LL=logLiks,
    AIC=aics,
    CAIC=caics,
    BIC=bics,
    ABIC=abics
)

save_as_docx(flextable(fitTab),path='../../../Manuscript/lcaTabs/lcaFit.docx')

bics2 <- npar2*log(sdat$nworked)-2*logLiks
aics2 <- 2*(npar2-logLiks)

############# parameters for full model
parList <- c('meanTime','sigTime','effErr','effHint')
SEs1 <- seList(invHes1)

names(ppp1$beta) <- colnames(sdat$X)

as.data.frame(cbind(Est=ppp1$beta,SE=SEs1$beta,pval=2*pnorm(-abs(ppp1$beta/SEs1$beta)))%>%round(3))%>%
    rownames_to_column('Covariate')%>%
    flextable()%>%
    save_as_docx(path='../../../Manuscript/lcaTabs/LCAbetas.docx')

map_dfr(
    parList,
    classParFun,
    ppp=ppp1,SEs=SEs1,invHes=invHes1,digits=3
)%>%
    bind_rows(tibble(Parameter="Probability",
                     `Class 1`=sprintf("%.3f", round(mean(ppp1$nu),3)),
                     `Class 2`=sprintf("%.3f", round(1-mean(ppp1$nu),3))))%>%
    flextable()%>%
    save_as_docx(path='../../../Manuscript/lcaTabs/LCApars.docx')



################### parameters for other sig_stud

altTab1 <- altClassTab(parList=parList,
                       ppp=ppp1,SEs=SEs1,invHes=invHes1,sig='1')

altTab2 <- altCoefTab1(ppp=ppp1,SEs=SEs1,sig='1')

for(sig in c('0','.5','2','5')){
    ppp <- get(paste0('ppp',sig))
    invHes <- get(paste0('invHes',sig))
    SEs <- seList(invHes)
    names(ppp$beta) <- colnames(sdat$X)

    altTab1 <- bind_rows(
        altTab1,
        altClassTab(parList=parList,ppp=ppp,SEs=SEs,invHes=invHes,sig=sig)
    )

    altTab2 <- bind_rows(
        altTab2,
        altCoefTab1(ppp=ppp,SEs=SEs,sig=sig)
    )
}

altTab1 <- arrange(
    altTab1,
    factor(Parameter,levels=c(parList,'Probability')),sigStud)

altTab2 <- arrange(
    altTab2,
    factor(covariate,levels=colnames(sdat$X)),sigStud)


altTab1%>%
    pivot_longer(`Class 1`:SEdiff,names_to='par',values_to='value')%>%
    mutate(
        which=ifelse(endsWith(par,'1'),'1',
              ifelse(endsWith(par,'2'),'2','Diff')),
        estSE=ifelse(startsWith(par,'SE'),'SE','Est'),
        sigStud=ifelse(sigStud=='.5','0.5',sigStud)
    )%>%
    select(-par)%>%
    pivot_wider(names_from='estSE',values_from='value')%>%
    ggplot(aes(
        x=factor(sigStud,
                 levels=as.character(sort(as.numeric(unique(sigStud))))),
        y=Est,ymin=Est-2*SE,ymax=Est+2*SE))+
    geom_point()+
    geom_errorbar()+
    facet_grid(Parameter~which)


as.data.frame(cbind(Est=ppp$beta,SE=SEs$beta,pval=2*pnorm(-abs(ppp$beta/SEs$beta)))%>%round(3))%>%
        rownames_to_column('Covariate')%>%
        flextable()%>%
        save_as_docx(path=paste0('../../../Manuscript/lcaTabs/LCAbetas',sig,'.docx'))
    map_dfr(
        c('effErr','effHint','meanTime','sigTime'),
        classParFun,
        ppp=ppp,SEs=SEs,invHes=invHes,digits=3
    )%>%
        bind_rows(tibble(Parameter="Probability",
                         `Class 1`=sprintf("%.3f", round(mean(ppp$nu),3)),
                         `Class 2`=sprintf("%.3f", round(1-mean(ppp$nu),3))))%>%
        flextable()%>%
        save_as_docx(path=paste0('../../../Manuscript/lcaTabs/LCApars',sig,'.docx'))
}

################### WAIC
draws <- fit1$theta_tilde

looDat <- as.data.frame(sdat[sapply(sdat,length)==sdat$nworked])

loo10 <-   loo_approximate_posterior(x=llikLoo,
                       draws=draws,
                       data=looDat,
                       log_p=fit1$log_p,
                       log_g=fit1$log_g,
                       cores=parallel::detectCores()
                       )

save(loo10,file='loosubsample.RData')
