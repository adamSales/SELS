load('mod.RData')
load('fit1class.RData')
load('prelimMods.RData')

ppp1 <- parList(fit1)

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
