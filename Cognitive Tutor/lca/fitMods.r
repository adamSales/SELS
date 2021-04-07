

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

for(sig in c('0','.5','1','2','5')){
    fit <- get(paste0('fit',sig))
    assign(paste0('ppp',sig),parList(fit$par))
    assign(paste0('invHes',sig),solve(fit$hessian))
}
save(list=outer(c('ppp','invHes'),c('0','.5','1','2','5'),paste0),file='parSE.RData')

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
