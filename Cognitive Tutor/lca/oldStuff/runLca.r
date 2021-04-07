library(rstan)
load('lcaDatMod.RData')

model <- stan_model('lcaPoisCov.stan')

for(i in 3:10){
    print(i)
    print(Sys.time())
    sdat$nclass <- i
    fit <- optimizing(model,data=sdat,draws=10,importance=TRUE)
    save(fit,model,sdat,file=paste0('fit',i,'class.RData'))
    rm(fit)
}
