library(rstan)

oneClassMod <- stan_model('lca1class.stan')

load('smallSdat.RData')

fit <- optimizing(oneClassMod,data=sdat,hessian=TRUE,
                  init=list(sigTime=runif(1,0.05,0.95),
                            sigProb=runif(3,0.05,0.95)
                            ),
                  draws=100)
