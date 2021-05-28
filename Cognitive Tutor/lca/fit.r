library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load('smallSdat.RData')

mod <- stan('lca2class.stan',data=sdat,iter=2000,chains=8)
save(mod,file='lcaStan.RData')
