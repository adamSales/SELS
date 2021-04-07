library(tidyverse)
library(rstan)


print(load('lcaDatMod.RData'))


fits <- list()
time <- NULL
lps <- list()
for(i in 2:10){
    load(paste0('lcaMLE',i,'.RData'))
    fits[[i-1]] <- fit
    time <- c(time,tdiff)
    lps[[i-1]] <- lp
}

converge <- map(fits,~.$return_code)

### look at different LP for each model



lp <- map_dbl(fits,~.$value)

ll <- map_dbl(fits,logLik,dat=sdat)

npar <- map_dbl(2:10,~(sdat$ncov+5+(.-1)/2)*.+5+sdat$ncov)

aic <- 2*npar-2*lp
bic <- npar*log(sdat$nworked)-2*lp
