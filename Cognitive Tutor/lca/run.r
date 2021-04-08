library(loo)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(flextable)
library(officer)
library(tidyverse)
library(knitr)
library(kableExtra)
library(texreg)

source('llik.r')
source('extractFunctions.r')

pars <- c('meanTime','sigTime','effErr','effHint')

source('makeData.r')
load('smallSdat.RData')


### fit models
source('fitMods.r')

### model comparison table
source('modelComparison.r')

### main parameter & coef tables
source('mainTables.r')

### stuff for appendix 2
source('appendix2.r')
