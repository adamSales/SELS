library(tidyverse)

if(!exists('small')) small <- TRUE

if(small) cat('\n WARNING: working with 1-day data\n') else cat('\n Working with full data\n')

if(small){
    dfs <- load('../CTA_OneUnit_SDay.RData')
    probDat <- bind_rows(data_prob_es01_sday)
    rm(list=dfs)
} else load('../ctData.RData')

#pretest, post test, race, sex, free-lunch, ESL

probDat <- probDat%>%
#    filter(nerrs1<103,nhints1<50)%>%
    mutate(studID=as.numeric(as.factor(field_id)))%>%
    mutate(across(c(nhints1,nerrs1),~cut(.,c(-1,0,1,Inf),labels=FALSE)))

X <- probDat%>%
    group_by(studID)%>%
    summarize(across(c(pretest,gainscore,race,sex,frl),~.[1]))%>%
    mutate(race=as.factor(race))%>%
    arrange(studID)%>%
    mutate(across(c(pretest,gainscore),scale))

X <- model.matrix(studID~.,data=X)[,-1]




sdat <- with(probDat,
             list(
               prob=as.numeric(as.factor(kcs)),
               err=nerrs1,#as.numeric(err>0),
               hint=nhints1,#as.numeric(hint>0),
               ltime=scale(ltime)[,1],
               stud=studID,
               nprob=n_distinct(kcs),
               nstud=n_distinct(field_id),
               nworked=length(kcs),
               zeros=c(0,0,0),
               X=X,
               ncov=ncol(X)
             )
           )

#mod <- stan_model('lcaPoisCov.stan')

if(!small) save(sdat,file='lcaDatMod.RData')
if(small) save(sdat,file='smallSdat.RData')
