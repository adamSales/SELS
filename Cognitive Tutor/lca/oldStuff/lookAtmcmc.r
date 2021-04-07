pars <- extract(fit3,permute=FALSE,inc_warmup=TRUE)
lp <- extract(fit3,par='lp__',permuted=FALSE,inc_warmup=TRUE)
lp[,,1]%>%as.data.frame()%>%mutate(iter=1:n())%>%pivot_longer(-iter,names_to='chain',values_to='lp')%>%ggplot(aes(iter,lp,color=chain))+geom_line()
lp[100:2000,,1]%>%as.data.frame()%>%mutate(iter=1:n())%>%pivot_longer(-iter,names_to='chain',values_to='lp')%>%ggplot(aes(iter,lp,color=chain))+geom_line()
lp[1500:2000,,1]%>%as.data.frame()%>%mutate(iter=1:n())%>%pivot_longer(-iter,names_to='chain',values_to='lp')%>%ggplot(aes(iter,lp,color=chain))+geom_line()+geom_smooth(se=FALSE)

inits <- fit3$inits

head(which(startsWith(colnames(pars),'nu')))
unique(substr(colnames(pars)[6263:29100],1,2))
colnames(pars)[6262]

pars <- pars[,1:6262]
parNames <- map_chr(strsplit(colnames(pars),'[',fixed=TRUE),~.[1])
table(parNames)

randPars <- pars[,endsWith(parNames,'Eff')]
fixedPars <- pars[,!endsWith(parNames,'Eff')]

#par(ask=TRUE)
for(pn in unique(parNames[!endsWith(parNames,'Eff')])){
    fp <- fixedPars[,startsWith(colnames(fixedPars),pn)]
    print(
        round(
            cbind(
                tau=apply(fp,2,function(x) cor(x,1:501,method='kendall')),
                p=apply(fp,2,function(x) cor.test(x,1:501,method='kendall')$p.value)
            ),3
        )
    )

    p <- fp%>%
        as.data.frame()%>%
        mutate(iter=1:n())%>%
        pivot_longer(-iter,names_to='par')%>%
        ggplot(aes(iter,value,color=par))+geom_line()+geom_smooth(se=FALSE)+geom_smooth(method='lm')
    print(p)



}


mcmcPar <- list()
for(pn in unique(parNames[!endsWith(parNames,'Eff')]))
    mcmcPar[[pn]]  <- colMeans(pars[,parNames==pn])

ord <- order(mcmcPar[[1]])
for(i in 1:4) mcmcPar[[i]] <- mcmcPar[[i]][ord]

### compare to mle
sdat0 <- sdat
print(load('lca/fit3class.RData'))

parNames2 <- map_chr(strsplit(names(fit$par),'[',fixed=TRUE),~.[1])

mlePar <- list()
for(i in 1:12)
    mlePar[[unique(parNames2)[i]]] <- fit$par[parNames2==unique(parNames2)[i]]

ord <- order(mlePar[[1]])
for(i in 1:4) mlePar[[i]] <- mlePar[[i]][ord]


head(which(startsWith(names(fit$par),'nu')))
unique(substr(names(fit$par)[6530:length],1,2))
names(fit$par)[6262]

pars2 <- fit$par[,1:6262]

table(parNames)

randPars2 <- pars2[,endsWith(parNames2,'Eff')]
fixedPars2 <- pars2[,!endsWith(parNames2,'Eff')]
