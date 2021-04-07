
load('parSE.RData')

################### parameters for other sig_stud
SEs1 <- seList(invHes1)

altTab1 <- altClassTab(parList=pars,
                       ppp=ppp1,SEs=SEs1,invHes=invHes1,sig='1')



for(sig in c('0','.5','2','5')){
    ppp <- get(paste0('ppp',sig))
    invHes <- get(paste0('invHes',sig))
    SEs <- seList(invHes)
    names(ppp$beta) <- colnames(sdat$X)

    altTab1 <- bind_rows(
        altTab1,
        altClassTab(parList=pars,ppp=ppp,SEs=SEs,invHes=invHes,sig=sig)
    )

}

estSE <- function(est,se,digits)
    paste0(sprintf(paste0("%.",digits,"f"), round(est,digits)),' (',
                         sprintf(paste0("%.",digits,"f"), round(se,digits)),')')

ind <- rep(5,length(pars)+1)
names(ind) <- c(c(
            meanTime='Time (mean)',
            sigTime='Time (SD)',
            effErr='Error',
            effHint='Hint')[pars],'Probability')

sink('appendix2table1.tex')
altTab1%>%
    mutate(        sigStud=ifelse(sigStud=='.5','0.5',sigStud))%>%
    arrange(factor(Parameter,levels=c(pars,'Probability')),sigStud)%>%
    transmute(
        "$SD(\\eta)$"=sigStud,
        `Class 1`=ifelse(Parameter=='Probability',sprintf("%.3f",round(`Class 1`,3)),estSE(`Class 1`,SE1,3)),
        `Class 2`=ifelse(Parameter=='Probability',sprintf("%.3f",round(`Class 2`,3)),estSE(`Class 2`,SE2,3)),
        Difference=ifelse(Parameter=='Probability','',estSE(Difference,SEdiff,3))
    )%>%
    kbl(format='html',booktabs=TRUE,escape=FALSE)%>%
    pack_rows(index=ind)

sink()

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
    filter(which!='Diff')%>%
    mutate(
        which=paste('Class',which),
        Parameter=c(
            meanTime='Time (mean)',
            sigTime='Time (SD)',
            effErr='Error',
            effHint='Hint',
            Probability='Probability')[Parameter],
        Parameter=factor(Parameter,levels=unique(Parameter))
    )%>%
    ggplot(aes(
        x=which,
        y=Est,ymin=Est-2*SE,ymax=Est+2*SE,color=sigStud,group=sigStud))+
    geom_point(position=position_dodge(width=0.2))+
    geom_errorbar(position=position_dodge(width=0.2),width=0)+#"dodge")+
    facet_wrap(~Parameter,ncol=1,scales="free_y")+
    labs(x=NULL,y="Estimate",color=expression(paste("SD(",eta,")")))
ggsave('measurementParCompare.jpg')






altTab2 <- altCoefTab1(ppp=ppp1,SEs=SEs1,sig='1')

for(sig in c('0','.5','2','5')){
    ppp <- get(paste0('ppp',sig))
    invHes <- get(paste0('invHes',sig))
    SEs <- seList(invHes)
    names(ppp$beta) <- colnames(sdat$X)

    altTab2 <- bind_rows(
        altTab2,
        altCoefTab1(ppp=ppp,SEs=SEs,sig=sig)
    )
}

tr <- altTab2%>%
    mutate(        sigStud=ifelse(sigStud=='.5','0.5',sigStud))%>%
    group_by(sigStud)%>%
    group_map(~createTexreg(coef.names=.$covariate,coef=.$Est,se=.$SE,pvalues=2*pnorm(-abs(.$Est/.$SE))))

names(tr) <- paste0("$SD(\\eta)=$",sort(as.numeric(unique(altTab2$sigStud))))

texreg(tr,file='Appendix2coefTable.tex',digits=3)
htmlreg(tr,file='Appendix2coefTable.html',digits=3)

ind <- rep(5,length(ppp1$beta))
names(ind) <- names(ppp1$beta)


altTab2%>%
    mutate(sigStud=ifelse(sigStud=='.5','0.5',sigStud))%>%
ggplot(aes(factor(covariate,levels=colnames(sdat$X)),Est,color=sigStud,group=sigStud,ymin=Est-2*SE,ymax=Est+2*SE))+
    geom_point(position=position_dodge(width=0.2))+
    geom_errorbar(position=position_dodge(width=0.2),width=0)+
    geom_hline(yintercept=0)+
    labs(x=NULL,y='Coefficient Estimate',color=expression(paste("SD(",eta,")")))
ggsave('coefCompare.jpg')

