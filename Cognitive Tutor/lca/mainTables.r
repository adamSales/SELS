
load('parSE.RData')


############# parameters for full model

SEs1 <- seList(invHes1)

names(ppp1$beta) <- colnames(sdat$X)

as.data.frame(cbind(Est=ppp1$beta,SE=SEs1$beta,pval=2*pnorm(-abs(ppp1$beta/SEs1$beta)))%>%round(3))%>%
    rownames_to_column('Covariate')%>%
    flextable()%>%
    save_as_docx(path='../../../Manuscript/lcaTabs/LCAbetas.docx')

map_dfr(
    pars,
    classParFun,
    ppp=ppp1,SEs=SEs1,invHes=invHes1,digits=3
)%>%
    bind_rows(tibble(Parameter="Probability",
                     `Class 1`=sprintf("%.3f", round(mean(ppp1$nu),3)),
                     `Class 2`=sprintf("%.3f", round(1-mean(ppp1$nu),3))))%>%
    flextable()%>%
    save_as_docx(path='../../../Manuscript/lcaTabs/LCApars.docx')



