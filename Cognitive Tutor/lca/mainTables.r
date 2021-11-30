
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



#### posterior probabilities for example students
## "10390705014x" 57 observations:
stud643 <- 1-postStud( "852",sdat,ppp1)
## 30010304012 79 obs
stud2117 <- 1-postStud("1207",sdat,ppp1)
## 50010102012 50 obs
stud1081 <- 1-postStud("1771",sdat,ppp1)

save(stud643,stud2117,stud1081,file='exampleStudents.RData')
