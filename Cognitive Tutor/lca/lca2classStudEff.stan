// adapted from https://discourse.mc-stan.org/t/latent-class-model-estimation-with-long-format-data/1799
// "Stan code for wide format"

data {
  int<lower=1> nworked;	//number of worked items--rows of data
  int<lower=1> nprob;                                     // number of items
  int<lower=1> nstud;                                     // number of respondents
  int<lower=1> ncov;					  // number of person-level covariates
  int<lower=0> hint[nworked];                                     // score for obs n
  int<lower=0> err[nworked];
  real ltime[nworked];
  int<lower=1,upper=nprob> prob[nworked];
  int<lower=1,upper=nstud> stud[nworked];

  //matrix[nstud,ncov] X;

  vector[3] zeros;
}

parameters {

  real meanTime[2];
  real<lower=0> sigTime[2];
  real effHint[2];
  real effErr[2];

  //vector[3] probEff[nprob]; // hint, err, time

  //corr_matrix[3] OmegaProb;
  //vector[3] sigProb;
  
  real alpha;
  
  real<lower=0> sigStud;

  //vector[ncov] beta;

  ordered[2] cHint;
  ordered[2] cErr;

vector[nstud] studEff;

}
//transformed parameters {

  
  //cov_matrix[3] SigmaProb=quad_form_diag(OmegaProb, sigProb);
//}

model{
  
  //vector[nstud] nu=inv_logit(studEff);

// priors
 meanTime~normal(0,5);
 sigTime~normal(0,5);
 effHint~normal(0,5);
 effErr~normal(0,5);
 //sigProb~normal(0,1);
 sigStud~normal(2,1);
 //to_vector(beta)~normal(0,1);
 alpha~normal(0,5);

 studEff~normal(alpha,sigStud);

 //probEff~multi_normal(zeros,SigmaProb);

 for(w in 1:nworked){
   real nu=inv_logit(studEff[stud[w]]);
   target += log_sum_exp(
    log(nu)+
    ordered_logistic_lpmf(hint[w]|effHint[1],cHint)+
    ordered_logistic_lpmf(err[w]|effErr[1],cErr)+
    normal_lpdf(ltime[w]| meanTime[1],sigTime[1]),
    log(1-nu)+
    ordered_logistic_lpmf(hint[w]|effHint[2],cHint)+
    ordered_logistic_lpmf(err[w]|effErr[2],cErr)+
    normal_lpdf(ltime[w]| meanTime[2],sigTime[2])
   );
  }
}  