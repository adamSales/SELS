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

  matrix[nstud,ncov] X;

  vector[3] zeros;
}

parameters {

  real meanTime[2];
  real<lower=0> sigTime[2];
  real effHint[2];
  real effErr[2];

  vector[3] probEff[nprob]; // hint, err, time

  corr_matrix[3] OmegaProb;
  vector<lower = 0>[3] sigProb;
  
  real alpha;
  
  real studEff[nstud];
  real<lower=0> sigStud;
  
  vector[ncov] beta;
}
transformed parameters {
  real<lower=0,upper=1> nu[nstud,2];
  vector[nstud] x_beta=X*beta;

  cov_matrix[3] SigmaProb=quad_form_diag(OmegaProb, sigProb);
 
  for(i in 1:nstud){
   nu[i,1]=inv_logit(alpha+x_beta[i]+studEff[i]);
   nu[i,2]=1-nu[i,1];
  }
}

model{
 real ps[2];

 studEff~normal(0,sigStud);

 probEff~multi_normal(zeros,SigmaProb);

 for(w in 1:nworked){
  for(c in 1:2){
   ps[c]=
    log(nu[stud[w],c])+
    poisson_log_lpmf(hint[w]|probEff[prob[w]][1]+effHint[c])+
    poisson_log_lpmf(err[w]|probEff[prob[w]][2]+effErr[c])+
    normal_lpdf(ltime[w]| probEff[prob[w]][3]+meanTime[c],sigTime[c]);
  }
  target += log_sum_exp(ps);
 } 
}  