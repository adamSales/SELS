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
  int<lower=1> nclass;
}

parameters {

  real meanTime[nclass];
  real sigTime[nclass];
  real effHint[nclass];
  real effErr[nclass];

  vector[3] probEff[nprob]; // hint, err, time

  corr_matrix[3] OmegaProb;
  vector[3] sigProb;
  
  vector[nclass-1] alpha;
  
  vector[nclass-1] studEff[nstud];
  corr_matrix[nclass-1] OmegaStud;
  vector[nclass-1] sigStud;

  matrix[ncov,nclass-1] beta;

}
transformed parameters {
  vector[nclass] nuRaw[nstud];
  simplex[nclass] nu[nstud];
  matrix[nstud,nclass-1] x_beta=X*beta;
  
  cov_matrix[3] SigmaProb=quad_form_diag(OmegaProb, sigProb);
  cov_matrix[nclass-1] SigmaStud=quad_form_diag(OmegaStud,sigStud);
 
  for(i in 1:nstud){
   nuRaw[i]=append_row(0,alpha+x_beta[i]'+studEff[i]);
   nu[i]=softmax(nuRaw[i]);
  }
}

model{

  real ps[nclass];
  
// priors
 meanTime~normal(0,5);
 sigTime~normal(0,5);
 effHint~normal(0,5);
 effErr~normal(0,5);
 sigProb~normal(0,1);
 sigStud~normal(0,1);
 to_vector(beta)~normal(0,1);

 studEff~multi_normal(rep_vector(0,nclass-1),SigmaStud);

 probEff~multi_normal(zeros,SigmaProb);

 for(w in 1:nworked){
  for(c in 1:nclass){
   ps[c]=
    log(nu[stud[w]][c])+
    poisson_log_lpmf(hint[w]|probEff[prob[w]][1]+effHint[c])+
    poisson_log_lpmf(err[w]|probEff[prob[w]][2]+effErr[c])+
    normal_lpdf(ltime[w]| probEff[prob[w]][3]+meanTime[c],sigTime[c]);
  }
  target += log_sum_exp(ps);
 } 
}  