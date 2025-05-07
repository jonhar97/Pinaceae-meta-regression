data {
  int<lower=1> N; // number of ID estimates 
  //int<lower=1> K; // number of Inbreeding levels
  int<lower=1> J; // number of life stages
  //int<lower=1> M; // number of population fragmantation groups
  //int<lower=1> JM; // number of interaction life stage x population levels
  int<lower=1> Nsp; // number of species
  int<lower=1> L; // number of studies 
  real F[N]; // inbreeding coeficient levels
  //int<lower=1,upper=K> incidence_F[N]; // incidence vector relating ID est (y) to inbreeding levels (alpha) 
  int<lower=1,upper=J> incidence_L[N]; // incidence vector relating ID est (y) to life stages (beta)
  //int<lower=1,upper=M> incidence_G[N]; // incidence vector relating ID est (y) to degree of population fragmentation
  //int<lower=1,upper=JM> incidence_LG[N]; // incidence vector for interaction group effect 
  int<lower=1,upper=L> group[N]; // incidence vector relating ID est to study (group effect) 
  int<lower=1,upper=Nsp> group_sp[N]; // incidence vector for species 
  real Fst[Nsp]; //Fst covariate
  real<lower=-1,upper=2> y[N]; // estimated ID 
  real<lower=0> SE[N]; // s.e. of effect estimates 
}
parameters {
  real<lower=2> nu;
  real mu; 
  //vector[K-1] alpha0;
  real alpha;
  vector[J-1] beta0;
  //vector[M-1] gamma0;
  vector[J-1] delta0;
  vector[L] a;
  vector[Nsp] b;
//  vector[JM] c;
  real<lower=0.0001,upper=100> sigma_y;	
  real<lower=0.0001,upper=100> sigma_g;
  real<lower=0.0001,upper=100> sigma_sp;
 // real<lower=0.0001,upper=100> sigma_lg;
  real gamma;

}
transformed parameters {
  real Fst_std[Nsp];
  real theta[N];
  //vector[K] alpha;
  vector[J] beta;
  //vector[M] gamma;
  vector[J] delta;

  //alpha[1] = 0; // zero contrast for baseline Inbreeding level
  //for (j in 1:(K-1))
    //alpha[j+1] = alpha0[j];
  beta[1] = beta0[1];
  beta[2] = 0; // zero contrast for baseline life-stage level
  for (j in 2:(J-1))
    beta[j+1] = beta0[j];
  //gamma[2] = 0; // zero contrast for baseline population fragmentation level
  //for (j in 1:(M-1))
  //gamma[1] = gamma0[1];
 delta[1] = delta0[1];
 delta[2] = 0;// zero contrast for baseline interaction between population fragmentation level and life-history stage
  for (j in 2:(J-1))
    delta[j+1] = delta0[j];
 for (j in 1:Nsp)
  Fst_std[j] = (Fst[j]-mean(Fst))/sd(Fst);

 for (j in 1:N)
    theta[j] = mu + a[group[j]] + b[group_sp[j]] + alpha*F[j] + beta[incidence_L[j]] + gamma*Fst_std[group_sp[j]] + delta[incidence_L[j]]*Fst_std[group_sp[j]]; //*delta0[incidence_L[j]]; //  b[group_sp[j]] + c[group_lg[j]] + 

}
model {
  // priors
 // for (k in 1 : (K-1)) 
   // alpha0[k] ~ normal(0, 10); // vague priors
  //alpha ~ normal();	
  for (k in 1 : (J-1))
    beta0[k] ~ normal(0, 10); // vague priors 
 // for (k in 1 : (M-1))
 //   gamma0[k] ~ normal(0, 10); // vague priors 
  for (k in 1 : (J-1))
    delta0[k] ~ normal(0, 10); // vague priors 
  for (j in 1:L)
    a[j] ~ normal(0,sigma_g); // group coefficient
  for (j in 1:Nsp)
    b[j] ~ normal(0,sigma_sp); // group coefficient for species 
  //for (j in 1:JM)
    //c[j] ~ normal(0,sigma_lg); // group coefficient for interaction   
  mu ~ normal(0,10); // overall intercept
  nu ~ gamma(2,0.1); // degrees of freedom parameter
  // likelihood	
  for (j in 1:N)
    y[j] ~ student_t(nu,theta[j], sigma_y*SE[j]); // weighted regression
}
generated quantities {
  //real alpha_n[K];
  real beta_n[J];
  real delta_n[J];
  //real mean_alpha;
  real mean_beta;
  real mean_delta;
  vector[N] log_lik;

 // mean_alpha = mean(alpha);
 // for (i in 1:K) 
   // alpha_n[i] = alpha[i] - mean_alpha;
  mean_beta = mean(beta);
  for (i in 1:J) 
    beta_n[i] = beta[i] - mean_beta;  
 mean_delta = mean(delta0);
 for (i in 1:J) 
   delta_n[i] = delta[i] - mean_delta; 
  for (n in 1:N)
    log_lik[n] = student_t_lpdf(y[n] | nu,theta[n], sigma_y*SE[n]);

}

