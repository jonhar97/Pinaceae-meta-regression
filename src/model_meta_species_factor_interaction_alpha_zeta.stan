data {
  int<lower=1> N; // number of ID estimates 
  //int<lower=1> K; // number of Inbreeding levels
  int<lower=1> J; // number of life stages
  int<lower=1> M; // number of population fragmantation groups
  int<lower=1> JM; // number of interaction life stage x population levels
  int<lower=1> Nsp; // number of species
  int<lower=1> L; // number of studies 
  real F[N]; // inbreeding coeficient levels
  //int<lower=1,upper=K> incidence_F[N]; // incidence vector relating ID est (y) to inbreeding levels (alpha) 
  int<lower=1,upper=J> incidence_L[N]; // incidence vector relating ID est (y) to life stages (beta)
  int<lower=1,upper=M> incidence_G[N]; // incidence vector relating ID est (y) to degree of population fragmentation
  int<lower=1,upper=JM> incidence_LG[N]; // incidence vector for interaction group effect 
  int<lower=1,upper=L> group[N]; // incidence vector relating ID est to study (group effect) 
  int<lower=1,upper=Nsp> group_sp[N]; // incidence vector for species 
  real<lower=-1,upper=2> y[N]; // estimated ID 
  real<lower=0> SE[N]; // s.e. of effect estimates 
}
parameters {
  real<lower=2> nu;
  real mu; 
  real alpha;
  vector[J-1] beta0;
  vector[M-1] gamma0;
  vector[JM-1] delta0;
  vector[J-1] epsilon0;
  vector[M-1] zeta0;
  vector[L] a;
  vector[Nsp] b;
  real<lower=0.0001,upper=100> sigma_y;	
  real<lower=0.0001,upper=100> sigma_g;
  real<lower=0.0001,upper=100> sigma_sp;

}
transformed parameters {
  //real Fst_std[Nsp];
  real theta[N];
  vector[J] beta;
  vector[M] gamma;
  vector[JM] delta;
  vector[J] epsilon;
  vector[M] zeta;

  for (j in 1:(J-1))
    beta[j] = beta0[j];
  beta[4] = -sum(beta0);
   // population fragmentation level

  gamma[1] = gamma0[1];
  gamma[2] = -gamma0[1];
// interaction population fragmentation level and life-history 
  for (j in 1:(JM-1))
    delta[j] = delta0[j];
  delta[JM] = -sum(delta0);

  for (j in 1:(J-1))
    epsilon[j] = epsilon0[j];
  epsilon[4] = -sum(epsilon0);
  
  // interaction population fragmentation and inbreeding level
  zeta[1] = zeta0[1];
  zeta[2] = -zeta0[1];
 //for (j in 1:Nsp)
  //Fst_std[j] = (Fst[j]-mean(Fst))/sd(Fst);

 for (j in 1:N)
    theta[j] = mu + a[group[j]] + b[group_sp[j]] + alpha*F[j] + beta[incidence_L[j]] + gamma[incidence_G[j]] + delta[incidence_LG[j]] + epsilon[incidence_L[j]]*F[j] + zeta[incidence_G[j]]*F[j]; // //*delta0[incidence_L[j]]; //  b[group_sp[j]] + c[group_lg[j]] + 

}
model {
  // priors
 // for (k in 1 : (K-1)) 
   // alpha0[k] ~ normal(0, 10); // vague priors
  alpha ~ normal(0, 10);	
  for (k in 1 : (J-1))
    beta0[k] ~ normal(0, 1); // vague priors 
  for (k in 1 : (M-1))
    gamma0[k] ~ normal(0, 1); // vague priors 
  for (k in 1 : (JM-1))
    delta0[k] ~ normal(0, 1); // vague priors 
  for (k in 1 : (J-1))
    epsilon0[k] ~ normal(0, 1); // vague priors 
  for (k in 1 : (M-1))
    zeta0[k] ~ normal(0, 1); // vague priors 
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

  vector[N] log_lik;
  real Ftmp[4] = {0.125, 0.25, 0.5, 0.75};
  vector[4] ypred;
  real mean_beta;
  real mean_SE;
  real mean_a;
  real mean_b;
  int j;
  real pred;
  //real pred1;
 // real pred2;


  mean_beta = mean(beta);
  mean_SE = mean(SE);
  mean_a = mean(a);
  mean_b = mean(b);

  for (n in 1:N)
    log_lik[n] = student_t_lpdf(y[n] | nu,theta[n], sigma_y*SE[n]);
  for (i in 1:4)  {
    //j = i + 4;
    pred = mean_a + mean_b + alpha*Ftmp[i] + beta[1] + mu + gamma[1] + delta[1] + epsilon[1]*Ftmp[i];
   // pred2 = pred + gamma[2] + delta[5];
    ypred[i] = student_t_rng(nu,pred, sigma_y*mean_SE);
   // ypred[j] = student_t_rng(nu,pred2, sigma_y*mean_SE);
  }

}
