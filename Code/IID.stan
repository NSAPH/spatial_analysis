data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[p] beta;
  vector[N] s;         // latent effects
  real<lower=0> sigma_s;        // marginal standard deviation of spatial effects
}

transformed parameters {
  vector[N] b; // latent effect
  b =  sigma_s*s;
}

model {
  y ~ poisson_log(log_E + beta0 + X*beta + b); 
  s ~ normal(0.0, 1.0);
  
  for(j in 1:p){
    beta[j] ~ normal(0.0, 1.0);
  }
  
  beta0 ~ normal(0.0, 10.0);
  sigma_s ~ normal(0.0,1.0);
}

generated quantities {
  vector[N] mu=exp(log_E + beta0 + X*beta + b);
  vector[N] lik;
  vector[N] log_lik;
  
  for(i in 1:N){
    lik[i] = exp(poisson_lpmf(y[i] | mu[i] ));
    log_lik[i] = poisson_lpmf(y[i] | mu[i] );
  }
}
