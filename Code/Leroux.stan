data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=1> N_edges;
  int<lower=1, upper=N> node1[N_edges]; 
  int<lower=1, upper=N> node2[N_edges];  
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[p] beta;
  real<lower=0> sigma;     // conditional std dev (=1/sqrt(tau)=sqrt(delta))
  real<lower=0, upper=1> lambda; // mixing parameter
  vector[N] s;         // spatial effects
}

transformed parameters {
  vector[N] b; 
  b = sigma*s;
}

model {
  y ~ poisson_log(log_E + beta0 + X*beta + b);

  for(j in 1:p){
    beta[j] ~ normal(0.0, 10.0);
  }

  beta0 ~ normal(0.0, 10.0);
  sigma ~ normal(0.0,1.0);
  lambda ~ uniform(0.0,1.0);
  target += -0.5 *( (1-lambda)*dot_self(s) + lambda*dot_self(s[node1] - s[node2]) );
}

generated quantities {
  vector[N] mu = exp(log_E + beta0 + X*beta + b);
  vector[N] log_lik ;
  
  for(i in 1:N){
    log_lik[i]= poisson_lpmf(y[i] | mu[i]);
  }
}