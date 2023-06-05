data {
  int<lower=1> N;
  int<lower=1> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int<lower=1> p;
  matrix[N,p] X;
  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] E;           // exposure 
  real<lower=0> scaling_factor; // scales the variance of the spatial effects
}

transformed data {
  vector[N] log_E = log(E);
}

parameters {
  real beta0;            // intercept
  vector[p] beta;
  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  vector[N] theta;       // heterogeneous effects
  vector[N] phi;         // spatial effects
}
transformed parameters {
  vector[N] convolved_re;
  // variance of each component should be approximately equal to 1
  convolved_re =  sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;
}
model {
  y ~ poisson_log(log_E + beta0 + X*beta + convolved_re*sigma); 

  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  sum(phi) ~ normal(0.0, 0.001*N);

  beta0 ~ normal(0.0, 1.0);
  for(j in 1:p){
    beta[j] ~ normal(0.0, 1.0);
  }
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  rho ~ beta(1.0, 1.0);
}
generated quantities {
  real logit_rho = log(rho / (1.0 - rho));
  vector[N] mu = log_E + beta0 + X * beta + convolved_re * sigma; // co-variates
  vector[N] fitted;
  
  for(i in 1:N){
    fitted[i] = poisson_rng( mu[i] );
  }
}
