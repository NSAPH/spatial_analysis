data {
  int<lower=1> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
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
  real<lower=0> tau_theta;   // precision of heterogeneous effects
  real<lower=0> tau_phi;     // precision of spatial effects
  vector[N] theta;       // heterogeneous effects
  vector[N] phi;         // spatial effects
}
transformed parameters {
  real<lower=0> sigma_theta = inv(sqrt(tau_theta));  // convert precision to sigma
  real<lower=0> sigma_phi = inv(sqrt(tau_phi));      // convert precision to sigma
}
model {
  
  y ~ poisson_log(log_E + beta0 + X*beta + sigma_phi*phi+ theta * sigma_theta); 

  // NOTE:  no prior on phi_raw, it is used to construct phi
  // the following computes the prior on phi on the unit scale with sd = 1
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)  
  
  beta0 ~ normal(0, 5);
  for(j in 1:p){
    beta[j] ~ normal(0.0, 1.0);
  }
  theta ~ normal(0, 1);
  tau_theta ~ gamma(3.2761, 1.81);  // Carlin WinBUGS priors
  tau_phi ~ gamma(1, 1);            // Carlin WinBUGS priors 
}
generated quantities {
  vector[N] mu=exp(log_E + beta0 + X*beta + sigma_phi*phi)+sigma_theta*theta;
  vector[N] fitted;
  vector[N] log_lik;

  for(i in 1:N){
    fitted[i] = poisson_rng( mu[i] );
    log_lik[i] = poisson_lpmf(y[i] | mu[i] );
  }
}
