library(nimble)
#######################
#    NIMBLE: Simple   #
#######################
Simple_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    log(proba[i]) <- log(E[i]) + beta0 + inprod(X[i,1:p],beta[1:p])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
#######################
#    NIMBLE: IID   #
#######################
IID_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    log(proba[i]) <- log(E[i]) + beta0 + b[i] + inprod(X[i,1:p],beta[1:p])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  for(i in 1:N){
    b[i] ~ dnorm(beta0, sd=sigma_b)
  }
  sigma_b ~ T(dnorm(0,sd=1),0,)
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
#######################
#     NIMBLE: CAR     #
#######################
CAR_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    log(proba[i]) <- log(E[i]) + beta0 + b[i] + inprod(X[i,1:p],beta[1:p])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  
  b[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean=1)
  tau <- 1/(sigma_b^2)
  sigma_b ~ T(dnorm(0,sd=1),0,)
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
#######################
#     NIMBLE: BYM     #
#######################
BYM_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    log(proba[i]) <- log(E[i]) + beta0 + b[i] + inprod(X[i,1:p],beta[1:p])
    b[i] <- theta[i] + s[i]
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  
  s[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean=1)
  tau <- 1/(sigma_s^2)
  sigma_s ~ T(dnorm(0,sd=1),0,)
  for(i in 1:N){
    theta[i] ~ dnorm(beta0, sd=sigma_theta)
  }
  sigma_theta ~ T(dnorm(0,sd=1),0,)
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
#######################
#     NIMBLE: BYM2    #
#######################
BYM2_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    proba[i] <- exp(log(E[i])+ beta0 + b[i] + inprod(X[i,1:p],beta[1:p]) )
    b[i] <- sigma*(sqrt(1-lambda)*theta[i] + sqrt(lambda/scaling_factor)*s[i])
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  
  s[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean=1)
  sigma ~ T(dnorm(0,sd=1),0,)
  for(i in 1:N){
    theta[i] ~ dnorm(0, sd=1)
  }
  lambda ~ dunif(0,1)
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
#######################
#   NIMBLE: Leroux    #
#######################
Leroux_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    y[i] ~ dpois(proba[i])
    proba[i] <- exp(log(E[i]) + beta0 + b.centred[i] + inprod(X[i,1:p],beta[1:p]) )
    # proba[i] <- exp(log(E[i])+ beta0 + b[i] + inprod(X[i,1:p],beta[1:p]) )
    b.centred[i] <- b[i] - b.mean
  }
  
  # Priors
  beta0 ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta[k] ~ dnorm(0, sd=1)
  }
  
  PrecMat[1:N, 1:N] <- (lambda*Q[1:N, 1:N] + (1 - lambda)*diag(N))/(pow(sigma,2))
  mean_centred[1:N] <- beta0*mu_b[1:N]
  b[1:N] ~ dmnorm(mean = mean_centred[1:N], prec = PrecMat[1:N,1:N])
  b.mean <- sum(b[1:N])/N
  
  sigma ~ T(dnorm(0,sd=1),0,)
  
  lambda ~ dunif(0,1)
  
  # Fitted values and likelihood for WAIC
  for(i in 1:N){
    fitted[i] ~ dpois(proba[i])
  }
})
