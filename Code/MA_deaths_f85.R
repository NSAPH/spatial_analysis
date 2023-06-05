library(nimble)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###########################
#        SIMPLE project       #
###########################
fitSimple_deaths_f85 = stan("Simple.stan",
                            data=list(N=N,y=y_deaths_f85,E=E_f85, X=X_agg, p=p),
                            warmup=1e3, iter=1e4, chains=4, thin=1)
traceplot(fitSimple_deaths_f85, pars = c("beta"))
summary(fitSimple_deaths_f85, pars = c("beta"))
plot(fitSimple_deaths_f85, pars=c("beta"))
samples_Simple=extract(fitSimple_deaths_f85)

log_lik_simple_deaths_f85 <- extract_log_lik(fitSimple_deaths_f85)
loo_simple_deaths_f85 <- loo(log_lik_simple_deaths_f85)$estimate
WAIC_simple_deaths_f85 <- waic(log_lik_simple_deaths_f85)$estimate

saveRDS(fitSimple_deaths_f85, "fitSimple_f85.rds")


fitSimple_deaths_f85_thin = stan("Simple.stan",
                                 data=list(N=N,y=y_deaths_f85,E=E_f85, X=X_agg, p=p),
                                 warmup=1e3, iter=2e3, chains=4, thin=1)
traceplot(fitSimple_deaths_f85_thin, pars = c("beta"))
summary(fitSimple_deaths_f85_thin)$summary[,"Rhat"]

###########################
#        IID project      #
###########################
fitIID_deaths_f85 = stan("IID.stan",
                         data=list(N=N,y=y_deaths_f85,E=E_f85, X=X_agg, p=p),
                         warmup=1e3, iter=1e4, chains=4, thin=1)

# traceplot(fitIID_deaths_f85, pars = c("beta"))
# plot(fitIID_deaths_f85, pars = c("beta"))

saveRDS(fitIID_deaths_f85, "fitIID_deaths_f85.rds")

###########################
#        CAR stan         #
###########################
fitCAR_deaths_f85 <- stan('CAR.stan', data = list(n = N,         
                                                  p = p,         
                                                  X=X_agg,               
                                                  y = y_deaths_f85,             
                                                  E=E_f85, 
                                                  W_n = sum(W) / 2,   
                                                  W = W),
                          warmup=1e3, iter=1e4, chains = 4)
# traceplot(fitCAR_deaths_f85, par = c("beta", "tau"))
# plot(fitCAR_deaths_f85, par = c('beta'))

saveRDS(fitCAR_deaths_f85, "fitCAR_deaths_f85.rds")

###########################
#        BYM stan        #
###########################
fitBYM_deaths_f85 = stan("BYM.stan",
                         data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y_deaths_f85,p=p,X=X_agg,E=E_f85),
                         control=list(adapt_delta = 0.97, stepsize = 0.1),
                         chains=4, warmup=1e3, iter=1e4);

# traceplot(fitBYM_deaths_f85, par = c("beta"))
# summary(fitBYM_deaths_f85, par = c("beta"))
# plot(fitBYM_deaths_f85, par = c("beta"))

saveRDS(fitBYM_deaths_f85, "fitBYM_deaths_f85.rds")

###########################
#        BYM2 stan        #
###########################
fitBYM2_deaths_f85 = stan("BYM2.stan",
                          data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y_deaths_f85,X=X_agg,E=E_f85,scaling_factor=scaling_factor,K=p),
                          warmup=1e3, iter=1e4, chains=4);
# traceplot(fitBYM2_deaths_f85, pars=c("beta", "rho"))
# plot(fitBYM2_deaths_f85, par=c("beta"))
# summary(fitBYM2_deaths_f85, par=c("beta", "rho"))

saveRDS(fitBYM2_deaths_f85, "fitBYM2_deaths_f85.rds")


###########################
#        Leroux stan      #
###########################
fitLeroux_deaths_f85 = stan("Leroux.stan",
                            data=list(N=N,y=y_deaths_f85,E=E_f85, X=X_agg, p=p, N_edges=N_edges, node1=node1, node2=node2),
                            warmup=1e3, iter=1e4, chains=4, thin=1)
# traceplot(fitLeroux_deaths_f85, pars = c("beta", "lambda"))
# summary(fitLeroux_deaths_f85, pars = c("beta", "lambda"))
# plot(fitLeroux_deaths_f85, pars = c("beta", "lambda"))

saveRDS(fitLeroux_deaths_f85, "fitLeroux_deaths_f85.rds")


#######################
#    NIMBLE: Simple   #
#######################
Consts_Simple_deaths_f85 = list(N=N,
                                E=E_f85,
                                X=X_agg,
                                p=p)
Data_Simple_deaths_f85 = list(y=y_deaths_f85)
Inits_Simple_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), y=rep(1,N))
Simple_deaths_f85_Model = nimbleModel(code=Simple_Code, 
                                      constants=Consts_Simple_deaths_f85, 
                                      data=Data_Simple_deaths_f85, 
                                      inits=Inits_Simple_deaths_f85)
Simple_deaths_f85_conf <- configureMCMC(Simple_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y"), enableWAIC=T)
Simple_deaths_f85_MCMC=buildMCMC(Simple_deaths_f85_conf)
Simple_deaths_f85_Cmodel=compileNimble(Simple_deaths_f85_Model)
Simple_deaths_f85_Cmcmc=compileNimble(Simple_deaths_f85_MCMC, project = Simple_deaths_f85_Model)
Simple_deaths_f85_runMCMC=runMCMC(Simple_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# plot(Simple_deaths_f85_runMCMC$samples[,c("beta[1]")])
# beta_Simple_deaths_f85 <- Simple_deaths_f85_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                  "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                  "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(Simple_deaths_f85_runMCMC)) {
#   colnames(beta_Simple_deaths_f85[[i]]) <- names
# }
# summary(beta_Simple_deaths_f85)

saveRDS(Simple_deaths_f85_runMCMC, "Simple_deaths_f85_runMCMC.rds")

#######################
#    NIMBLE: IID   #
#######################
Consts_IID_deaths_f85 = list(N=N,
                             E=E_f85,
                             X=X_agg,
                             p=p)
Data_IID_deaths_f85 = list(y=y_deaths_f85)
Inits_IID_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), sigma_b=1, b=rep(0,N), y=rep(1,N))

IID_deaths_f85_Model = nimbleModel(code=IID_Code, 
                                   constants=Consts_IID_deaths_f85, 
                                   data=Data_IID_deaths_f85, 
                                   inits=Inits_IID_deaths_f85)
IID_deaths_f85_conf <- configureMCMC(IID_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_b"), enableWAIC=T)
IID_deaths_f85_MCMC=buildMCMC(IID_deaths_f85_conf)
IID_deaths_f85_Cmodel=compileNimble(IID_deaths_f85_Model)
IID_deaths_f85_Cmcmc=compileNimble(IID_deaths_f85_MCMC, project = IID_deaths_f85_Model)
IID_deaths_f85_runMCMC=runMCMC(IID_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_IID_deaths_f85 <- IID_deaths_f85_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                            "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                            "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(IID_deaths_f85_runMCMC)) {
#   colnames(beta_IID_deaths_f85[[i]]) <- names
# }
# summary(beta_IID_deaths_f85)

saveRDS(IID_deaths_f85_runMCMC, "IID_deaths_f85_runMCMC.rds")


#######################
#     NIMBLE: CAR     #
#######################
Consts_CAR_deaths_f85 = list(N=N,
                             L=length(neigh),
                             adj=neigh,
                             weights=rep(1,length(neigh)),
                             num=numneigh,
                             E=E_f85,
                             X=X_agg,
                             p=p)
Data_CAR_deaths_f85 = list(y=y_deaths_f85)
Inits_CAR_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=y_deaths_f85, sigma_b=1, b=rep(0,N), y=y_deaths_f85)
CAR_deaths_f85_Model = nimbleModel(code=CAR_Code, 
                                   constants=Consts_CAR_deaths_f85, 
                                   data=Data_CAR_deaths_f85, 
                                   inits=Inits_CAR_deaths_f85)
CAR_deaths_f85_conf <- configureMCMC(CAR_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_b"), enableWAIC=T)
CAR_deaths_f85_MCMC=buildMCMC(CAR_deaths_f85_conf)
CAR_deaths_f85_Cmodel=compileNimble(CAR_deaths_f85_Model)
CAR_deaths_f85_Cmcmc=compileNimble(CAR_deaths_f85_MCMC, project = CAR_deaths_f85_Model)
CAR_deaths_f85_runMCMC=runMCMC(CAR_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_CAR_deaths <- CAR_deaths_f85_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                    "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                    "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(CAR_deaths_f85_runMCMC)) {
#   colnames(beta_CAR_deaths[[i]]) <- names
# }
# summary(beta_CAR_deaths)

saveRDS(CAR_deaths_f85_runMCMC, "CAR_deaths_f85_runMCMC.rds")


#######################
#     NIMBLE: BYM     #
#######################
Consts_BYM_deaths_f85 = list(N=N,
                             L=length(neigh),
                             adj=neigh,
                             weights=rep(1,length(neigh)),
                             num=numneigh,
                             E=E_f85,
                             X=X_agg,
                             p=p)
Data_BYM_deaths_f85 = list(y=y_deaths_f85)
Inits_BYM_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=y_deaths_f85, sigma_s=1, s=rep(0,N),sigma_theta=1, theta=rep(0,N),y=y_deaths_f85)
BYM_deaths_f85_Model = nimbleModel(code=BYM_Code, 
                                   constants=Consts_BYM_deaths_f85, 
                                   data=Data_BYM_deaths_f85, 
                                   inits=Inits_BYM_deaths_f85)
BYM_deaths_f85_conf <- configureMCMC(BYM_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_s", "sigma_theta", "s", "theta"), enableWAIC=T)
BYM_deaths_f85_MCMC=buildMCMC(BYM_deaths_f85_conf)
BYM_deaths_f85_Cmodel=compileNimble(BYM_deaths_f85_Model)
BYM_deaths_f85_Cmcmc=compileNimble(BYM_deaths_f85_MCMC, project = BYM_deaths_f85_Model)
BYM_deaths_f85_runMCMC=runMCMC(BYM_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_BYM_deaths_f85 <- BYM_deaths_f85_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                            "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                            "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(BYM_deaths_f85_runMCMC)) {
#   colnames(beta_BYM_deaths_f85[[i]]) <- names
# }
# summary(beta_BYM_deaths_f85)

saveRDS(BYM_deaths_f85_runMCMC, "BYM_deaths_f85_runMCMC.rds")

#######################
#     NIMBLE: BYM2    #
#######################
Consts_BYM2_deaths_f85 = list(N=N,
                              L=length(neigh),
                              adj=neigh,
                              weights=rep(1,length(neigh)),
                              num=numneigh,
                              scaling_factor=scaling_factor,
                              E=E_f85,
                              X=X_agg,
                              p=p)
Data_BYM2_deaths_f85 = list(y=y_deaths_f85)
Inits_BYM2_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=y_deaths_f85, sigma=1, s=rep(0,N), theta=rep(0,N), lambda=0.5, y=y_deaths_f85)
BYM2_deaths_f85_Model = nimbleModel(code=BYM2_Code, 
                                    constants=Consts_BYM2_deaths_f85, 
                                    data=Data_BYM2_deaths_f85, 
                                    inits=Inits_BYM2_deaths_f85)
BYM2_deaths_f85_conf <- configureMCMC(BYM2_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma", "s", "theta", "lambda"),enableWAIC=T)
BYM2_deaths_f85_MCMC=buildMCMC(BYM2_deaths_f85_conf)
BYM2_deaths_f85_Cmodel=compileNimble(BYM2_deaths_f85_Model)
BYM2_deaths_f85_Cmcmc=compileNimble(BYM2_deaths_f85_MCMC, project = BYM2_deaths_f85_Model)
BYM2_deaths_f85_runMCMC=runMCMC(BYM2_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_BYM2_deaths_f85 <- BYM2_deaths_f85_runMCMC$samples[,c("lambda", "beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                              "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                              "beta[11]","beta[12]","beta[13]")]
# names <- c("lambda", "intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(beta_BYM2_deaths_f85)) {
#   colnames(beta_BYM2_deaths_f85[[i]]) <- names
# }
# summary(beta_BYM2_deaths_f85)

saveRDS(BYM2_deaths_f85_runMCMC, "BYM2_deaths_f85_runMCMC.rds")

#######################
#   NIMBLE: Leroux    #
#######################
Consts_Leroux_deaths_f85 = list(N=N,
                                Q=as.matrix(diag(numneigh) - W),
                                E=E_f85,
                                X=X_agg,
                                mu_b=rep(1,N),
                                p=p)
Data_Leroux_deaths_f85 = list(y=y_deaths_f85)
Inits_Leroux_deaths_f85 = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), sigma=1, b=rep(0,N), lambda=0.5, y=rep(1,N))

Leroux_deaths_f85_Model = nimbleModel(code=Leroux_Code, 
                                      constants=Consts_Leroux_deaths_f85, 
                                      data=Data_Leroux_deaths_f85, 
                                      inits=Inits_Leroux_deaths_f85)
Leroux_deaths_f85_conf <- configureMCMC(Leroux_deaths_f85_Model, monitors = c("beta0", "beta","fitted", "y", #"b.centred", 
                                                                              "b", 
                                                                              "sigma", "lambda"), enableWAIC=T)
Leroux_deaths_f85_MCMC=buildMCMC(Leroux_deaths_f85_conf)
Leroux_deaths_f85_Cmodel=compileNimble(Leroux_deaths_f85_Model)
Leroux_deaths_f85_Cmcmc=compileNimble(Leroux_deaths_f85_MCMC, project = Leroux_deaths_f85_Model)
Leroux_deaths_f85_runMCMC=runMCMC(Leroux_deaths_f85_Cmcmc, nburnin=5e4, niter=1e5, nchains=4, thin=10, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_Leroux_deaths_f85 <- Leroux_deaths_f85_runMCMC$samples[,c("lambda","beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                  "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                  "beta[11]","beta[12]","beta[13]")]
# names <- c("lambda","intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(Leroux_deaths_f85_runMCMC)) {
#   colnames(beta_Leroux_deaths_f85[[i]]) <- names
# }
# summary(beta_Leroux_deaths_f85)
# 
# plot(Leroux_deaths_f85_runMCMC$samples[,c("lambda","beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                    "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                    "beta[11]","beta[12]","beta[13]")])
# plot(Leroux_deaths_f85_runMCMC$samples[,c("lambda")])
# 
# Leroux_deaths_f85_runMCMC$WAIC$WAIC

saveRDS(Leroux_deaths_f85_runMCMC, "Leroux_deaths_f85_runMCMC.rds")