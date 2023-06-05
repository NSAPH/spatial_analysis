library(nimble)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###########################
#        SIMPLE Stan     #
###########################
fitSimple_deaths_overall = stan("Simple.stan",
                            data=list(N=N,y=y_deaths,E=E, X=X_overall, p=p),
                            warmup=1e3, iter=1e4, chains=4, thin=1)
# traceplot(fitSimple_deaths_overall, pars = c("beta"))
# summary(fitSimple_deaths_overall, pars = c("beta"))
# plot(fitSimple_deaths_overall, pars=c("beta"))
saveRDS(fitSimple_deaths_overall, "fitSimple_deaths_overall.rds")

###########################
#        IID Stan        #
###########################
fitIID_deaths_overall = stan("IID.stan",
                         data=list(N=N,y=y_deaths,E=E, X=X_overall, p=p),
                         warmup=1e3, iter=1e4, chains=4, thin=1)

# traceplot(fitIID_deaths_overall, pars = c("beta"))
# summary(fitSimple_deaths_overall, pars = c("beta"))
# plot(fitIID_deaths_overall, pars = c("beta"))
saveRDS(fitIID_deaths_overall, "fitIID_deaths_overall.rds")


###########################
#        CAR Stan        #
###########################
fitCAR_deaths_overall <- stan('CAR.stan', data = list(n = N,        
                                                         p = p,       
                                                         X=X_overall,              
                                                         y = y_deaths,             
                                                         E=E, 
                                                         W_n = sum(W) / 2,   
                                                         W = W),
                          warmup=1e3, iter=1e4, chains=4)
# traceplot(fitCAR_deaths_overall, par = c("beta", "tau"))
# plot(fitCAR_deaths_overall, par = c('beta'))
saveRDS(fitCAR_deaths_overall, "fitCAR_deaths_overall.rds")

###########################
#        BYM stan         #
###########################
fitBYM_deaths_overall = stan("BYM.stan",
                         data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y_deaths,p=p,X=X_overall,E=E),
                         control=list(adapt_delta = 0.97, stepsize = 0.1),
                         warmup=1e3, chains=4, iter=1e4);
# traceplot(fitBYM_deaths_overall, par = c("beta"))
# summary(fitBYM_deaths_overall, par = c("beta"))
# plot(fitBYM_deaths_overall, par = c("beta"))
# summary(fitBYM_deaths_overall, par=c("beta"))
saveRDS(fitBYM_deaths_overall, "fitBYM_deaths_overall.rds")

###########################
#        BYM2 stan        #
###########################
fitBYM2_deaths_overall = stan("BYM2.stan",
                          data=list(N=N,N_edges=N_edges,node1=node1,node2=node2,y=y_deaths,X=X_overall,E=E,scaling_factor=scaling_factor,K=p),
                          warmup=5e3, iter=1e4, chains=4, thin=1);
# traceplot(fitBYM2_deaths_overall, pars=c("phi[1]", "theta[1]"))
# plot(fitBYM2_deaths_overall, par=c("beta"))
# summary(fitBYM2_deaths_overall, par=c("beta", "rho"))
saveRDS(fitBYM2_deaths_overall, "fitBYM2_deaths_overall.rds")

###########################
#        Leroux stan      #
###########################
fitLeroux_deaths_overall = stan("Leroux.stan",
                            data=list(N=N,y=y_deaths,E=E, X=X_overall, p=p, N_edges=N_edges, node1=node1, node2=node2),
                            warmup=1e3, iter=1e4, chains=4, thin=1)
traceplot(fitLeroux_deaths_overall, pars = c("beta", "lambda"))
summary(fitLeroux_deaths_overall, pars = c("beta", "lambda"))
plot(fitLeroux_deaths_overall, pars = c("beta", "lambda"))
saveRDS(fitLeroux_deaths_overall, "fitLeroux_deaths_overall.rds")

#######################
#    NIMBLE: Simple   #
#######################
Consts_Simple_deaths_overall = list(N=N,
                                    E=E,
                                    X=X_overall,
                                    p=p)
Data_Simple_deaths_overall = list(y=y_deaths)
Inits_Simple_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), y=rep(1,N))
Simple_deaths_overall_Model = nimbleModel(code=Simple_Code, 
                                          constants=Consts_Simple_deaths_overall, 
                                          data=Data_Simple_deaths_overall, 
                                          inits=Inits_Simple_deaths_overall)
Simple_deaths_overall_conf <- configureMCMC(Simple_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y"), enableWAIC=T)
Simple_deaths_overall_MCMC=buildMCMC(Simple_deaths_overall_conf)
Simple_deaths_overall_Cmodel=compileNimble(Simple_deaths_overall_Model)
Simple_deaths_overall_Cmcmc=compileNimble(Simple_deaths_overall_MCMC, project = Simple_deaths_overall_Model)
Simple_deaths_overall_runMCMC=runMCMC(Simple_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_Simple_deaths_overall <- Simple_deaths_overall_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                                "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                                "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(Simple_deaths_overall_runMCMC$samples)) {
#   colnames(beta_Simple_deaths_overall[[i]]) <- names
# }
# summary(beta_Simple_deaths_overall)

saveRDS(Simple_deaths_overall_runMCMC, "Simple_deaths_overall_runMCMC.rds")

#######################
#    NIMBLE: IID   #
#######################
Consts_IID_deaths_overall = list(N=N,
                                 E=E,
                                 X=X_overall,
                                 p=p)
Data_IID_deaths_overall = list(y=y_deaths)
Inits_IID_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), sigma_b=1, b=rep(0,N), y=rep(1,N))

IID_deaths_overall_Model = nimbleModel(code=IID_Code, 
                                       constants=Consts_IID_deaths_overall, 
                                       data=Data_IID_deaths_overall, 
                                       inits=Inits_IID_deaths_overall)
IID_deaths_overall_conf <- configureMCMC(IID_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_b"), enableWAIC=T)
IID_deaths_overall_MCMC=buildMCMC(IID_deaths_overall_conf)
IID_deaths_overall_Cmodel=compileNimble(IID_deaths_overall_Model)
IID_deaths_overall_Cmcmc=compileNimble(IID_deaths_overall_MCMC, project = IID_deaths_overall_Model)
IID_deaths_overall_runMCMC=runMCMC(IID_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_IID_deaths_overall <- IID_deaths_overall_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                          "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                          "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(IID_deaths_overall_runMCMC$samples)) {
#   colnames(beta_IID_deaths_overall[[i]]) <- names
# }
# summary(beta_IID_deaths_overall)

saveRDS(IID_deaths_overall_runMCMC, "IID_deaths_overall_runMCMC.rds")

#######################
#     NIMBLE: CAR     #
#######################
Consts_CAR_deaths_overall = list(N=N,
                                 L=length(neigh),
                                 adj=neigh,
                                 weights=rep(1,length(neigh)),
                                 num=numneigh,
                                 E=E,
                                 X=X_overall,
                                 p=p)
Data_CAR_deaths_overall = list(y=y_deaths)
Inits_CAR_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=y_deaths, sigma_b=1, b=rep(0,N), y=y_deaths)
CAR_deaths_overall_Model = nimbleModel(code=CAR_Code, 
                                       constants=Consts_CAR_deaths_overall, 
                                       data=Data_CAR_deaths_overall, 
                                       inits=Inits_CAR_deaths_overall)
CAR_deaths_overall_conf <- configureMCMC(CAR_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_b"), enableWAIC=T)
CAR_deaths_overall_MCMC=buildMCMC(CAR_deaths_overall_conf)
CAR_deaths_overall_Cmodel=compileNimble(CAR_deaths_overall_Model)
CAR_deaths_overall_Cmcmc=compileNimble(CAR_deaths_overall_MCMC, project = CAR_deaths_overall_Model)
CAR_deaths_overall_runMCMC=runMCMC(CAR_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_CAR_deaths_overall <- CAR_deaths_overall_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                    "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                    "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(CAR_deaths_overall_runMCMC$samples)) {
#   colnames(beta_CAR_deaths_overall[[i]]) <- names
# }
# summary(beta_CAR_deaths_overall)

saveRDS(CAR_deaths_overall_runMCMC, "CAR_deaths_overall_runMCMC.rds")

#######################
#     NIMBLE: BYM     #
#######################
Consts_BYM_deaths_overall = list(N=N,
                                 L=length(neigh),
                                 adj=neigh,
                                 weights=rep(1,length(neigh)),
                                 num=numneigh,
                                 E=E,
                                 X=X_overall,
                                 p=p)
Data_BYM_deaths_overall = list(y=y_deaths)
Inits_BYM_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=y_deaths, sigma_s=1, s=rep(0,N),sigma_theta=1, theta=rep(0,N),y=y_deaths)
BYM_deaths_overall_Model = nimbleModel(code=BYM_Code, 
                                       constants=Consts_BYM_deaths_overall, 
                                       data=Data_BYM_deaths_overall, 
                                       inits=Inits_BYM_deaths_overall)
BYM_deaths_overall_conf <- configureMCMC(BYM_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma_s", "sigma_theta", "s", "theta"), enableWAIC=T)
BYM_deaths_overall_MCMC=buildMCMC(BYM_deaths_overall_conf)
BYM_deaths_overall_Cmodel=compileNimble(BYM_deaths_overall_Model)
BYM_deaths_overall_Cmcmc=compileNimble(BYM_deaths_overall_MCMC, project = BYM_deaths_overall_Model)
BYM_deaths_overall_runMCMC=runMCMC(BYM_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_BYM_deaths_overall <- BYM_deaths_overall_runMCMC$samples[,c("beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                          "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                          "beta[11]","beta[12]","beta[13]")]
# names <- c("intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(BYM_deaths_overall_runMCMC$samples)) {
#   colnames(beta_BYM_deaths_overall[[i]]) <- names
# }
# summary(beta_BYM_deaths_overall)

saveRDS(BYM_deaths_overall_runMCMC, "BYM_deaths_overall_runMCMC.rds")

#######################
#     NIMBLE: BYM2    #
#######################
Consts_BYM2_deaths_overall = list(N=N,
                                  L=length(neigh),
                                  adj=neigh,
                                  weights=rep(1,length(neigh)),
                                  num=numneigh,
                                  scaling_factor=scaling_factor,
                                  E=E,
                                  X=X_overall,
                                  p=p)
Data_BYM2_deaths_overall = list(y=y_deaths)
Inits_BYM2_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=y_deaths, sigma=1, s=rep(0,N), theta=rep(0,N), lambda=0.5, y=y_deaths)
BYM2_deaths_overall_Model = nimbleModel(code=BYM2_Code, 
                                        constants=Consts_BYM2_deaths_overall, 
                                        data=Data_BYM2_deaths_overall, 
                                        inits=Inits_BYM2_deaths_overall)
BYM2_deaths_overall_conf <- configureMCMC(BYM2_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y", "b", "sigma", "s", "theta", "lambda"),enableWAIC=T)
BYM2_deaths_overall_MCMC=buildMCMC(BYM2_deaths_overall_conf)
BYM2_deaths_overall_Cmodel=compileNimble(BYM2_deaths_overall_Model)
BYM2_deaths_overall_Cmcmc=compileNimble(BYM2_deaths_overall_MCMC, project = BYM2_deaths_overall_Model)
BYM2_deaths_overall_runMCMC=runMCMC(BYM2_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, thin=10, nchains=4, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_BYM2_deaths_overall <- BYM2_deaths_overall_runMCMC$samples[,c("lambda", "beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                            "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                            "beta[11]","beta[12]","beta[13]")]
# names <- c("lambda", "intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(beta_BYM2_deaths_overall)) {
#   colnames(beta_BYM2_deaths_overall$samples[[i]]) <- names
# }
# summary(beta_BYM2_deaths_overall)
# 
# plot(BYM2_deaths_overall_runMCMC$samples[,c("s[1]", "s[2]")])

saveRDS(BYM2_deaths_overall_runMCMC, "BYM2_deaths_overall_runMCMC.rds")


#######################
#   NIMBLE: Leroux    #
#######################
Consts_Leroux_deaths_overall = list(N=N,
                                    Q=as.matrix(diag(numneigh) - W),
                                    E=E,
                                    X=X_overall,
                                    mu_b=rep(1,N),
                                    p=p)
Data_Leroux_deaths_overall = list(y=y_deaths)
Inits_Leroux_deaths_overall = list(beta0=0, beta=rep(0,p), fitted=rep(0,N), sigma=1, b=rep(0,N), lambda=0.5, y=rep(1,N))

Leroux_deaths_overall_Model = nimbleModel(code=Leroux_Code, 
                                          constants=Consts_Leroux_deaths_overall, 
                                          data=Data_Leroux_deaths_overall, 
                                          inits=Inits_Leroux_deaths_overall)
Leroux_deaths_overall_conf <- configureMCMC(Leroux_deaths_overall_Model, monitors = c("beta0", "beta","fitted", "y", #"b.centred", 
                                                                                      "b", 
                                                                                      "sigma", "lambda"), enableWAIC=T)
Leroux_deaths_overall_MCMC=buildMCMC(Leroux_deaths_overall_conf)
Leroux_deaths_overall_Cmodel=compileNimble(Leroux_deaths_overall_Model)
Leroux_deaths_overall_Cmcmc=compileNimble(Leroux_deaths_overall_MCMC, project = Leroux_deaths_overall_Model)
Leroux_deaths_overall_runMCMC=runMCMC(Leroux_deaths_overall_Cmcmc, nburnin=5e4, niter=1e5, nchains=4, thin=10, samplesAsCodaMCMC = TRUE, WAIC=T)

# beta_Leroux_deaths_overall <- Leroux_deaths_overall_runMCMC$samples[,c("lambda","beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                                                "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                                                "beta[11]","beta[12]","beta[13]")]
# names <- c("lambda","intercept", "pm25_ensemble","poverty","popdensity",
#            "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
#            "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
# for (i in 1:length(Leroux_deaths_overall_runMCMC$samples)) {
#   colnames(beta_Leroux_deaths_overall[[i]]) <- names
# }
# summary(beta_Leroux_deaths_overall)
# 
# plot(Leroux_deaths_overall_runMCMC$samples[,c("lambda","beta0", "beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
#                                           "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
#                                           "beta[11]","beta[12]","beta[13]")])
# plot(Leroux_deaths_overall_runMCMC$samples[,c("lambda")])

saveRDS(Leroux_deaths_overall_runMCMC, "Leroux_deaths_overall_runMCMC.rds")