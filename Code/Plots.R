library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstan)
library(coda)
library(loo)

###########################
### Overdispersion
###########################
pdf(file="Plots/hist_overdispersion")
deaths_shp_MA_overall$hist_deaths <- deaths_shp_MA_overall$sum_deaths
deaths_shp_MA_overall$hist_deaths[deaths_shp_MA_overall$hist_deaths>250] <- 250
ggplot(data = data.frame(deaths = deaths_shp_MA_overall$hist_deaths), aes(x = deaths)) +
  geom_histogram(binwidth = 55, color = "black", fill = "white") +
  labs(x = "Deaths", y = "Count", title = "")
dev.off()
mean(deaths_shp_MA_overall$sum_deaths)
var(deaths_shp_MA_overall$sum_deaths)

###########################
### Weight matrix 
###########################
MA.contig_overall <- st_contiguity(deaths_shp_MA_overall)
MA.contig_sf <- st_as_edges(st_centroid(deaths_shp_MA_overall$geometry), nb = MA.contig_overall)
weights.contig.B <- st_weights(MA.contig_overall, style = "B")
W <- wt_as_matrix(MA.contig_overall, weights.contig.B)

pdf("Plots/binmatrix")
# global_moran_test(y_deaths, nb = MA.contig_overall, wt = weights.contig.B)
image(W, axes = FALSE, main = "Binary Contiguous")
dev.off()

pdf("Plots/binneigh")
ggplot() + geom_sf(data = deaths_shp_MA_overall) +
  geom_sf(data = MA.contig_sf) + ggtitle("Contiguous Neighbors")+  theme_void()
dev.off()


# Symmetric K-nearest neighbors
MA.distk4.sym <- st_knn(st_centroid(deaths_shp_MA_overall$geometry), k = 4, symmetric = TRUE)
MA.distk4_sf <- st_as_edges(st_centroid(deaths_shp_MA_overall$geometry), nb = MA.distk4.sym)
pdf("Plots/4nnmneigh")
ggplot(deaths_shp_MA_overall) + geom_sf() + geom_sf(data = MA.distk4_sf) + ggtitle("4th NN Neighbor")+  theme_void()
dev.off()
weights.dist4.B.sym <- st_weights(MA.distk4.sym, style = "B") # binary weights
pdf("Plots/4nnmatrix")
image(wt_as_matrix(nb = MA.distk4.sym, wt = weights.dist4.B.sym),
      axes = FALSE, main = "Binary KNN")
dev.off()
# global_moran_test(deaths_shp_MA_overall$sum_deaths, nb = MA.distk4.sym, wt = weights.dist4.B.sym)

# Distance-based neighbors
MA.dist.range <- st_dist_band(st_centroid(deaths_shp_MA_overall$geometry), 
                              lower = 0, # specify lower bound of range
                              upper = 10 # specify upper bound of range
)
MA.dist.range_sf <- st_as_edges(st_centroid(deaths_shp_MA_overall$geometry), nb = MA.dist.range)
pdf("Plots/distneigh")
ggplot(deaths_shp_MA_overall) + geom_sf() + geom_sf(data = MA.dist.range_sf) + ggtitle("Distance Cutoff Neighbors")+  theme_void()
dev.off()
weights.inv.dist <- st_inverse_distance(nb = MA.dist.range, 
                                        geometry = st_centroid(deaths_shp_MA_overall$geometry)
) # inverse distance weights
pdf("Plots/distmatrix")
image(log(1+wt_as_matrix(nb = MA.dist.range, wt = weights.inv.dist)),
      axes = FALSE, main = "Inverse distance weights")
dev.off()
# global_moran_test(deaths_shp_MA_overall$sum_deaths, nb = MA.dist.range, wt = weights.inv.dist)

# Cumulative neighbors
MA.contig.lag2 <- st_nb_lag_cumul(MA.contig_overall, order = 2)
MA.contig.lag2_sf <- st_as_edges(st_centroid(deaths_shp_MA_overall$geometry), nb = MA.contig.lag2)
pdf("Plots/lagneigh")
ggplot() + geom_sf(data = deaths_shp_MA_overall) +
  geom_sf(data = MA.contig.lag2_sf) + ggtitle("Lagged Contiguous Neighbors")+  theme_void()
dev.off()
pdf("Plots/lagmatrix")
weights.contig.lag2.B <- st_weights(MA.contig.lag2, style = "B") # binary weights
image(wt_as_matrix(nb = MA.contig.lag2, wt = weights.contig.lag2.B),
      axes = FALSE, main = "Binary Lagged")
dev.off()
# global_moran_test(deaths_shp_MA_overall$sum_deaths, nb = MA.contig.lag2, wt = weights.contig.lag2.B)

###########################
### Coefficients stratified
###########################
Simple_m65_samples <- do.call("rbind", readRDS("models/Simple_deaths_m65.rds")$samples)
IID_m65_samples <- do.call("rbind", readRDS("models/IID_deaths_m65.rds")$samples)
CAR_m65_samples <- do.call("rbind", readRDS("models/CAR_deaths_m65.rds")$samples)
BYM_m65_samples <- do.call("rbind", readRDS("models/BYM_deaths_m65.rds")$samples)
BYM2_m65_samples <- do.call("rbind", readRDS("models/BYM2_deaths_m65.rds")$samples)
Leroux_m65_samples <- do.call("rbind", readRDS("models/Leroux_deaths_m65.rds")$samples)
Simple_m75_samples <- do.call("rbind", readRDS("models/Simple_deaths_m75.rds")$samples)
IID_m75_samples <- do.call("rbind", readRDS("models/IID_deaths_m75.rds")$samples)
CAR_m75_samples <- do.call("rbind", readRDS("models/CAR_deaths_m75.rds")$samples)
BYM_m75_samples <- do.call("rbind", readRDS("models/BYM_deaths_m75.rds")$samples)
BYM2_m75_samples <- do.call("rbind", readRDS("models/BYM2_deaths_m75.rds")$samples)
Leroux_m75_samples <- do.call("rbind", readRDS("models/Leroux_deaths_m75.rds")$samples)
Simple_m85_samples <- do.call("rbind", readRDS("models/Simple_deaths_m85.rds")$samples)
IID_m85_samples <- do.call("rbind", readRDS("models/IID_deaths_m85.rds")$samples)
CAR_m85_samples <- do.call("rbind", readRDS("models/CAR_deaths_m85.rds")$samples)
BYM_m85_samples <- do.call("rbind", readRDS("models/BYM_deaths_m85.rds")$samples)
BYM2_m85_samples <- do.call("rbind", readRDS("models/BYM2_deaths_m85.rds")$samples)
Leroux_m85_samples <- do.call("rbind", readRDS("models/Leroux_deaths_m85.rds")$samples)
Simple_f65_samples <- do.call("rbind", readRDS("models/Simple_deaths_f65.rds")$samples)
IID_f65_samples <- do.call("rbind", readRDS("models/IID_deaths_f65.rds")$samples)
CAR_f65_samples <- do.call("rbind", readRDS("models/CAR_deaths_f65.rds")$samples)
BYM_f65_samples <- do.call("rbind", readRDS("models/BYM_deaths_f65.rds")$samples)
BYM2_f65_samples <- do.call("rbind", readRDS("models/BYM2_deaths_f65.rds")$samples)
Leroux_f65_samples <- do.call("rbind", readRDS("models/Leroux_deaths_f65.rds")$samples)
Simple_f75_samples <- do.call("rbind", readRDS("models/Simple_deaths_f75.rds")$samples)
IID_f75_samples <- do.call("rbind", readRDS("models/IID_deaths_f75.rds")$samples)
CAR_f75_samples <- do.call("rbind", readRDS("models/CAR_deaths_f75.rds")$samples)
BYM_f75_samples <- do.call("rbind", readRDS("models/BYM_deaths_f75.rds")$samples)
BYM2_f75_samples <- do.call("rbind", readRDS("models/BYM2_deaths_f75.rds")$samples)
Leroux_f75_samples <- do.call("rbind", readRDS("models/Leroux_deaths_f75.rds")$samples)
Simple_f85_samples <- do.call("rbind", readRDS("models/Simple_deaths_f85.rds")$samples)
IID_f85_samples <- do.call("rbind", readRDS("models/IID_deaths_f85.rds")$samples)
CAR_f85_samples <- do.call("rbind", readRDS("models/CAR_deaths_f85.rds")$samples)
BYM_f85_samples <- do.call("rbind", readRDS("models/BYM_deaths_f85.rds")$samples)
BYM2_f85_samples <- do.call("rbind", readRDS("models/BYM2_deaths_f85.rds")$samples)
Leroux_f85_samples <- do.call("rbind", readRDS("models/Leroux_deaths_f85.rds")$samples)
Simple_overall_samples <- do.call("rbind", readRDS("models/Simple_deaths_overall.rds")$samples)
IID_overall_samples <- do.call("rbind", readRDS("models/IID_deaths_overall.rds")$samples)
CAR_overall_samples <- do.call("rbind", readRDS("models/CAR_deaths_overall.rds")$samples)
BYM_overall_samples <- do.call("rbind", readRDS("models/BYM_deaths_overall.rds")$samples)
BYM2_overall_samples <- do.call("rbind", readRDS("models/BYM2_deaths_overall.rds")$samples)
Leroux_overall_samples <- do.call("rbind", readRDS("models/Leroux_deaths_overall.rds")$samples)
lasso_m65 <- readRDS("models/fitlasso_deaths_m65.rds")
lasso_m75 <- readRDS("models/fitlasso_deaths_m75.rds")
lasso_m85 <- readRDS("models/fitlasso_deaths_m85.rds")
lasso_f65 <- readRDS("models/fitlasso_deaths_f65.rds")
lasso_f75 <- readRDS("models/fitlasso_deaths_f75.rds")
lasso_f85 <- readRDS("models/fitlasso_deaths_f85.rds")
lasso_overall <- readRDS("models/fitlasso_deaths_overall.rds")


mcmc_draws_list <- list(Simple_m65_samples, IID_m65_samples, CAR_m65_samples, BYM_m65_samples, BYM2_m65_samples, Leroux_m65_samples,
                        Simple_f65_samples, IID_f65_samples, CAR_f65_samples, BYM_f65_samples, BYM2_f65_samples, Leroux_f65_samples,
                        Simple_m75_samples, IID_m75_samples, CAR_m75_samples, BYM_m75_samples, BYM2_m75_samples, Leroux_m75_samples,
                        Simple_f75_samples, IID_f75_samples, CAR_f75_samples, BYM_f75_samples, BYM2_f75_samples, Leroux_f75_samples,
                        Simple_m85_samples, IID_m85_samples, CAR_m85_samples, BYM_m85_samples, BYM2_m85_samples, Leroux_m85_samples,
                        Simple_f85_samples, IID_f85_samples, CAR_f85_samples, BYM_f85_samples, BYM2_f85_samples, Leroux_f85_samples)
n_iterations <- dim(mcmc_draws_list[[1]])[1]
data_summary <- data.frame()
model_names <- c("Simple", "IID", "CAR", "BYM", "BYM2", "Leroux")
covariate_names <- c("pm25", "poverty", "popdensity",
                     "medhousevalue", "pct_black", "medhouseincome", "pct_owner_occ",
                     "pct_hispanic", "education", "smoke_rate", "mean_bmi", "amb_visit_pct", "a1c_exm_pct")
for (i in seq_along(mcmc_draws_list)) { #compute summaries
  matrix_data <- mcmc_draws_list[[i]][,paste0("beta[", seq(1,13), "]")]
summary_data <- data.frame(
  model = factor(rep(model_names[(i - 1) %% 6 + 1], 13), levels = model_names),
  group = rep(floor((i - 1) / 6) + 1, 13),
  age_group = rep(c("65-74", "75-84", "85+"), each = 6, times = 2)[i],
  sex = rep(c("male", "Unstructured"), each = 6 * 3)[i],
  covariate = 1:13,
  covariate_name = covariate_names,
  mean_value = colMeans(matrix_data),
  lower_ci = apply(matrix_data, 2, function(col) quantile(col, 0.025)),
  upper_ci = apply(matrix_data, 2, function(col) quantile(col, 0.975))
)
  data_summary <- rbind(data_summary, summary_data)
}

data_summary$age_group <- factor(data_summary$age_group, levels = c("65-74", "75-84", "85+"))
data_summary$sex <- factor(data_summary$sex, levels = c("male", "Unstructured"))
plot <- ggplot(data_summary, aes(x = covariate + 0.1 * as.numeric(model), y = mean_value, color = model)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  facet_grid(age_group ~ sex, labeller = labeller(sex = c("male" = "Male", "Unstructured" = "Unstructured"), age_group = c("65-74" = "65-74", "75-84" = "75-84", "85+" = "85+"))) +
  scale_x_continuous(
    breaks = seq(1, 13, by = 1),
    labels = covariate_names,
    name = ""
  ) +
  scale_color_discrete(name = "Model", labels = model_names) +
  labs(
    # title = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pdf(file="Plots/coefficients_stratified.pdf")
print(plot)
dev.off()

###########################
### Coefficients overall
###########################
mcmc_draws_list <-mcmc_draws_list <- list(Simple_overall_samples, IID_overall_samples, CAR_overall_samples, BYM_overall_samples, BYM2_overall_samples, Leroux_overall_samples)
n_iterations <- dim(mcmc_draws_list[[1]])[1]
data_summary <- data.frame()
model_names <- c("Simple", "IID", "CAR", "BYM", "BYM2", "Leroux")
for (i in seq_along(mcmc_draws_list)) { #compute summaries
  matrix_data <- mcmc_draws_list[[i]][,paste0("beta[", seq(1,13), "]")]
  summary_data <- data.frame(
    model = factor(rep(model_names[i], 13), levels = model_names),
    covariate = 1:13,
    mean_value = colMeans(matrix_data),
    lower_ci = apply(matrix_data, 2, function(x) quantile(x, 0.025)),
    upper_ci = apply(matrix_data, 2, function(x) quantile(x, 0.975))
  )
    data_summary <- rbind(data_summary, summary_data)
}
plot <- ggplot(data_summary, aes(x = covariate + 0.1 * as.numeric(model), y = mean_value, color = model)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  scale_x_continuous(
    breaks = seq(1, 13, by = 1),
    labels = covariate_names,
    name = ""
  ) +
  scale_color_discrete(name = "Model", labels = model_names) +
  labs(
    # title = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(plot)
pdf(file="Plots/coefficients_overall.pdf")
print(plot)
dev.off()

###########################
### Coefficients lasso stratified
###########################
lasso_m65_samples <- extract(lasso_m65)$beta
lasso_m75_samples <- extract(lasso_m75)$beta
lasso_m85_samples <- extract(lasso_m85)$beta
lasso_f65_samples <- extract(lasso_f65)$beta
lasso_f75_samples <- extract(lasso_f75)$beta
lasso_f85_samples <- extract(lasso_f85)$beta
mcmc_draws_list <- list(lasso_m65_samples, lasso_m75_samples, lasso_m85_samples, lasso_f65_samples, lasso_f75_samples, lasso_f85_samples)
n_iterations <- dim(mcmc_draws_list[[1]])[1]
data_summary <- data.frame()
group_names <- c("male_65-74", "male_75-84", "male_85+", "Unstructured_65-74", "Unstructured_75-84", "Unstructured_85+")
for (i in seq_along(mcmc_draws_list)) { #compute summaries
  matrix_data <- mcmc_draws_list[[i]]
  summary_data <- data.frame(
    group = factor(rep(group_names[i], 13), levels = group_names),
    covariate = 1:13,
    mean_value = colMeans(matrix_data),
    lower_ci = apply(matrix_data, 2, function(x) quantile(x, 0.025)),
    upper_ci = apply(matrix_data, 2, function(x) quantile(x, 0.975))
  )
    data_summary <- rbind(data_summary, summary_data)
}
data_summary <- data_summary %>%
  separate(group, c("sex", "age_group"), sep = "_", remove = FALSE)
data_summary$age_group <- factor(data_summary$age_group, levels = c("65-74", "75-84", "85+"))
data_summary$sex <- factor(data_summary$sex, levels = c("male", "Unstructured"))
plot <- ggplot(data_summary, aes(x = covariate, y = mean_value)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  facet_grid(age_group ~ sex, labeller = labeller(sex = c("male" = "Male", "Unstructured" = "Unstructured"), age_group = c("65-74" = "65-74", "75-84" = "75-84", "85+" = "85+"))) +
  scale_x_continuous(
    breaks = seq(1, 13, by = 1),
    labels = covariate_names,
    name = ""
  ) +
  labs(
    y = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pdf(file="Plots/coefficients_lasso_stratified.pdf")
print(plot)
dev.off()

###########################
### Coefficients overall
###########################
matrix_data <- extract(lasso_overall_samples)$beta
n_iterations <- dim(matrix_data)[1]
data_summary <- data.frame( #compute summaries
  covariate = 1:13,
  mean_value = colMeans(matrix_data),
  lower_ci = apply(matrix_data, 2, function(x) quantile(x, 0.025)),
  upper_ci = apply(matrix_data, 2, function(x) quantile(x, 0.975))
)
plot <- ggplot(data_summary, aes(x = covariate, y = mean_value)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0) +
  scale_x_continuous(
    breaks = seq(1, 13, by = 1),
    labels = covariate_names,
    name = ""
  ) +
  labs(
    y = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pdf(file="Plots/coefficients_lasso_overall.pdf")
print(plot)
dev.off()

###########################
### Maps of PM25 and covariates
###########################

pdf(file="Plots/pm25")
deaths_shp_MA_overall$SMR <- deaths_shp_MA_overall$sum_deaths/deaths_shp_MA_overall$people

ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = pm25_ensemble))+
  scale_fill_gradientn(name = "PM2.5",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/SMR")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = SMR))+
  scale_fill_gradientn(name = "SMR",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()


pdf(file="Plots/poverty")
names
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = poverty))+
  scale_fill_gradientn(name = "poverty",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/popdensity")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = popdensity))+
  scale_fill_gradientn(name = "popdensity",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/medianhousevalue")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = medianhousevalue))+
  scale_fill_gradientn(name = "medianhousevalue",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/pct_blk")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = pct_blk))+
  scale_fill_gradientn(name = "pct_blk",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/medhouseholdincome")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = medhouseholdincome))+
  scale_fill_gradientn(name = "medhouseholdincome",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/pct_owner_occ")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = pct_owner_occ))+
  scale_fill_gradientn(name = "pct_owner_occ",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/pct_owner_occ")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = hispanic))+
  scale_fill_gradientn(name = "hispanic",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/education")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = education))+
  scale_fill_gradientn(name = "education",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/smoke_rate")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = smoke_rate))+
  scale_fill_gradientn(name = "smoke_rate",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/mean_bmi")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = mean_bmi))+
  scale_fill_gradientn(name = "mean_bmi",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/amb_visit_pct")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = amb_visit_pct))+
  scale_fill_gradientn(name = "amb_visit_pct",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

pdf(file="Plots/a1c_exm_pct")
ggplot() + 
  geom_sf(data = deaths_shp_MA_overall, aes(fill = a1c_exm_pct))+
  scale_fill_gradientn(name = "a1c_exm_pct",
                       colours = colors,
                       na.value = "grey50")+
  theme_void()
dev.off()

###########################
### Map relative risk
###########################
fitted_Simple = Simple_overall_samples[,paste0("fitted[", seq(1,N), "]")]
fitted_IID = IID_overall_samples[,paste0("fitted[", seq(1,N), "]")] 
fitted_CAR = CAR_overall_samples[,paste0("fitted[", seq(1,N), "]")]
fitted_BYM = BYM_overall_samples[,paste0("fitted[", seq(1,N), "]")]
fitted_BYM2 = BYM2_overall_samples[,paste0("fitted[", seq(1,N), "]")]
fitted_Leroux = Leroux_overall_samples[,paste0("fitted[", seq(1,N), "]")]
RR = rbind(apply(fitted_Simple,2,mean)/E, apply(fitted_IID,2,mean)/E, apply(fitted_CAR,2,mean)/E, apply(fitted_BYM,2,mean)/E, apply(fitted_BYM2,2,mean)/E, apply(fitted_Leroux,2,mean)/E)
deaths_shp_MA_overall_RR <- deaths_shp_MA_overall
deaths_shp_MA_overall_RR$simple <- RR[1,]
deaths_shp_MA_overall_RR$IID <- RR[2,]
deaths_shp_MA_overall_RR$CAR <- RR[3,]
deaths_shp_MA_overall_RR$BYM <- RR[4,]
deaths_shp_MA_overall_RR$BYM2 <- RR[5,]
deaths_shp_MA_overall_RR$leroux <- RR[6,]
  
pdf(file="Plots/simpleRR")
colors <- c("#045a8d", "#7fb7d2", "#b7d5e5", "#ffffbf", "#f46d43", "#d73027", "#a50026")
simpleRR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = simple))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("Simple")+
theme_void()
dev.off()

pdf(file="Plots/IIDRR")
IIDRR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = IID))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("IID")+
theme_void()
dev.off()

pdf(file="Plots/CARRR")
CARRR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = CAR))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("CAR")+
theme_void()
dev.off()

pdf(file="Plots/BYMRR")
BYMRR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = BYM))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("BYM")+
theme_void()
dev.off()

pdf(file="Plots/BYM2RR")
BYM2RR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = BYM2))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("BYM2")+
theme_void()
dev.off()

pdf(file="Plots/LerouxRR")
LerouxRR <- ggplot() + 
geom_sf(data = deaths_shp_MA_overall_RR, aes(fill = leroux))+
scale_fill_gradientn(name = "RR",
                 colours = colors,
                 na.value = "grey50")+
ggtitle("Leroux")+
theme_void()
dev.off()

pdf(file="Plots/RR")
ggarrange(simpleRR,IIDRR,CARRR, BYMRR, BYM2RR, LerouxRR, ncol=2, nrow=3)
dev.off()

###########################
### WAIC
###########################
fi_bar_Simple= apply(extract(lasso_overall)$lik, 2,mean) ; p_w_Simple = sum(apply(log(extract(lasso_overall)$lik), 2,var))
waic_Simple = -2*sum(log(fi_bar_Simple)) + 2*p_w_Simple #compute by hand

paste(round(c(Simple_overall$WAIC$WAIC,Simple_overall$WAIC$pWAIC,IID_overall$WAIC$WAIC,IID_overall$WAIC$pWAIC,
              CAR_overall$WAIC$WAIC,CAR_overall$WAIC$pWAIC,BYM_overall$WAIC$WAIC,BYM_overall$WAIC$pWAIC,
              BYM2_overall$WAIC$WAIC,BYM2_overall$WAIC$pWAIC,Leroux_overall$WAIC$WAIC,Leroux_overall$WAIC$pWAIC,
              waic(extract_log_lik(lasso_overall))$waic,waic(extract_log_lik(lasso_overall))$p_waic),1), collapse = " & ")

paste(round(c(Simple_m65$WAIC$WAIC,Simple_m65$WAIC$pWAIC,IID_m65$WAIC$WAIC,IID_m65$WAIC$pWAIC,
              CAR_m65$WAIC$WAIC,CAR_m65$WAIC$pWAIC,BYM_m65$WAIC$WAIC,BYM_m65$WAIC$pWAIC,
              BYM2_m65$WAIC$WAIC,BYM2_m65$WAIC$pWAIC,Leroux_m65$WAIC$WAIC,Leroux_m65$WAIC$pWAIC,
              waic(extract_log_lik(lasso_m65))$waic,waic(extract_log_lik(lasso_m65))$p_waic),1), collapse = " & ")

paste(round(c(Simple_f65$WAIC$WAIC,Simple_f65$WAIC$pWAIC,IID_f65$WAIC$WAIC,IID_f65$WAIC$pWAIC,
              CAR_f65$WAIC$WAIC,CAR_f65$WAIC$pWAIC,BYM_f65$WAIC$WAIC,BYM_f65$WAIC$pWAIC,
              BYM2_f65$WAIC$WAIC,BYM2_f65$WAIC$pWAIC,Leroux_f65$WAIC$WAIC,Leroux_f65$WAIC$pWAIC,
              waic(extract_log_lik(lasso_f65))$waic,waic(extract_log_lik(lasso_f65))$p_waic),1), collapse = " & ")

paste(round(c(Simple_m75$WAIC$WAIC,Simple_m75$WAIC$pWAIC,IID_m75$WAIC$WAIC,IID_m75$WAIC$pWAIC,
              CAR_m75$WAIC$WAIC,CAR_m75$WAIC$pWAIC,BYM_m75$WAIC$WAIC,BYM_m75$WAIC$pWAIC,
              BYM2_m75$WAIC$WAIC,BYM2_m75$WAIC$pWAIC,Leroux_m75$WAIC$WAIC,Leroux_m75$WAIC$pWAIC,
              waic(extract_log_lik(lasso_m75))$waic,waic(extract_log_lik(lasso_m75))$p_waic),1), collapse = " & ")

paste(round(c(Simple_f75$WAIC$WAIC,Simple_f75$WAIC$pWAIC,IID_f75$WAIC$WAIC,IID_f75$WAIC$pWAIC,
              CAR_f75$WAIC$WAIC,CAR_f75$WAIC$pWAIC,BYM_f75$WAIC$WAIC,BYM_f75$WAIC$pWAIC,
              BYM2_f75$WAIC$WAIC,BYM2_f75$WAIC$pWAIC,Leroux_f75$WAIC$WAIC,Leroux_f75$WAIC$pWAIC,
              waic(extract_log_lik(lasso_f75))$waic,waic(extract_log_lik(lasso_f75))$p_waic),1), collapse = " & ")

paste(round(c(Simple_m85$WAIC$WAIC,Simple_m85$WAIC$pWAIC,IID_m85$WAIC$WAIC,IID_m85$WAIC$pWAIC,
              CAR_m85$WAIC$WAIC,CAR_m85$WAIC$pWAIC,BYM_m85$WAIC$WAIC,BYM_m85$WAIC$pWAIC,
              BYM2_m85$WAIC$WAIC,BYM2_m85$WAIC$pWAIC,Leroux_m85$WAIC$WAIC,Leroux_m85$WAIC$pWAIC,
              waic(extract_log_lik(lasso_m85))$waic,waic(extract_log_lik(lasso_m85))$p_waic),1), collapse = " & ")

paste(round(c(Simple_f85$WAIC$WAIC,Simple_f85$WAIC$pWAIC,IID_f85$WAIC$WAIC,IID_f85$WAIC$pWAIC,
              CAR_f85$WAIC$WAIC,CAR_f85$WAIC$pWAIC,BYM_f85$WAIC$WAIC,BYM_f85$WAIC$pWAIC,
              BYM2_f85$WAIC$WAIC,BYM2_f85$WAIC$pWAIC,Leroux_f85$WAIC$WAIC,Leroux_f85$WAIC$pWAIC,
              waic(extract_log_lik(lasso_f85))$waic,waic(extract_log_lik(lasso_f85))$p_waic),1), collapse = " & ")

###########################
### Convergence
###########################
Simple_m65 <- readRDS("models/Simple_deaths_m65.rds")
IID_m65 <- readRDS("models/IID_deaths_m65.rds")
CAR_m65 <- readRDS("models/CAR_deaths_m65.rds")
BYM_m65 <- readRDS("models/BYM_deaths_m65.rds")
BYM2_m65 <- readRDS("models/BYM2_deaths_m65.rds")
Leroux_m65 <- readRDS("models/Leroux_deaths_m65.rds")
Simple_m75 <- readRDS("models/Simple_deaths_m75.rds")
IID_m75 <- readRDS("models/IID_deaths_m75.rds")
CAR_m75 <- readRDS("models/CAR_deaths_m75.rds")
BYM_m75 <- readRDS("models/BYM_deaths_m75.rds")
BYM2_m75 <- readRDS("models/BYM2_deaths_m75.rds")
Leroux_m75 <- readRDS("models/Leroux_deaths_m75.rds")
Simple_m85 <- readRDS("models/Simple_deaths_m85.rds")
IID_m85 <- readRDS("models/IID_deaths_m85.rds")
CAR_m85 <- readRDS("models/CAR_deaths_m85.rds")
BYM_m85 <- readRDS("models/BYM_deaths_m85.rds")
BYM2_m85 <- readRDS("models/BYM2_deaths_m85.rds")
Leroux_m85 <- readRDS("models/Leroux_deaths_m85.rds")
Simple_f65 <- readRDS("models/Simple_deaths_f65.rds")
IID_f65 <- readRDS("models/IID_deaths_f65.rds")
CAR_f65 <- readRDS("models/CAR_deaths_f65.rds")
BYM_f65 <- readRDS("models/BYM_deaths_f65.rds")
BYM2_f65 <- readRDS("models/BYM2_deaths_f65.rds")
Leroux_f65 <- readRDS("models/Leroux_deaths_f65.rds")
Simple_f75 <- readRDS("models/Simple_deaths_f75.rds")
IID_f75 <- readRDS("models/IID_deaths_f75.rds")
CAR_f75 <- readRDS("models/CAR_deaths_f75.rds")
BYM_f75 <- readRDS("models/BYM_deaths_f75.rds")
BYM2_f75 <- readRDS("models/BYM2_deaths_f75.rds")
Leroux_f75 <- readRDS("models/Leroux_deaths_f75.rds")
Simple_f85 <- readRDS("models/Simple_deaths_f85.rds")
IID_f85 <- readRDS("models/IID_deaths_f85.rds")
CAR_f85 <- readRDS("models/CAR_deaths_f85.rds")
BYM_f85 <- readRDS("models/BYM_deaths_f85.rds")
BYM2_f85 <- readRDS("models/BYM2_deaths_f85.rds")
Leroux_f85 <- readRDS("models/Leroux_deaths_f85.rds")
Simple_overall <- readRDS("models/Simple_deaths_overall.rds")
IID_overall <- readRDS("models/IID_deaths_overall.rds")
CAR_overall <- readRDS("models/CAR_deaths_overall.rds")
BYM_overall <- readRDS("models/BYM_deaths_overall.rds")
BYM2_overall <- readRDS("models/BYM2_deaths_overall.rds")
Leroux_overall <- readRDS("models/Leroux_deaths_overall.rds")

pdf(file="Plots/convergence_nimble_BYM2s1")
plot(BYM2_m65$samples[,c("s[1]")])
dev.off()


pdf(file="Plots/convergence_nimble_CARb12")
plot(BYM2_m65$samples[,c("beta[12]")])
dev.off()

Simple_overall_rhat <- gelman.diag(Simple_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_overall_rhat <- gelman.diag(IID_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_overall_rhat <- gelman.diag(CAR_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_overall_rhat <- gelman.diag(BYM_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_overall_rhat <- gelman.diag(BYM2_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_overall_rhat <- gelman.diag(Leroux_overall$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_m65_rhat <- gelman.diag(Simple_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_m65_rhat <- gelman.diag(IID_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_m65_rhat <- gelman.diag(CAR_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_m65_rhat <- gelman.diag(BYM_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_m65_rhat <- gelman.diag(BYM2_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_m65_rhat <- gelman.diag(Leroux_m65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_m65$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_m65$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_m65$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_m65$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_m65$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_m65$samples[,paste0("beta[", seq(1,13), "]")])

Simple_m75_rhat <- gelman.diag(Simple_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_m75_rhat <- gelman.diag(IID_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_m75_rhat <- gelman.diag(CAR_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_m75_rhat <- gelman.diag(BYM_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_m75_rhat <- gelman.diag(BYM2_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_m75_rhat <- gelman.diag(Leroux_m75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_m75$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_m75$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_m75$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_m75$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_m75$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_m75$samples[,paste0("beta[", seq(1,13), "]")])

Simple_m85_rhat <- gelman.diag(Simple_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_m85_rhat <- gelman.diag(IID_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_m85_rhat <- gelman.diag(CAR_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_m85_rhat <- gelman.diag(BYM_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_m85_rhat <- gelman.diag(BYM2_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_m85_rhat <- gelman.diag(Leroux_m85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_m85$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_m85$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_m85$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_m85$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_m85$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_m85$samples[,paste0("beta[", seq(1,13), "]")])

Simple_f65_rhat <- gelman.diag(Simple_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_f65_rhat <- gelman.diag(IID_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_f65_rhat <- gelman.diag(CAR_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_f65_rhat <- gelman.diag(BYM_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_f65_rhat <- gelman.diag(BYM2_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_f65_rhat <- gelman.diag(Leroux_f65$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_f65$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_f65$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_f65$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_f65$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_f65$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_f65$samples[,paste0("beta[", seq(1,13), "]")])

Simple_f75_rhat <- gelman.diag(Simple_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_f75_rhat <- gelman.diag(IID_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_f75_rhat <- gelman.diag(CAR_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_f75_rhat <- gelman.diag(BYM_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_f75_rhat <- gelman.diag(BYM2_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_f75_rhat <- gelman.diag(Leroux_f75$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_f75$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_f75$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_f75$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_f75$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_f75$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_f75$samples[,paste0("beta[", seq(1,13), "]")])

Simple_f85_rhat <- gelman.diag(Simple_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
IID_f85_rhat <- gelman.diag(IID_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
CAR_f85_rhat <- gelman.diag(CAR_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM_f85_rhat <- gelman.diag(BYM_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
BYM2_f85_rhat <- gelman.diag(BYM2_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]
Leroux_f85_rhat <- gelman.diag(Leroux_f85$samples[,paste0("beta[", seq(1,13), "]")])$psrf[,1]

Simple_overall_neff <- effectiveSize(Simple_f85$samples[,paste0("beta[", seq(1,13), "]")])
IID_overall_neff <- effectiveSize(IID_f85$samples[,paste0("beta[", seq(1,13), "]")])
CAR_overall_neff <- effectiveSize(CAR_f85$samples[,paste0("beta[", seq(1,13), "]")])
BYM_overall_neff <- effectiveSize(BYM_f85$samples[,paste0("beta[", seq(1,13), "]")])
BYM2_overall_neff <- effectiveSize(BYM2_f85$samples[,paste0("beta[", seq(1,13), "]")])
Leroux_overall_neff <- effectiveSize(Leroux_f85$samples[,paste0("beta[", seq(1,13), "]")])

summary(lasso_overall, par=c("beta"))
summary(lass_m65, par=c("beta"))
summary(lass_m75, par=c("beta"))
summary(lass_m85, par=c("beta"))
summary(lass_f65, par=c("beta"))
summary(lass_f75, par=c("beta"))
summary(lass_f85, par=c("beta"))