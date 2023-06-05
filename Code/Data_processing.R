###
# Run script below first to create shapefile "deaths_shp_MA_agg.rds" and "deaths_shp_MA_overall.rds", then only use these shapefiles processed in the Model_data.R files
###

library(dplyr)
library(lubridate)
library(fst)
library(ggplot2)
library(data.table)
library(sfdep)
library(sf)
library(sp)
library(spdep)
library(spatialreg)
library(MASS)
library(gridExtra)
library(sfdep)
library(mapview)

### deaths data set 2016
deaths <- read.fst("/n/dominici_nsaph_l3/Lab/projects/analytic/denom_by_year/confounder_exposure_merged_nodups_health_2016.fst")

# Filter data by state
deaths_MA <- filter(deaths, statecode == "MA")
deaths_MA_m65 <- filter(deaths_MA, age >= 65 & age < 75, sex == 1)
deaths_MA_m75 <- filter(deaths_MA, age >= 75 & age < 85, sex == 1)
deaths_MA_m85 <- filter(deaths_MA, age >= 85, sex == 1)
deaths_MA_f65 <- filter(deaths_MA, age >= 65 & age < 75, sex == 2)
deaths_MA_f75 <- filter(deaths_MA, age >= 75 & age < 85, sex == 2)
deaths_MA_f85 <- filter(deaths_MA, age >= 85, sex == 2)

# Sum deaths within zipcode
deathssum <- deaths_MA %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_m65 <- deaths_MA_m65 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_m75 <- deaths_MA_m75 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_m85 <- deaths_MA_m85 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_f65 <- deaths_MA_f65 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_f75 <- deaths_MA_f75 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()
deathssum_f85 <- deaths_MA_f85 %>% group_by(zip) %>% 
  summarise(sum_deaths=sum(dead),
            .groups = 'drop') %>%
  as.data.frame()

#Aggregate deaths
deaths_MA_zip <- deaths_MA[!duplicated(deaths_MA$zip), ][,c(1,22:64)]
deaths_MA_sums <- data.frame(zip=deaths_MA[!duplicated(deaths_MA$zip), ]$zip, 
                             sum_deaths=deathssum$sum_deaths, 
                             people=as.numeric(table(deaths_MA$zip)))
deaths_MA_m65_sums <- data.frame(zip=deaths_MA_m65[!duplicated(deaths_MA_m65$zip), ]$zip, 
                                 sum_deathsm65=deathssum_m65$sum_deaths, 
                                 peoplem65=as.numeric(table(deaths_MA_m65$zip)))
deaths_MA_m75_sums <- data.frame(zip=deaths_MA_m75[!duplicated(deaths_MA_m75$zip), ]$zip, 
                                 sum_deathsm75=deathssum_m75$sum_deaths, 
                                 peoplem75=as.numeric(table(deaths_MA_m75$zip)))
deaths_MA_m85_sums <- data.frame(zip=deaths_MA_m85[!duplicated(deaths_MA_m85$zip), ]$zip, 
                                 sum_deathsm85=deathssum_m85$sum_deaths, 
                                 peoplem85=as.numeric(table(deaths_MA_m85$zip)))
deaths_MA_f65_sums <- data.frame(zip=deaths_MA_f65[!duplicated(deaths_MA_f65$zip), ]$zip, 
                                 sum_deathsf65=deathssum_f65$sum_deaths, 
                                 peoplef65=as.numeric(table(deaths_MA_f65$zip)))
deaths_MA_f75_sums <- data.frame(zip=deaths_MA_f75[!duplicated(deaths_MA_f75$zip), ]$zip, 
                                 sum_deathsf75=deathssum_f75$sum_deaths, 
                                 peoplef75=as.numeric(table(deaths_MA_f75$zip)))
deaths_MA_f85_sums <- data.frame(zip=deaths_MA_f85[!duplicated(deaths_MA_f85$zip), ]$zip, 
                                 sum_deathsf85=deathssum_f85$sum_deaths, 
                                 peoplef85=as.numeric(table(deaths_MA_f85$zip)))

deaths_MA_merged <- plyr::join_all(list(deaths_MA_zip,deaths_MA_sums,deaths_MA_m65_sums,deaths_MA_m75_sums,deaths_MA_m85_sums,
                                        deaths_MA_f65_sums,deaths_MA_f75_sums,deaths_MA_f85_sums), by='zip', type='left')

### Shapefile
USA <- st_read("/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/Zipcode_Info/polygon/ESRI16USZIP5_POLY_WGS84.shp")
names(USA)[names(USA) == "ZIP"] <- "zip"
USA$zip <- as.numeric(USA$zip)
MA <- filter(USA, STATE=="MA")

MA$zip_medicare <- MA$zip %in% deaths_MA_merged$zip
ggplot() + geom_sf(data = MA)+ 
  geom_sf(data = MA, aes(fill = zip_medicare))

### merge shapefile with deaths
deaths_shp_MA <- left_join(MA, deaths_MA_merged, by = "zip")

# saveRDS(deaths_shp_MA, file = "deaths_shp_MA.rds")

### Calculate rates
# deaths_shp_MA$rate_deaths <- deaths_shp_MA$sum_deathsf85/deaths_shp_MA$people
# ggplot() + geom_sf(data = deaths_shp_MA, aes(fill = rate))

# cvdzipNA <- deaths_MA_zip$zip[!deaths_MA_zip$zip %in% cvd_MA_zip$zip]
# deaths_shp_MA$cvd[deaths_shp_MA$zip == cvdzipNA] <- 0
# deaths_shp_MA <- filter(deaths_shp_MA, is.na(deaths) ==  FALSE, is.na(population) ==  FALSE)
# deaths_shp_MA$rate <-  deaths_shp_MA$cvd/deaths_shp_MA$people

# Investigate missing zips
# missing <- deaths_shp_MA$zip[unique(c(which(is.na(deaths_shp_MA$sum_deathsm65)),which(is.na(deaths_shp_MA$sum_deathsm75)),which(is.na(deaths_shp_MA$sum_deathsm85)),
#              which(is.na(deaths_shp_MA$sum_deathsf65)),which(is.na(deaths_shp_MA$sum_deathsf75)),which(is.na(deaths_shp_MA$sum_deathsf85))))]
# deaths_shp_MA$zip_missing <- deaths_shp_MA$zip %in% missing
# mapview(deaths_shp_MA, zcol = "zip_missing")

# deaths_shp[missing,]$population
# deaths_shp_MA$zip_missing <- !deaths_shp_MA$zip %in% deaths_MA_merged$zip
# ggplot() + geom_sf(data = deaths_shp_MA, aes(fill = zip_missing))

# Input 0 events if there are people in that zip
# which(is.na(deaths_shp_MA$sum_deathsm65)==TRUE & is.na(deaths_shp_MA$peoplem65)==FALSE)
# which(is.na(deaths_shp_MA$sum_deathsm75)==TRUE & is.na(deaths_shp_MA$peoplem75)==FALSE)
# which(is.na(deaths_shp_MA$sum_deathsm85)==TRUE & is.na(deaths_shp_MA$peoplem85)==FALSE)
# which(is.na(deaths_shp_MA$sum_deathsf65)==TRUE & is.na(deaths_shp_MA$peoplef65)==FALSE)
# which(is.na(deaths_shp_MA$sum_deathsf75)==TRUE & is.na(deaths_shp_MA$peoplef75)==FALSE)
# which(is.na(deaths_shp_MA$sum_deathsf85)==TRUE & is.na(deaths_shp_MA$peoplef85)==FALSE)

### Contiguity Matrix Weights
zip_remove <- deaths_shp_MA$zip[c(24,71,88)]
deaths_shp_MA[which(deaths_shp_MA$zip%in%zip_remove),]$PO_NAME
deaths_shp_MA <- filter(deaths_shp_MA, !zip %in% zip_remove)

data_agg <- cbind(deaths_shp_MA$sum_deathsm65, deaths_shp_MA$sum_deathsm75,deaths_shp_MA$sum_deathsm85,
                  deaths_shp_MA$sum_deathsf65,deaths_shp_MA$sum_deathsf75,deaths_shp_MA$sum_deathsf85,
                  deaths_shp_MA$peoplem65, deaths_shp_MA$peoplem75, deaths_shp_MA$peoplem85, 
                  deaths_shp_MA$peoplef65, deaths_shp_MA$peoplef75, deaths_shp_MA$peoplef85,
                  deaths_shp_MA$zip, deaths_shp_MA$pm25_ensemble, deaths_shp_MA$poverty,
                  deaths_shp_MA$popdensity,deaths_shp_MA$medianhousevalue,deaths_shp_MA$pct_blk,
                  deaths_shp_MA$medhouseholdincome,deaths_shp_MA$pct_owner_occ,deaths_shp_MA$hispanic,
                  deaths_shp_MA$education,deaths_shp_MA$population,deaths_shp_MA$smoke_rate,
                  deaths_shp_MA$mean_bmi,deaths_shp_MA$amb_visit_pct,deaths_shp_MA$a1c_exm_pct)
data.model_agg <- model.matrix(~data_agg-1)
X_agg <- scale(data.model_agg[,-(1:13)])
colnames(X_agg) <- c("pm25_ensemble","poverty","popdensity",
                     "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
                     "hispanic","education","population","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
y_deaths_m65 <-  data.model_agg[,1]
y_deaths_m75 <-  data.model_agg[,2]
y_deaths_m85 <-  data.model_agg[,3]
y_deaths_f65 <-  data.model_agg[,4]
y_deaths_f75 <-  data.model_agg[,5]
y_deaths_f85 <-  data.model_agg[,6]
E_m65 <-  data.model_agg[,7]
E_m75 <-  data.model_agg[,8]
E_m85 <-  data.model_agg[,9]
E_f65 <-  data.model_agg[,10]
E_f75 <-  data.model_agg[,11]
E_f85 <-  data.model_agg[,12]
zips.model_agg <- data.model_agg[,13]
deaths_shp_MA_agg <- filter(deaths_shp_MA, zip %in% zips.model_agg)
saveRDS(deaths_shp_MA_agg, file = "deaths_shp_MA_agg.rds")

data_overall <- cbind(deaths_shp_MA$sum_deaths, deaths_shp_MA$people,
                      deaths_shp_MA$zip, deaths_shp_MA$pm25_ensemble, deaths_shp_MA$poverty,
                      deaths_shp_MA$popdensity,deaths_shp_MA$medianhousevalue,deaths_shp_MA$pct_blk,
                      deaths_shp_MA$medhouseholdincome,deaths_shp_MA$pct_owner_occ,deaths_shp_MA$hispanic,
                      deaths_shp_MA$education,deaths_shp_MA$population,deaths_shp_MA$smoke_rate,
                      deaths_shp_MA$mean_bmi,deaths_shp_MA$amb_visit_pct,deaths_shp_MA$a1c_exm_pct)
data.model_overall <- model.matrix(~data_overall-1)
X_overall <- scale(data.model_overall[,-(1:3)])
colnames(X_overall) <- c("pm25_ensemble","poverty","popdensity",
                         "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
                         "hispanic","education","population","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
y_deaths <-  data.model_overall[,1]
E <-  data.model_overall[,2]
zips.model_overall <- data.model_overall[,3]
deaths_shp_MA_overall <- filter(deaths_shp_MA, zip %in% zips.model_overall)
saveRDS(deaths_shp_MA_overall, file = "deaths_shp_MA_overall.rds")