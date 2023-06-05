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
library(tidyverse)
library(INLA)

### Data
deaths_shp_MA_agg <- readRDS("deaths_shp_MA_agg.rds")
data_agg <- cbind(deaths_shp_MA_agg$sum_deathsm65, deaths_shp_MA_agg$sum_deathsm75,deaths_shp_MA_agg$sum_deathsm85,
                  deaths_shp_MA_agg$sum_deathsf65,deaths_shp_MA_agg$sum_deathsf75,deaths_shp_MA_agg$sum_deathsf85,
                  deaths_shp_MA_agg$peoplem65, deaths_shp_MA_agg$peoplem75, deaths_shp_MA_agg$peoplem85,
                  deaths_shp_MA_agg$peoplef65, deaths_shp_MA_agg$peoplef75, deaths_shp_MA_agg$peoplef85,
                  deaths_shp_MA_agg$zip, deaths_shp_MA_agg$pm25_ensemble, deaths_shp_MA_agg$poverty,
                  deaths_shp_MA_agg$popdensity,deaths_shp_MA_agg$medianhousevalue,deaths_shp_MA_agg$pct_blk,
                  deaths_shp_MA_agg$medhouseholdincome,deaths_shp_MA_agg$pct_owner_occ,deaths_shp_MA_agg$hispanic,
                  deaths_shp_MA_agg$education,deaths_shp_MA_agg$smoke_rate,
                  deaths_shp_MA_agg$mean_bmi,deaths_shp_MA_agg$amb_visit_pct,deaths_shp_MA_agg$a1c_exm_pct)
data.model_agg <- model.matrix(~data_agg-1)
X_agg <- scale(data.model_agg[,-(1:13)])
colnames(X_agg) <- c("pm25_ensemble","poverty","popdensity",
                     "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
                     "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
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

### Functions
`%||%` <- function (x, y) {
  if (rlang::is_null(x))
    y
  else x
}
recreate_listw <- function(nb, wt) {
  which_style <- c(attr(wt, "W") %||% NA,
                   attr(wt, "B") %||% NA,
                   attr(wt, "C") %||% NA,
                   attr(wt, "U") %||% NA,
                   attr(wt, "minmax") %||% NA,
                   attr(wt, "S") %||% NA)
  
  possible_styles <- c("W", "B", "C", "U", "minmax", "S")
  
  if (!inherits(nb, "nb")) nb <- class_modify(nb, "nb")
  
  listw <- list(style = possible_styles[!is.na(which_style)],
                neighbours = nb,
                weights = wt)
  
  class(listw) <- c("listw", "nb", "list")
  
  listw
}
wt_as_matrix <- function(nb, wt) {
  listw <- recreate_listw(nb, wt)
  spdep::listw2mat(listw)
}
adjlist = function(W,N){ 
  adj=0
  for(i in 1:N){  
    for(j in 1:N){  
      if(W[i,j]==1){adj = append(adj,j)}
    }
  }
  adj = adj[-1]
  return(adj)}
mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}

### Contiguity Matrix Weights
MA.contig_agg <- st_contiguity(deaths_shp_MA_agg)
MA.contig_sf <- st_as_edges(st_centroid(deaths_shp_MA_agg$geometry), nb = MA.contig_agg)
weights.contig.B <- st_weights(MA.contig_agg, style = "B")
W <- wt_as_matrix(MA.contig_agg, weights.contig.B)

N = nrow(data.model_agg)
neigh = adjlist(W, N)
numneigh = apply(W,2,sum)
nbs = mungeCARdata4stan(neigh, numneigh)
N = nbs$N ; node1 = nbs$node1 ; node2 = nbs$node2 ; N_edges = nbs$N_edges
n = nrow(X_agg); p = ncol(X_agg);y_deaths = data.model_agg[,1]; y_cvd = data.model_agg[,2] = W_n = sum(W) / 2;W = W; E=data.model_agg[,3]

#Build the adjacency matrix using INLA library functions
adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
#The ICAR precision matrix (note! This is singular)
Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix
#Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)

# Compute the diagonal elements of the covariance matrix subject to the 
# constraint that the entries of the ICAR sum to zero.
#See the inla.qinv function help for further details.
Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))

#Compute the geometric mean of the variances, which are on the diagonal of Q.inv
scaling_factor = exp(mean(log(diag(Q_inv))))