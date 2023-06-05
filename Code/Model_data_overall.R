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
deaths_shp_MA_overall <- readRDS("deaths_shp_MA_overall.rds")
data_overall <- cbind(deaths_shp_MA_overall$sum_deaths, deaths_shp_MA_overall$people,
                      deaths_shp_MA_overall$zip, deaths_shp_MA_overall$pm25_ensemble, deaths_shp_MA_overall$poverty,
                      deaths_shp_MA_overall$popdensity,deaths_shp_MA_overall$medianhousevalue,deaths_shp_MA_overall$pct_blk,
                      deaths_shp_MA_overall$medhouseholdincome,deaths_shp_MA_overall$pct_owner_occ,deaths_shp_MA_overall$hispanic,
                      deaths_shp_MA_overall$education,deaths_shp_MA_overall$smoke_rate,
                      deaths_shp_MA_overall$mean_bmi,deaths_shp_MA_overall$amb_visit_pct,deaths_shp_MA_overall$a1c_exm_pct)
data.model_overall <- model.matrix(~data_overall-1)
X_overall <- scale(data.model_overall[,-(1:3)])
colnames(X_overall) <- c("pm25_ensemble","poverty","popdensity",
                         "medianhousevalue","pct_blk","medhouseholdincome","pct_owner_occ",
                         "hispanic","education","smoke_rate","mean_bmi","amb_visit_pct","a1c_exm_pct")
y_deaths <-  data.model_overall[,1]
E <-  data.model_overall[,2]
zips.model_overall <- data.model_overall[,3]
deaths_shp_MA_overall_overall <- filter(deaths_shp_MA_overall, zip %in% zips.model_overall)

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
MA.contig_overall <- st_contiguity(deaths_shp_MA_overall)
MA.contig_sf <- st_as_edges(st_centroid(deaths_shp_MA_overall$geometry), nb = MA.contig_overall)
weights.contig.B <- st_weights(MA.contig_overall, style = "B")
W <- wt_as_matrix(MA.contig_overall, weights.contig.B)

N = nrow(data.model_overall)
neigh = adjlist(W, N)
numneigh = apply(W,2,sum)
nbs = mungeCARdata4stan(neigh, numneigh)
N = nbs$N ; node1 = nbs$node1 ; node2 = nbs$node2 ; N_edges = nbs$N_edges
n = nrow(X_overall); p = ncol(X_overall);y_deaths = data.model_overall[,1]; y_cvd = data.model_overall[,2] = W_n = sum(W) / 2;W = W; E=data.model_overall[,3]

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