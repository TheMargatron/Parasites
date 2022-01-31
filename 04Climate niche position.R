# Climate niche functions and analysis
# Written by Margaret Bolton mb804(at)exeter.ac.uk (some isn't mine, check what)

# Libraries and data ##########################################################################################
library(ade4)
library(adehabitatMA)
library(here)
library(raster)
library(tidyverse)
source(here::here("Climate niche functions.R"))

GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data.csv"), header = TRUE, stringsAsFactors = FALSE)

# IUCN and GBIF names are not always matching
GBIF_Data <- GBIF_Data %>%
  mutate(species = case_when(GBIFName == "Neovison vison"  ~ "Mustela vison",
                            GBIFName == "Tragelaphus oryx" ~ "Taurotragus oryx",
                            GBIFName == "Martes pennanti"  ~ "Pekania pennanti",
                            TRUE                           ~ species))

bio.dat <- getData('worldclim', var = 'bio', res = 10)

BIO_050612 <- subset(bio.dat, c(5, 6, 12))

# PCA and climatic niche ##########################################################################################

gmpd.cols <- c("ParasiteCorrectedName", 
               "HostsSampled", 
               "Longitude", 
               "Latitude", 
               "Prevalence", 
               "EquatorwardsProp",
               "EquatorwardsDist",
               "MedianDist",
               "MedianProp")

# extract all bio data
clim.xy <- xyFromCell(BIO_050612, 1:ncell(BIO_050612))
clim <- as.data.frame(na.omit(cbind(clim.xy, raster::extract(BIO_050612, clim.xy)))) 
xVar <- c(3:5)
rm(clim.xy)

pca.cal <- dudi.pca(clim[xVar], center = T, scale = T, scannf = F, nf = 2)

# species.name <- "Genetta genetta" #example for test runs
# spat.dat <- GBIF_Data
# clim.dat <- BIO_050612
# par.dat <- GMPD_Data

ks.test <- function(species.name, spat.dat, clim.dat, par.dat){
  
  spat.dat <- spat.dat %>%
    filter(species == species.name) %>%
    dplyr::select(c("decimalLongitude", "decimalLatitude"))
  
  # (partially) account for sampling bias by gridding species occurence data and extracting filled cells
  rastr10 <- raster(resolution = 10/60)
  rast.sp <- rasterize(spat.dat, rastr10, fun = 'count')
  rast.filled <- Which(rast.sp, cells = TRUE)
  
  # extract bioclim variables for points where species occurs
  occ.xy <- xyFromCell(rast.sp, rast.filled)
  occ.xy <- as.data.frame(na.omit(cbind(occ.xy, raster::extract(x = clim.dat, y = occ.xy))))
  
  nvar <- length(xVar)     
  R = 100     
  
  #rbind occ data to glob data and add binary weighting column/vector
  #row.w.env <- c(rep(1, 1:nrow(clim)), rep(0, (nrow(clim) + 1):(nrow(clim) + nrow(occ.sp))))   #1 for background clim, zero for occ clim
  #Is this necessary? 
  
  scores.clim <- suprow(pca.cal, clim[, xVar])$lisup     #The pca scores for all climate
  scores.occ <- suprow(pca.cal, occ.xy[, xVar])$lisup    #The pca scores for current species
  
  gmpd.sp <- par.dat %>%
    filter(HostCorrectedName == species.name) %>%
    dplyr::select(gmpd.cols)
  
  gmpd.xy <- as.data.frame(na.omit(cbind(gmpd.sp, raster::extract(clim.dat, gmpd.sp[, c("Longitude", "Latitude")]))))
  scores.gmpd <- cbind(suprow(pca.cal, gmpd.xy[, c("bio5", "bio6", "bio12")])$lisup, gmpd.xy)      # pca scores for current species from gmpd
  
  z2 <- grid.clim(scores.clim, scores.occ, scores.gmpd, R)    #Niche_dyn_funcs_myversion
  #z2$sp.scores <- scores.gmpd
  return(z2)
}

Hostlist <- unique(GMPD_Data$HostCorrectedName)

ks.out <- lapply(Hostlist, ks.test, spat.dat = GBIF_Data, clim.dat = BIO_050612, par.dat = GMPD_Data)
names(ks.out) <- Hostlist

saveRDS(ks.out, file = "KS_data")
saveRDS(ks.out2, file = "KS_data2")

ks.out <- ks.out2
