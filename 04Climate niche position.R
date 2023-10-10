# Calculate pca, kernel density, distances from samples to peak density, and angle of vector between sample points and peak
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ##########################################################################################
library(ade4)
# library(adehabitatMA)
library(here)
# library(raster)
library(geodata) # To replace raster getData
library(terra)
library(tidyverse)
source(here::here("Functions.R"))

GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Subgroups_02.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), header = TRUE, stringsAsFactors = FALSE)

# IUCN and GBIF names are not always matching
# TODO: This will be fixed in 01.1
GBIF_Data <- GBIF_Data %>%
  dplyr::mutate(species = case_when(species == "Mustela vison"    ~ "Neovison vison",
                                    species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                                    species == "Pekania pennanti" ~ "Martes pennanti",
                                    TRUE                          ~ species))

# Prepping climate data:
## load all climate data, subset variables I'm interested in:
### 5: Max Temperature of Warmest Month
### 6: Min Temperature of Coldest Month
### 12: Annual Precipitation
## extract climate data into dataframe for use in functions

# old:
## bio.dat_old <- raster::getData('worldclim', var = 'bio', res = 10)
## BIO_050612_old <- raster::subset(bio.dat_old, c(5, 6, 12))
## clim_old <- as.data.frame(na.omit(cbind(clim.xy, raster::extract(BIO_050612, clim.xy)))) 

BIO_050612 <- geodata::worldclim_global(var = "bio", res = 10, path = here::here("Data/WorldClim")) %>% 
  terra::subset(c(5, 6, 12))

Bioclim_DF <- terra::extract(x = BIO_050612,
                             y = terra::cells(BIO_050612),
                             xy = TRUE)

Clim_Variables <- c("wc2.1_10m_bio_5", "wc2.1_10m_bio_6", "wc2.1_10m_bio_12")

# PCA and climatic niche ##########################################################################################

PCA_Full <- ade4::dudi.pca(Bioclim_DF[, Clim_Variables], center = T, scale = T, scannf = F, nf = 2)

# columns to keep during analysis
gmpd.cols <- c("HostCorrectedName",
               "ParasiteCorrectedName", 
               "HostsSampled", 
               "Longitude", 
               "Latitude", 
               "Prevalence", 
               "EquatorwardsProp",
               "EquatorwardsDist",
               "MedianDist",
               "MedianProp",
               "AboveMedn",
               "subgroup",
               "CleanAll",
               "CleanSub",
               "RestrAll",
               "RestrSub",
               "RangeMethod",
               "RangeTaxonLvl")

GMPD_Kernel_Data <- kd_prep(clim.raw = BIO_050612,
                            spat.dat = GBIF_Data,
                            samp.dat = GMPD_Data[, gmpd.cols],
                            pca.full = PCA_Full,
                            bioclim.full.df = Bioclim_DF)

GMPD_Kernel_Data <- clim_density(climate.pca.scores = GMPD_Kernel_Data$pca.climate,
                                 species.pca.scores = GMPD_Kernel_Data$pca.species,
                                 samples.pca.scores = GMPD_Kernel_Data$pca.samples)

# Save output
saveRDS(KS_Output, here::here("Data/Data back ups/GMPD_Kernel_Data_04.rds"))
