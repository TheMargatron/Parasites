# Calculate pca, kernel density, distances from samples to peak density, and angle of vector between sample points and peak
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ##########################################################################################
library(ade4)
library(here)
library(geodata) # To replace raster getData
library(terra)
library(tidyverse)
source(here::here("Functions.R"))

GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data_02.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Distances_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), header = TRUE, stringsAsFactors = FALSE)

# Prepping climate data:
## load all climate data, (resolution = 10arcmin) subset variables I'm interested in:
### 5: Max Temperature of Warmest Month
### 6: Min Temperature of Coldest Month
### 12: Annual Precipitation
## extract climate data into dataframe for use in functions


BIO_050612 <- geodata::worldclim_global(var = "bio", res = 10, path = here::here("Data/WorldClim")) %>% 
  terra::subset(c(5, 6, 12))

Bioclim_DF <- terra::extract(x = BIO_050612,
                             y = terra::cells(BIO_050612),
                             xy = TRUE)

Clim_Variables <- c("wc2.1_10m_bio_5", "wc2.1_10m_bio_6", "wc2.1_10m_bio_12")

# PCA and climatic niche ##########################################################################################

PCA_Full <- ade4::dudi.pca(Bioclim_DF[, Clim_Variables], center = T, scale = T, scannf = F, nf = 2)

# columns to keep during analysis
# TODO: update this list
gmpd.cols <- c("HostCorrectedName",
               "ParasiteCorrectedName", 
               "ParType",
               "ParPhylum",
               "Group",
               "HostsSampled", 
               "SampleSize",
               "Longitude", 
               "Latitude", 
               "Prevalence", 
               "EquatorwardsProp",
               "EquatorwardsDist",
               "MedianDist",
               "MedianProp",
               "AboveMedn",
               # "subgroup",
               "CleanAll",
               # "CleanSub",
               "RestrAll",
               # "RestrSub",
               "RangeMethod",
               "RangeTaxonLvl")
# Urocyon littoralis doesn't have a big enough range to occupy two 20arcmin raster cells
# Also each island population is very separated so I should remove it anyway

GMPD_Distances_Data <- GMPD_Distances_Data %>% filter(HostCorrectedName != "Urocyon littoralis")

# Running at 10arcmin res
GMPD_Kernel_Data <- kd_prep(clim.raw = BIO_050612,
                            spat.dat = GBIF_Data,
                            samp.dat = GMPD_Distances_Data[, gmpd.cols], # doesn't matter that I'm not filtering because of rasterisation
                            pca.full = PCA_Full,
                            bioclim.full.df = Bioclim_DF,
                            density.resolution = 10/60) # to match resolution of worldclim data

GMPD_Kernel_Data <- clim_density(climate.pca.scores = GMPD_Kernel_Data$pca.climate,
                                 species.pca.scores = GMPD_Kernel_Data$pca.species,
                                 samples.pca.scores = GMPD_Kernel_Data$pca.samples)

GMPD_Kernel_Data_half <- kd_prep(clim.raw = BIO_050612,
                                 spat.dat = GBIF_Data,
                                 samp.dat = GMPD_Distances_Data[, gmpd.cols],
                                 pca.full = PCA_Full,
                                 bioclim.full.df = Bioclim_DF,
                                 density.resolution = 5/60) # half resolution of worldclim data

GMPD_Kernel_Data_half <- clim_density(climate.pca.scores = GMPD_Kernel_Data_half$pca.climate,
                                      species.pca.scores = GMPD_Kernel_Data_half$pca.species,
                                      samples.pca.scores = GMPD_Kernel_Data_half$pca.samples)

GMPD_Kernel_Data_double <- kd_prep(clim.raw = BIO_050612,
                                   spat.dat = GBIF_Data,
                                   samp.dat = GMPD_Distances_Data[, gmpd.cols],
                                   pca.full = PCA_Full,
                                   bioclim.full.df = Bioclim_DF,
                                   density.resolution = 20/60) # double resolution of worldclim data

GMPD_Kernel_Data_double <- clim_density(climate.pca.scores = GMPD_Kernel_Data_double$pca.climate,
                                        species.pca.scores = GMPD_Kernel_Data_double$pca.species,
                                        samples.pca.scores = GMPD_Kernel_Data_double$pca.samples)


# Save output
saveRDS(list(GMPD_Kernel_Data, GMPD_Kernel_Data_half, GMPD_Kernel_Data_double), here::here("Data/Data back ups/GMPD_Kernel_Data_04.rds"))

GMPD_Climate_Data <- GMPD_Kernel_Data$samples.out %>% 
  dplyr::bind_rows() 

GMPD_Climate_Data <- GMPD_Kernel_Data_half$samples.out %>% 
  dplyr::bind_rows() %>% 
  full_join(GMPD_Climate_Data, 
            by = names(.)[!grepl("distance|angle|density", names(.), ignore.case = TRUE)],
            suffix = c("_5", ""))

GMPD_Climate_Data <- GMPD_Kernel_Data_double$samples.out %>% 
  dplyr::bind_rows() %>% 
  full_join(GMPD_Climate_Data, 
            by = names(.)[!grepl("distance|angle|density", names(.), ignore.case = TRUE)],
            suffix = c("_20", ""))

write.csv(GMPD_Climate_Data, here::here("Data/Data back ups/GMPD_Climate_Data_half_04.csv"), row.names = FALSE)

