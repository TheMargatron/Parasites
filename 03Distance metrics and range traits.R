# Calculating distance metrics (two methods) and extracting range traits
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ##########################################################################################

library(here)               #
library(tidyverse)          # beware of conflicts (mainly with raster)

source(here::here("Functions.R"))

GMPD_Data_res_sub <- read.csv(here::here("Data/Data back ups/GMPD_Data_res_sub_01.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data_res_all <- read.csv(here::here("Data/Data back ups/GMPD_Data_res_all_01.csv"), header = TRUE, stringsAsFactors = FALSE)

GMPD_Data_cln_sub <- read.csv(here::here("Data/Data back ups/GMPD_Data_cln_sub_01.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data_cln_all <- read.csv(here::here("Data/Data back ups/GMPD_Data_cln_all_01.csv"), header = TRUE, stringsAsFactors = FALSE)

#GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data.csv"), header = TRUE, stringsAsFactors = FALSE)

# Prep ####
# Creating Trait dataframe for later and simplifying GMPD_Data to essentials

## can no longer remember why I created this, was before I started using git..
# GMPD_Data_cln_all$method <- "clean"
# 
# GMPD_Trait_Data <- GMPD_Data_res_all %>%
#   mutate(method = "restricted") %>%
#   bind_rows(GMPD_Data_cln_all) %>%
#   dplyr::select(HostCorrectedName, ParasiteCorrectedName, Group, HostOrder, HostFamily, HostEnvironment, ParType, ParPhylum, ParClass) %>%
#   unique()

GMPD_Data_res_all <- GMPD_Data_res_all %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, 
                Citation, LocationName, Longitude, Latitude, 
                PopulationType, SamplingBasis, Prevalence, 
                HostsSampled, HostSex, HostAge, SamplingType, 
                subgroup)

GMPD_Data_res_sub <- GMPD_Data_res_sub %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, 
                Citation, LocationName, Longitude, Latitude, 
                PopulationType, SamplingBasis, Prevalence, 
                HostsSampled, HostSex, HostAge, SamplingType, subgroup)

GMPD_Data_cln_all <- GMPD_Data_cln_all %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, 
                Citation, LocationName, Longitude, Latitude, 
                PopulationType, SamplingBasis, Prevalence, 
                HostsSampled, HostSex, HostAge, SamplingType, 
                subgroup)

GMPD_Data_cln_sub <- GMPD_Data_cln_sub %>%
  dplyr::select(HostCorrectedName, ParasiteCorrectedName, 
                Citation, LocationName, Longitude, Latitude, 
                PopulationType, SamplingBasis, Prevalence, 
                HostsSampled, HostSex, HostAge, SamplingType, subgroup)

# Distance metrics and range traits ###########################################################################

Distances_Data_res_all <- range_distances(dat = GMPD_Data_res_all, range.pol = IUCN_Native_Data, method = "iucn")
Distances_Data_res_all$DistanceMetrics$method <- "restricted"
#Distances_Data_res_sub <- range_distances(dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn") # for when subgroup method is ready

Distances_Data_cln_all <- range_distances(dat = GMPD_Data_cln_all, range.dat = GBIF_Data, method = "gbif")
Distances_Data_cln_all$DistanceMetrics$method <- "cleaned"
#Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif") # for when subgroup method is ready

GMPD_Distances_Data <- bind_rows(Distances_Data_res_all$DistanceMetrics, Distances_Data_cln_all$DistanceMetrics)
Range_Traits <- bind_rows(Distances_Data_res_all$RangeTraits, Distances_Data_cln_all$RangeTraits)

write.csv(GMPD_Distances_Data, file = here::here("Data/Data back ups/GMPD_Distances_Data.csv"), row.names = FALSE)
write.csv(Range_Traits, file = here::here("Data/Data back ups/Range_Traits.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
