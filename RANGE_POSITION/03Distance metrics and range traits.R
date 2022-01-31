####################################################################################################
##### CALCULATE POSITION OF PARASITE SAMPLE LOCATIONS IN GEOGRAPHIC SPACE AND CALCULATE RANGE TRAITS ######################
##### Written by: Margaret Bolton mb804(at)exeter.ac.uk and Regan Early r.early@exeter.ac.uk #######
##### Written on: ... ##############################################################################
##### Modified on: 31 Jan 2022   ###################################################################
####################################################################################################

##### Outputs #####
# ?: ...

##### Libraries and data #####
library(here)               # I couldn't get this to work
library(tidyverse)          # beware of conflicts (mainly with raster)

wd.dat <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/DATA/31Jan2021"
wd.code <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/CODE_GITHUB_REGAN/RANGE_POSITION"
source(paste0(wd.code, "/range_functions.R")) # I couldn't get here to work

GMPD_Data_res_sub <- read.csv(paste0(wd.dat, "/GMPD_Data_res_sub_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite data associated with mammal sub-species. "Restricted"
GMPD_Data_res_all <- read.csv(paste0(wd.dat, "/GMPD_Data_res_all_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite data associated with mammal species. "Restricted"

GMPD_Data_cln_sub <- read.csv(paste0(wd.dat, "/GMPD_Data_cln_sub_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite data associated with mammal sub-species. "Cleaned"
GMPD_Data_cln_all <- read.csv(paste0(wd.dat, "/GMPD_Data_cln_all_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite data associated with mammal species. "Cleaned"

IUCN_Native_Data <- readRDS(file=paste0(wd.dat, "/IUCN_Native_Data_01.R")) # range polygons

GBIF_Data <- read.csv(paste0(wd.dat, "/GBIF_temp_subset.csv"))

##### Simplify data to essentials #####
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
                HostsSampled, HostSex, HostAge, SamplingType,
                subgroup)

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
                HostsSampled, HostSex, HostAge, SamplingType,
                subgroup)

GBIF_Data <- GBIF_Data %>%
  dplyr::select(species, taxonRank, scientificName, subgroup, decimalLongitude, decimalLatitude) %>%
  drop_na(decimalLongitude) %>%
  rename(lon = decimalLongitude, lat = decimalLatitude)

unique(GBIF_Data[,c("species", "scientificName", "subgroup", "taxonRank")])

# ### Temporarily altering the epithets to match those used in the GMPD file - just to make code work 
# GBIF_Data[GBIF_Data$infraspecificEpithet=="sonoriensis" & is.na(GBIF_Data$infraspecificEpithet)==F, "infraspecificEpithet"] <- "americana sonoriensis"

##### Distance metrics and range traits for species #####

Distances_Data_res_all <- range_distances(para.dat = GMPD_Data_res_all, range.pol = IUCN_Native_Data, method = "iucn", tax="species")
Distances_Data_res_all$DistanceMetrics$method <- "restricted"
#Distances_Data_res_sub <- range_distances(dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn") # for when subgroup method is ready

Distances_Data_cln_all <- range_distances(para.dat = GMPD_Data_cln_all, range.dat = GBIF_Data, method = "gbif", tax="species")
Distances_Data_cln_all$DistanceMetrics$method <- "cleaned"
#Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif") # for when subgroup method is ready

GMPD_Distances_Data <- bind_rows(Distances_Data_res_all$DistanceMetrics, Distances_Data_cln_all$DistanceMetrics)
Range_Traits <- bind_rows(Distances_Data_res_all$RangeTraits, Distances_Data_cln_all$RangeTraits)

write.csv(GMPD_Distances_Data, file = paste0(wd.dat, "/GMPD_Distances_Data.csv"), row.names = FALSE)
write.csv(Range_Traits, file = paste0(wd.dat, "/Range_Traits.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)

##### Distance metrics and range traits for subgroups #####

Distances_Data_res_all <- range_distances(para.dat = GMPD_Data_res_all, range.pol = IUCN_Native_Data, method = "iucn", tax="sub")
Distances_Data_res_all$DistanceMetrics$method <- "restricted"
#Distances_Data_res_sub <- range_distances(dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn") # for when subgroup method is ready

Distances_Data_cln_all <- range_distances(para.dat = GMPD_Data_cln_all, range.dat = GBIF_Data, method = "gbif", tax="sub")
Distances_Data_cln_all$DistanceMetrics$method <- "cleaned"
#Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif") # for when subgroup method is ready

GMPD_Distances_Data <- bind_rows(Distances_Data_res_all$DistanceMetrics, Distances_Data_cln_all$DistanceMetrics)
Range_Traits <- bind_rows(Distances_Data_res_all$RangeTraits, Distances_Data_cln_all$RangeTraits)

write.csv(GMPD_Distances_Data, file = paste0(wd.dat, "/GMPD_Distances_Data_sub.csv"), row.names = FALSE)
write.csv(Range_Traits, file = paste0(wd.dat, "/Range_Traits_sub.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
