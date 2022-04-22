####################################################################################################
##### CALCULATE POSITION OF PARASITE SAMPLE LOCATIONS IN GEOGRAPHIC SPACE AND CALCULATE RANGE TRAITS ######################
##### Written by: Margaret Bolton mb804(at)exeter.ac.uk and Regan Early r.early@exeter.ac.uk #######
##### Written on: ... ##############################################################################
##### Modified on: 22ndApril 2022   ###################################################################
####################################################################################################

##### Outputs #####
# GMPD_Distances_Data_all.csv and GMPD_Distances_Data_sub.csv: the location of each parasite sample point relative to the most equator-ward latitude in the mammal host's range, or relative to the median latitude in the mammal host's range, and each of these expressed in metres and as a proportion of the host mammals' latitudinal range span. 'All' refers to when a mammal's range is assigned to thea single species, and 'sub' when the ranges of a mammal's sub species are used. 
# Range_Traits_all.csv and Range_Traits_sub.csv: The max, min, span, area (when calculated using IUCN polygons), and median of the host mammal's range

##### Libraries and data #####
# library(here)               # I couldn't get this to work
library(tidyverse)          # beware of conflicts (mainly with raster)

wd.dat <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/DATA/22ndApril2022"
wd.code <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/CODE_GITHUB_REGAN/RANGE_POSITION"
source(paste0(wd.code, "/range_functions.R")) 

GMPD_Data <- read.csv(paste0(wd.dat, "/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite data associated with mammal sub-species. "Restricted"

GMPD_Data_res_sub <- GMPD_Data[GMPD_Data$RestrSub==T,] ## Parasite data associated with mammal species. "res_sub" means parasite samples and mammal distribution locations are restricted to within the each mammal species' subgroup IUCN polygons
GMPD_Data_res_all <- GMPD_Data[GMPD_Data$RestrAll==T,] ## Parasite data associated with mammal species. "res_all" means parasite samples and mammal distribution locations are restricted to within the whole species' IUCN polygons

GMPD_Data_cln_sub <- GMPD_Data[GMPD_Data$CleanSub==T,] ## Parasite data associated with mammal sub-species. "cln_sub" means data are cleaned, but includes parasite samples and mammal distribution locations outside of IUCN polygons, but records are assigned to the mammal species subgroup
GMPD_Data_cln_all <- GMPD_Data[GMPD_Data$CleanAll==T,] ## Parasite data associated with mammal species. "cln_all" means data are cleaned, but includes parasite samples and mammal distribution locations outside of IUCN polygons, and records are assigned to the mammal species as a whole

IUCN_Native_Data <- readRDS(file=paste0(wd.dat, "/IUCN_Native_Data_01.R")) # range polygons. I checked with Margaret on 6th April 2022 and this is still the most up to date version of the IUCN range polygons. 

GBIF_Data <- read.csv(paste0(wd.dat, "/GBIF_Subgroups_02.csv"))

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
Distances_Data_res_all <- range_distances(para.dat = GMPD_Data_res_all, range.pol = IUCN_Native_Data, method = "iucn", tax="species") ## Should be run with GBIF and IUCN
Distances_Data_res_all$DistanceMetrics$method <- Distances_Data_res_all$RangeTraits$method <- "restricted"
#Distances_Data_res_sub <- range_distances(dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn") # for when subgroup method is ready

Distances_Data_cln_all <- range_distances(para.dat = GMPD_Data_cln_all, range.dat = GBIF_Data, method = "gbif", tax="species") ## Can't use IUCN data because parasite data has not been restricted to the IUCN polygons.
Distances_Data_cln_all$DistanceMetrics$method <- Distances_Data_cln_all$RangeTraits$method <- "cleaned"
#Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif") # for when subgroup method is ready

GMPD_Distances_Data <- bind_rows(Distances_Data_res_all$DistanceMetrics, Distances_Data_cln_all$DistanceMetrics)
Range_Traits <- bind_rows(Distances_Data_res_all$RangeTraits, Distances_Data_cln_all$RangeTraits)

write.csv(GMPD_Distances_Data, file = paste0(wd.dat, "/GMPD_Distances_Data_all.csv"), row.names = FALSE)
write.csv(Range_Traits, file = paste0(wd.dat, "/Range_Traits_all.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)

##### Distance metrics and range traits for subgroups #####
Distances_Data_res_sub <- range_distances(para.dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn", tax="sub")
Distances_Data_res_sub$DistanceMetrics$method <- Distances_Data_res_sub$RangeTraits$method <- "restricted"
#Distances_Data_res_sub <- range_distances(dat = GMPD_Data_res_sub, range.pol = IUCN_Native_Data, method = "iucn") # for when subgroup method is ready

Distances_Data_cln_sub <- range_distances(para.dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif", tax="sub")
Distances_Data_cln_sub$DistanceMetrics$method <- Distances_Data_cln_sub$RangeTraits$method <- "cleaned"
#Distances_Data_cln_sub <- range_distances(dat = GMPD_Data_cln_sub, range.dat = GBIF_Data, method = "gbif") # for when subgroup method is ready

GMPD_Distances_Data <- bind_rows(Distances_Data_res_all$DistanceMetrics, Distances_Data_cln_all$DistanceMetrics)
Range_Traits <- bind_rows(Distances_Data_res_all$RangeTraits, Distances_Data_cln_all$RangeTraits)

write.csv(GMPD_Distances_Data, file = paste0(wd.dat, "/GMPD_Distances_Data_sub.csv"), row.names = FALSE)
write.csv(Range_Traits, file = paste0(wd.dat, "/Range_Traits_sub.csv"), row.names = FALSE)
#write.csv(GMPD_Trait_Data, file = here::here("Data/Data back ups/GMPD_Trait_Data.csv"), row.names = FALSE)
