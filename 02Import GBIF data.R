# GBIF import and data clean
# Written by Margaret Bolton mb804(at)exeter.ac.uk 
# Data citation: "GBIF.org (02 July 2021) GBIF Occurrence Download  https://doi.org/10.15468/dl.9hewb3"

############################################## Libraries and data #############################################

library(here)
library(rgbif)
library(taxize)
library(tidyverse)

GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data.csv"), header = TRUE, stringsAsFactors = FALSE)
Hostlist <- unique(GMPD_Data$HostCorrectedName)

############################################## Getting taxon keys #############################################
# Adding synonymous species names not picked up by taxize
Taxonlist <- c(Hostlist, "Pekania pennanti", "Taurotragus oryx", "Mustela vison")

Taxon_Keys <- taxize::get_gbifid_(Taxonlist, method = "backbone")
Taxon_Keys <- lapply(Taxonlist, function(name) {
  Taxon_Keys[[name]]["HostCorrectedName"] <- name
  return(Taxon_Keys[[name]])
})

Taxon_Keys <- Taxon_Keys %>%
  bind_rows() %>%
  filter(class == "Mammalia")

Taxon_Key_Search <- Taxon_Keys %>%                     # keeping full list separate to match misnamed species later
  filter(status == "ACCEPTED" & matchtype == "EXACT") 

############################################## actual download ################################################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

Download_Key <- occ_download(
  pred_in("taxonKey", Taxon_Key_Search$usagekey),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

Download_Get <- occ_download_get(Download_Key, path = here::here("Data/GBIF/"), overwrite = TRUE)

GBIF_Data <- occ_download_import(Download_Get, path = here::here("Data/GBIF/"))

saveRDS(Download_Key, here::here("Data/GBIF/Download_Key"))
saveRDS(Download_Get, here::here("Data/GBIF/Download_Get"))
