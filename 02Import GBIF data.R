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

Taxon_Key <- taxize::get_gbifid_(Taxonlist, method = "backbone")
Taxon_Key <- lapply(names(Taxon_Key), function(name) {
  Taxon_Key[[name]]["HostCorrectedName"] <- name
  return(Taxon_Key[[name]])
})

Taxon_Key <- Taxon_Key %>%
  bind_rows() %>%
  filter(class == "Mammalia")

############################################## Not using for now ##############################################
# Probably a better way of doing this but this works for now
# Using the taxonkey which gives the same number of search results as the species name 
# Taxon_Key$CountKey <- unlist(lapply(Taxon_Key$usagekey, function(key) occ_count(taxonKey = key, georeferenced = TRUE)))
# Taxon_Key$CountSpc <- unlist(lapply(Taxon_Key$HostCorrectedName, function(spc) occ_search(scientificName = spc, hasCoordinate = TRUE, limit = 0)$meta$count))

# Taxon_Key <- Taxon_Key %>%
#   filter(CountKey == CountSpc) %>%
#   select(usagekey, HostCorrectedName)

############################################## actual download ################################################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

Download_Key <- occ_download(
  pred_in("taxonKey", Taxon_Key$usagekey),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

Download_Get <- occ_download_get(Download_Key, path = here::here("Data/GBIF/"), overwrite = TRUE)

GBIF_Data <- occ_download_import(Download_Get, path = here::here("Data/GBIF/"))

