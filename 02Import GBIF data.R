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

Taxon_Key <- taxize::get_gbifid_(Hostlist, method = "backbone")
Taxon_Key <- lapply(names(Taxon_Key), function(name) {
  Taxon_Key[[name]]["HostCorrectedName"] <- name
  return(Taxon_Key[[name]])
})

Taxon_Key <- Taxon_Key %>%
  bind_rows() %>%
  filter(class == "Mammalia")

# Probably a better way of doing this but this works for now
# Using the taxonkey which gives the same number of search results as the species name 
Taxon_Key$CountKey <- unlist(lapply(Taxon_Key$usagekey, function(key) occ_count(taxonKey = key, georeferenced = TRUE)))
Taxon_Key$CountSpc <- unlist(lapply(Taxon_Key$HostCorrectedName, function(spc) occ_search(scientificName = spc, hasCoordinate = TRUE, limit = 0)$meta$count))

Taxon_Key <- Taxon_Key %>%
  filter(CountKey == CountSpc) %>%
  select(usagekey, HostCorrectedName)

################## actual download ##################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

occ_download(
  pred_in("taxonKey", Taxon_Key$usagekey),
  pred("classKey", ""),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

occ_download_get("0317209-200613084148143", path = here::here("Data/GBIF/"))

Occ2_Temp <- occ_download_import(key = "0317209-200613084148143", path = here::here("Data/GBIF/"))


########### old code ##############

gbif.import <- function(host){
  gbifsize <- occ_search(scientificName = host, hasCoordinate = TRUE, return = "meta")$count
  if(gbifsize < 100000){
    temp.gbif <- occ_data(scientificName = host, hasCoordinate = T, limit = 100000)
    #temp.gbif <- select(temp.gbif, one_of(gbif.cols))
    
  } else{
    warning(paste("Too many occurrences, manually download and import", host, "data"))
    #paste("Too many occurrences, manually download and import", host, "data")
    #NULL
  }
}


