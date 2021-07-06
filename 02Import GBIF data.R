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
  Taxon_Key[[name]]["SpeciesID"] <- name
  return(Taxon_Key[[name]])
})

Taxon_Key2_Temp <- Taxon_Key %>%
  bind_rows() %>%
  filter(class == "Mammalia")

# Probably a better way of doing this but this works for now
# Using the taxonkey which gives the same number of search results as the species name 
Taxon_Key2_Temp$CountKey <- unlist(lapply(Taxon_Key2_Temp$usagekey, function(key) occ_count(taxonKey = key, georeferenced = TRUE)))
Taxon_Key2_Temp$CountSpc <- unlist(lapply(Taxon_Key2_Temp$SpeciesID, function(spc) occ_search(scientificName = spc, hasCoordinate = TRUE, limit = 0)$meta$count))

Taxon_Key2_Temp <- Taxon_Key2_Temp %>%
  filter(CountKey == CountSpc) %>%
  select(usagekey, SpeciesID)

Occ_Temp <- unlist(lapply(Hostlist, function(host) occ_search(scientificName = host, hasCoordinate = TRUE, return = "meta")$count))
Occ_Temp <- data.frame(TaxonKey = Occ_Temp, SpeciesName = Hostlist)

Both_Temp <- merge(Occ_Temp, Taxon_Key2_Temp, by = "SpeciesID")

View(Both_Temp[which(Both_Temp$SpCount != Both_Temp$Count),])


################## actual download ##################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

occ_download(
  pred_in("taxonKey", Taxon_Key2_Temp[1, "usagekey"]),
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


