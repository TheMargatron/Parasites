# GBIF import and data clean
# Written by Margaret Bolton mb804(at)exeter.ac.uk 
# Data citation: "GBIF.org (6 October 2021) GBIF Occurrence Download  https://www.gbif.org/occurrence/download/0018590-210914110416597"

############################################## Libraries and data #############################################

library(CoordinateCleaner)
library(countrycode)
# library(ggplot2)
library(here)
#library(lubridate)
#library(maps)
library(rgbif)
# library(rgdal)              # read shapefiles
# library(rgeos)
library(rnaturalearthdata)
library(sf)
library(taxize)
library(tidyverse)
library(tmap)
library(beepr)
# source(here::here("Functions.R"))

sf::sf_use_s2(FALSE) # For "invalid spherical geometry" errors
tmap::tmap_mode("view")

# GMPD_Raw_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), header = TRUE, stringsAsFactors = FALSE) 
GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE)
Native_DF <- read.csv(here::here("Data/Data back ups/Native_DF_01.csv"), header = TRUE, stringsAsFactors = FALSE)

# Legend_Text <- sort(unique(read.csv(here::here("Data/Data back ups/Native_DF_01.csv"), header = TRUE)$Status))

# IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))
IUCN_Data <- readRDS(here::here("Data/Data back ups/IUCN_Data_01"))

data("World")
Projection_String <- sf::st_crs(IUCN_Data)

Hostlist <- sort(unique(GMPD_Data$HostCorrectedName))

Already_Prepped <- TRUE

# Adding synonymous species names not picked up by taxize
Host_Synonyms <- data.frame("IUCNName" = Hostlist, "GBIFName" = Hostlist, stringsAsFactors = FALSE)
Host_Synonyms <- Host_Synonyms %>%
  mutate(GBIFName = case_when(GBIFName == "Martes pennanti"  ~ "Pekania pennanti",
                              GBIFName == "Melogale subaurantiaca"  ~ "Melogale moschata",
                              GBIFName == "Neovison vison"   ~ "Mustela vison",
                              GBIFName == "Tragelaphus oryx" ~ "Taurotragus oryx",
                              TRUE                           ~ GBIFName))

if(!Already_Prepped) {
# Getting taxon keys ##########################################################################################

Taxon_Keys <- taxize::get_gbifid_(Host_Synonyms$GBIFName, method = "backbone")
Taxon_Keys <- lapply(Host_Synonyms$GBIFName, function(name) {
  Taxon_Keys[[name]]["GBIFName"] <- name
  return(Taxon_Keys[[name]])
})

Taxon_Keys <- Taxon_Keys %>%
  bind_rows() %>%
  filter(class == "Mammalia") %>%
  filter(status == "ACCEPTED" & matchtype == "EXACT")

# Data import #############################################################################################
warning("Need to provide GBIF credentials according to ?occ_download (under 'Authentication')")

Download_Key <- occ_download(
  pred_in("taxonKey", Taxon_Keys$usagekey),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

saveRDS(Download_Key, here::here("Data/GBIF/Download_Key"))
occ_download_wait(Download_Key)
Download_Get <- occ_download_get(Download_Key, path = here::here("Data/GBIF/"), overwrite = TRUE)
saveRDS(Download_Get, here::here("Data/GBIF/Download_Get"))

rm(Download_Key)

} else{
  Download_Get <- readRDS(here::here("Data/GBIF/Download_Get"))
  # GBIF.org (08 March 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.f24c8d
}

GBIF_Raw_Data <- occ_download_import(Download_Get, path = here::here("Data/GBIF/"))
nrow(GBIF_Raw_Data) #5121490
Host_Synonyms[!Host_Synonyms$GBIFName %in% unique(GBIF_Raw_Data$species),]
# Cervus canadensis is recorded as Cervus elaphus subsp canadensis so not actually missing data

# GBIF_Raw_Plots_00 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Raw_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_Raw_Plots_00) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_Raw_Plots_00.pdf'), width = 10, height = 7)
# GBIF_Raw_Plots_00
# dev.off()

rm(Download_Get)
# Cleaning data ###############################################################################################

GBIF_Data <- GBIF_Raw_Data %>% 
  filter(countryCode != "" & countryCode != "ZZ") %>% 
  mutate(countryCode = case_when(countryCode == "XK" ~ "RS",
                                 TRUE ~ countryCode)) 
GBIF_Data$countryCode <- countrycode(GBIF_Data$countryCode, origin = "iso2c", destination = "iso3c")
nrow(GBIF_Data) #5120849

GBIF_Data <- clean_coordinates(x = GBIF_Data,
                               lon = "decimalLongitude", 
                               lat = "decimalLatitude", 
                               countries = "countryCode",
                               tests = c("capitals", "centroids", "countries", "gbif", "institutions", "zeros"), 
                               capitals_rad = 10000,
                               centroids_rad = 1000,
                               centroids_detail = "country",
                               inst_rad = 100,
                               zeros_rad = 0.5,
                               value = "clean")

nrow(GBIF_Data) #4825721
gc()

# GBIF_bor_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "basisOfRecord", range.polygon = IUCN_Native_Data)
# names(GBIF_bor_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_bor_Plots_01.pdf'), width = 10, height = 7)
# GBIF_bor_Plots_01
# dev.off()
# 
# GBIF_base_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_base_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_base_Plots_01.pdf'), width = 10, height = 7)
# GBIF_base_Plots_01
# dev.off()
#
# GBIF_epithet_Plots_01 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "infraspecificEpithet", range.polygon = IUCN_Native_Data)
# names(GBIF_epithet_Plots_01) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_epithet_Plots_01.pdf'), width = 10, height = 7)
# GBIF_epithet_Plots_01
# dev.off()


## Issues #### 

### Coordinate uncertainty ####

# Raster resolution used later in analysis is 10 arcminutes (roughly 0.167 degrees)
# This corresponds to approx 18520m at the equator
# From histograms its clear that 5000m is a commonly estimated uncertainty for coordinates
# Proportionally few records have an uncertainty over 5000m
# Need to explain the logic jump a bit but it seems to make sense to use 5000m as the max uncertainty
# Wouldn't lose too much data and it means the true coordinates at worst will only be one raster cell over
# Because of spatial autocorrelation this shouldn't have too much effect on analyses

# GBIF_Data %>%
#   pull(coordinateUncertaintyInMeters) %>%
#   hist(main = "Histogram of coordinate uncertainty, no max")
# 
# GBIF_Data %>%
#   filter(coordinateUncertaintyInMeters <= 18520) %>%
#   pull(coordinateUncertaintyInMeters) %>%
#   hist(main = "Histogram of coordinate uncertainty, max 18520m")

### Coordinate precision ####

# Raster resolution used later in analysis is 10 arcminutes (roughly 0.167 degrees)
# I could use similar reasoning as I did for coordinate uncertainty and apply 0.05 as a cut off
# But using 0.02 only removes a marginal amount more while ensuring greater raster cell accuracy

# GBIF_Data %>%
#   pull(coordinatePrecision) %>%
#   hist(main = "Histogram of coordinate precision, no max")
# 
# GBIF_Data %>%
#   filter(coordinatePrecision <= 0.167) %>%
#   pull(coordinatePrecision) %>%
#   hist(main = "Histogram of coordinate precision, max 0.167") 
# 
# nrow(GBIF_Data[GBIF_Data$coordinatePrecision <= 0.02,]) / nrow(GBIF_Data)
# nrow(GBIF_Data[GBIF_Data$coordinatePrecision <= 0.05,]) / nrow(GBIF_Data)

### Event date ####
 ## Not sure whether to filter by date
# GBIF_Data %>%
#   filter(!is.na(eventDate)) %>%
#   pull(eventDate) %>%
#   as.Date() %>%
#   hist(breaks = "years",
#        main = "Histogram of event dates") 

### Basis of record ####

# https://data-blog.gbif.org/post/living-specimen-to-preserved-specimen-understanding-basis-of-record/

# Using Chrysocyon brachyurus as an illustrative example for specimen basis of record

# TODO: redo plots with tmap
# ggplot(data = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")) +
#   geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
#   xlab("Longitude") + ylab("Latitude") +
#   
#   geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Chrysocyon brachyurus",]), 
#                aes(x = long, y = lat, group = group),
#                colour = "black",
#                fill = NA) +
#   
#   geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
#              aes(x = decimalLongitude, y = decimalLatitude, colour = basisOfRecord)) +
#   
#   ggtitle("Chrysocyon brachyurus basis of record")

### Locality ####

# ggplot(data = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")) +
#   geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
#   xlab("Longitude") + ylab("Latitude") +
#   
#   geom_polypath(data = fortify(IUCN_Native_Data[IUCN_Native_Data$binomial == "Chrysocyon brachyurus",]), 
#                aes(x = long, y = lat, group = group),
#                colour = "black",
#                fill = NA) +
#   
#   geom_point(data = filter(GBIF_Data, species == "Chrysocyon brachyurus"),
#              aes(x = decimalLongitude, y = decimalLatitude, 
#                  colour = (species == "Chrysocyon brachyurus" & 
#                    str_detect(locality, regex("zoo", ignore_case = TRUE))))) +
#   
#   scale_colour_discrete(name = "Zoo") +
#   ggtitle("Chrysocyon brachyurus zoos")

## Filtering ####

GBIF_Data <- GBIF_Data %>%
  filter(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision <= 0.02 | is.na(coordinatePrecision)) %>%
  filter(!str_detect(basisOfRecord, "_SPECIMEN")) %>%
  filter(!str_detect(locality, regex("zoo", ignore_case = TRUE))) %>%
  filter(!basisOfRecord %in% c("MATERIAL_SAMPLE", "MATERIAL_CITATION")) %>%
  filter(!str_detect(basisOfRecord, "UNKNOWN")) %>%
  filter(infraspecificEpithet != "domesticus") %>%
  filter(collectionCode != "NationalInvasiveSpeciesDatabase")

nrow(GBIF_Data) #3888981

## Species names ####
GBIF_Data <- GBIF_Data %>% 
  mutate(species = case_when(species == "Pekania pennanti" ~ "Martes pennanti",
                             species == "Mustela vison"    ~ "Neovison vison",
                             species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                             TRUE ~ species))

# Cervus elaphus
# not enough samples have subspecies recorded to use as Cervus canadensis
GBIF_Data %>% filter(species == "Cervus elaphus") %>% pull(infraspecificEpithet) %>% as.factor() %>% summary()

# splitting them to make plots more readable
GBIF_Data <- GBIF_Data %>% 
  mutate(species = case_when(species == "Cervus elaphus" &
                               (decimalLongitude < -40 |
                               decimalLongitude > 68) ~ "Cervus canadensis",
                             TRUE ~ species))

# Canis lupus
# Lots of samples are Canis (lupus) familiaris or dingo
GBIF_Data <- GBIF_Data %>% 
  filter(!infraspecificEpithet %in% c("familiaris", "dingo")) %>% 
  filter(!str_detect(verbatimScientificName, "familiaris|dingo"))

# Capreolus capreolus
# Sample in Northern Ireland is invasive but caught by buffer
GBIF_Data <- GBIF_Data %>% 
  filter(!(species == "Capreolus capreolus" & stateProvince == "Northern Ireland"))

## Outliers ####

GBIF_Data <- sf::st_as_sf(GBIF_Data,
                          coords = c("decimalLongitude", "decimalLatitude"),
                          crs = Projection_String)

tmap_mode("plot")

# lapply(Hostlist, function(hostname, buff = 1){
#   pdf(here::here("GBIF cleaning", hostname, paste0(hostname, " outliers.pdf")),
#       width = 10, height = 7)
#   print(plot_native_gbif(hostname, buff = buff))
#   dev.off()
#   print(hostname)
# })

# Justifications
# Generally don't want to include vagrant animals migrating outside their usual range
# E.g. moose in Germany or golden jackals in Central Europe

## method: 
# If they are all within 1 buffer of the native IUCN range polygon, 
# & outside 1 buffer of the non-native I leave them be
# If there are any outside that I check the literature


### Aepyceros melampus ####
# IUCN SSC Antelope Specialist Group. 2016. Aepyceros melampus. The IUCN Red List of Threatened Species 2016: e.T550A50180828. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T550A50180828.en. Accessed on 11 March 2024.
# Introduced to Gabon and "numerous privately owned game ranches and small reserves throughout southern Africa"

### Alces alces ####
# https://rewildingeurope.com/rew-project/reintroducing-moose-to-lille-vildmose/
# Schönfeld, F. Presence of moose (Alces alces) in Southeastern Germany. Eur J Wildl Res 55, 449–453 (2009). https://doi.org/10.1007/s10344-009-0272-5
# Hundertmark, K. 2016. Alces alces. The IUCN Red List of Threatened Species 2016: e.T56003281A22157381. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T56003281A22157381.en. Accessed on 11 March 2024.
# Re-introduced in Denmark, migratory in parts of Germany
# No mention of introductions in IUCN so can be liberal 
# hard to avoid migratory in Germany while also including Eastern European samples
# Because of the large amount of data available, sticking with a limited buffer of 2
# Checking original scientific name against reported name

### Antidorcas marsupialis ####
#IUCN SSC Antelope Specialist Group. 2016. Antidorcas marsupialis (errata version published in 2017). The IUCN Red List of Threatened Species 2016: e.T1676A115056763. https://dx.doi.org/10.2305/IUCN.UK.2016-3.RLTS.T1676A50181753.en. Accessed on 11 March 2024.
# "Smaller, introduced populations occur widely in extralimital areas, e.g., on private land and provincial reserves in parts of KwaZulu-Natal and the northern bushveld (East 1999).
# Relatively frequent translocations to game reserves suggests the samples in Southern South Africa are introduced
# buffer of 0.5 stays restrictive
# 

### Antilocapra americana ####
# IUCN SSC Antelope Specialist Group. 2016. Antilocapra americana (errata version published in 2017). The IUCN Red List of Threatened Species 2016: e.T1677A115056938. https://dx.doi.org/10.2305/IUCN.UK.2016-3.RLTS.T1677A50181848.en. Accessed on 11 March 2024.
# Introduced to Hawaii but now extinct there
# Material citation in Mexico is surrounded by fossil specimens, so removing material citations
# No other issues in the data are apparent so can keep remainder

### Bison bison ####
# Aune, K., Jørgensen, D. & Gates, C. 2017. Bison bison (errata version published in 2018). The IUCN Red List of Threatened Species 2017: e.T2815A123789863. https://dx.doi.org/10.2305/IUCN.UK.2017-3.RLTS.T2815A45156541.en. Accessed on 11 March 2024.
# "About 97% of the continental population is managed for private captive commercial propagation; very few of these herds are managed primarily for species conservation and none is managed in the public interest for conservation."
# Assuming that most samples are captive bred populations
# Extreme outliers in Germany and Congo are also clearly invalid
# buffer of 1 excludes the majority of samples but keep those close to defined range

### Bison bonasus ####
# Didn't see any information on introductions, but they are also kept in managed or captive herds, or as livestock
# https://www.visitnordfyn.com/nordfyn/explore/ditlevsdal-bison-farm-gdk622638
# Buffer of 1 excludes the majority of samples but keep those close to defined range

### Canis adustus ####
# Hoffmann, M. 2014. Canis adustus. The IUCN Red List of Threatened Species 2014: e.T3753A46254734. https://dx.doi.org/10.2305/IUCN.UK.2014-1.RLTS.T3753A46254734.en. Accessed on 11 March 2024.
# No mention of introductions and all outliers are close to native range
# increase buffer to 1.5 to keep all samples

### Canis aureus ####
# Hoffmann, M., Arnold, J., Duckworth, J.W., Jhala, Y., Kamler, J.F. & Krofel, M. 2018. Canis aureus (errata version published in 2020). The IUCN Red List of Threatened Species 2018: e.T118264161A163507876. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T118264161A163507876.en. Accessed on 11 March 2024.
# "expanded their range into new countries where they probably mostly occur as vagrants for now, such as Switzerland, Poland and Germany."
# To avoid vagrants, keeping it to 1

### Canis latrans ####
# Kays, R. 2018. Canis latrans (errata version published in 2020). The IUCN Red List of Threatened Species 2018: e.T3745A163508579. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T3745A163508579.en. Accessed on 11 March 2024.
# Heavily sampled and all except extreme outliers in Europe appear valid
# Increase to 5 to keep all

### Canis lupus #### 
# Keeping it to 1 to avoid records of feral domestic dogs

### Canis mesomelas ####
# single outlier seems likely to be vagrant

### Capra ibex ####
# Toïgo, C., Brambilla, A., Grignolio, S. & Pedrotti, L. 2020. Capra ibex. The IUCN Red List of Threatened Species 2020: e.T42397A161916377. https://dx.doi.org/10.2305/IUCN.UK.2020-2.RLTS.T42397A161916377.en. Accessed on 11 March 2024.
# "All current populations originate from re-introductions or introductions, except for the population in the Gran Paradiso National Park (Italy)."
# " its distribution is still fragmented because of the very long recolonization time of the species."
# Long recolonisation time means outliers are probably misrecorded, introduced, or captive populations
# buffer of 1 rules out most outliers effectively

### Capreolus capreolus ####
# Most samples are near verified home range, aside from two outliers in Ural mountains region
# Extending the buffer to 3 after removing Irish which is invasive

### Cephalophus natalensis ####
# IUCN SSC Antelope Specialist Group. 2016. Cephalophus natalensis. The IUCN Red List of Threatened Species 2016: e.T4144A50183272. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T4144A50183272.en. Accessed on 11 March 2024.
# "Sometimes considered to include Harvey's Duiker C. harveyi (e.g., Grubb and Groves 2001, Grubb 2005), but the two species are here retained as distinct"
# IUCN SSC Antelope Specialist Group. 2016. Cephalophus harveyi. The IUCN Red List of Threatened Species 2016: e.T4154A50184807. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T4154A50184807.en. Accessed on 11 March 2024.
# C. harveyi's range covers the area where there are outliers
# Keeping it at 1 because there is some overlap in the northern part of C. natalenis's range

### Cerdocyon thous ####
# Lucherini, M. 2015. Cerdocyon thous. The IUCN Red List of Threatened Species 2015: e.T4248A81266293. https://dx.doi.org/10.2305/IUCN.UK.2015-4.RLTS.T4248A81266293.en. Accessed on 11 March 2024.
# Full range is unconfirmed so extending buffer to 2

### Cervus canadensis ####
# Brook, S.M., Pluháček, J., Lorenzini, R., Lovari, S., Masseti, M., Pereladova, O. & Mattioli, S. 2018. Cervus canadensis (errata version published in 2019). The IUCN Red List of Threatened Species 2018: e.T55997823A142396828. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T55997823A142396828.en. Accessed on 11 March 2024.
# "The distribution is much more patchy and fragmented than the apparent continuity suggested by the distribution map."
# Extending the buffer to 5 to capture samples in USA that are contiguous with those within the range polygon

### Cervus elaphus ####
# Lovari, S., Lorenzini, R., Masseti, M., Pereladova, O., Carden, R.F., Brook, S.M. & Mattioli, S. 2018. Cervus elaphus (errata version published in 2019). The IUCN Red List of Threatened Species 2018: e.T55997072A142404453. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T55997072A142404453.en. Accessed on 11 March 2024.
# "It is widely but somewhat patchily distributed throughout most of continental Europe,"
# "Whether the Red Deer is native to Ireland or introduced is still under debate"
# Extending to 2 to capture missed samples in Outer Hebrides and including Ireland as part of historic range

### Cervus nippon ####
# Dhakal, T., Kim, TS., Kim, SH. et al. Distribution of sika deer (Cervus nippon) and the bioclimatic impact on their habitats in South Korea. Sci Rep 13, 19040 (2023). https://doi.org/10.1038/s41598-023-45845-2
# Introduced to South Korea so although a couple of valid samples in Taiwan are missed, buffer needs to stay close
# Keeping it at 1

### Chrysocyon brachyurus ####
# Paula, R.C. & DeMatteo, K. 2015. Chrysocyon brachyurus (errata version published in 2016). The IUCN Red List of Threatened Species 2015: e.T4819A88135664. https://dx.doi.org/10.2305/IUCN.UK.2015-4.RLTS.T4819A82316878.en. Accessed on 11 March 2024.
# "increase in the Maned Wolf distribution is likely due to both genuine change (the species’ expansion into new regions) as well as non-genuine change (mostly new or additional records from the field)"
# Estimation of its range is still developing
# Only one outlier near native range
# Extending to 2

### Conepatus chinga ####
# Emmons, L., Schiaffini, M. & Schipper, J. 2016. Conepatus chinga. The IUCN Red List of Threatened Species 2016: e.T41630A45210528. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41630A45210528.en. Accessed on 11 March 2024. 
# Emmons, L. & Helgen, K. 2016. Conepatus humboldtii. The IUCN Red List of Threatened Species 2016: e.T41631A45210677. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41631A45210677.en. Accessed on 11 March 2024.
# Cuarón, A.D., Helgen, K. & Reid, F. 2016. Conepatus semistriatus. The IUCN Red List of Threatened Species 2016: e.T41633A45210987. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41633A45210987.en. Accessed on 11 March 2024.
# Could be confused with C. humboltii in South of range or C. semistriatus in the North East
# Decreasing the buffer to 0.5 to be conservative

### Connochaetes gnou ####
# Vrahimis, S., Grobler, P., Brink, J., Viljoen, P. & Schulze, E. 2017. Connochaetes gnou. The IUCN Red List of Threatened Species 2017: e.T5228A50184962. https://dx.doi.org/10.2305/IUCN.UK.2017-2.RLTS.T5228A50184962.en. Accessed on 11 March 2024.
# IUCN SSC Antelope Specialist Group. 2016. Connochaetes taurinus (errata version published in 2020). The IUCN Red List of Threatened Species 2016: e.T5229A163322525. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T5229A163322525.en. Accessed on 11 March 2024.
# Samples in Namibia might be C. taurinus
# There are also populations introduced for hunting so keeping buffer at 1

### Connochaetes taurinus ####
# IUCN SSC Antelope Specialist Group. 2016. Connochaetes taurinus (errata version published in 2020). The IUCN Red List of Threatened Species 2016: e.T5229A163322525. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T5229A163322525.en. Accessed on 11 March 2024.
# Similar to C. gnou, there are populations introduced for hunting so keeping it at 1

### Crocuta crocuta ####
# Remaining outliers are geographically isolated, keeping at 1

### Cynictis penicillata ####
# One isolated outlier in the North

### Damaliscus lunatus ####
# IUCN SSC Antelope Specialist Group. 2016. Damaliscus lunatus. The IUCN Red List of Threatened Species 2016: e.T6235A50185422. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T6235A50185422.en. Accessed on 11 March 2024.
# Probably has introduced populations for hunting
# No good justification for extending so leaving it as is

### Equus quagga ####
# King, S.R.B. & Moehlman, P.D. 2016. Equus quagga. The IUCN Red List of Threatened Species 2016: e.T41013A45172424. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T41013A45172424.en. Accessed on 11 March 2024.
# Probably has introduced populations for hunting
# No good justification for extending so leaving it as is

### Equus zebra ####
# Gosling, L.M., Muntifering, J., Kolberg, H., Uiseb, K. & King, S.R.B. 2019. Equus zebra (amended version of 2019 assessment). The IUCN Red List of Threatened Species 2019: e.T7960A160755590. https://dx.doi.org/10.2305/IUCN.UK.2019-1.RLTS.T7960A160755590.en. Accessed on 11 March 2024.
# Probably has introduced populations for hunting
# No good justification for extending so leaving it as is

### Felis silvestris ####
# Gerngross, P., Ambarli, H., Angelici, F.M., Anile, S., Campbell, R., Ferreras de Andres, P., Gil-Sanchez, J.M., Götz, M., Jerosch, S., Mengüllüoglu, D., Monterroso, P. & Zlatanova, D. 2023. Felis silvestris (amended version of 2022 assessment). The IUCN Red List of Threatened Species 2023: e.T181049859A224982454. https://dx.doi.org/10.2305/IUCN.UK.2023-1.RLTS.T181049859A224982454.en. Accessed on 11 March 2024.
# "Several wildcat populations occurring in southern and eastern Anatolia and the Lesser Caucasus (Azerbaijan, Armenia and Iran) are F. lybica rather than F. silvestris
# Felis lybica and Felis silvestris recognised as different species in 2017 so lots of incorrectly assigned samples in south
# Extending to 1.5 to catch samples in Portugal and France

### Giraffa camelopardalis ####
# Reducing non-native to avoid catching valid samples
# Keeping at 1 to avoid catching possible vagrant in Benin

### Hippotragus niger ####
# IUCN SSC Antelope Specialist Group. 2017. Hippotragus niger. The IUCN Red List of Threatened Species 2017: e.T10170A50188654. https://dx.doi.org/10.2305/IUCN.UK.2017-2.RLTS.T10170A50188654.en. Accessed on 11 March 2024.
# "this spectacular antelope’s aesthetic appeal and its high value as a trophy animal"
# Likely introductions outside native range for hunting so keeping at 1

### Kobus ellipsiprymnus ####
# IUCN SSC Antelope Specialist Group. 2016. Kobus ellipsiprymnus. The IUCN Red List of Threatened Species 2016: e.T11035A50189324. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T11035A50189324.en. Accessed on 11 March 2024.
# Another tha's probably introduced for hunting elsewhere
# keep at 1

### Leopardus geoffroyi ####
# One outlier in North America

### Leopardus pardalis ####
# Extending to 2 to include contiguous samples without including likely vagrant in North America

### Leopardus tigrinus ####
# Payan, E. & de Oliveira, T. 2016. Leopardus tigrinus. The IUCN Red List of Threatened Species 2016: e.T54012637A50653881. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T54012637A50653881.en. Accessed on 11 March 2024.
# de Oliveira, T., Paviolo, A., Schipper, J., Bianchi, R., Payan, E. & Carvajal, S.V. 2015. Leopardus wiedii. The IUCN Red List of Threatened Species 2015: e.T11511A50654216. https://dx.doi.org/10.2305/IUCN.UK.2015-4.RLTS.T11511A50654216.en. Accessed on 11 March 2024.
# de Oliveira, T., Trigo, T., Tortato, M., Paviolo, A., Bianchi, R. & Leite-Pitman, M.R.P. 2016. Leopardus guttulus. The IUCN Red List of Threatened Species 2016: e.T54010476A54010576. https://dx.doi.org/10.2305/IUCN.UK.2016-2.RLTS.T54010476A54010576.en. Accessed on 11 March 2024.
# Overlaps with lots of other lovely little cats, especially the Margay (L. weidii)
# Southern samples are probably L. guttulus
# Restricting to 0.5 because of Margay

### Lontra canadensis ####
# Only far outliers are excluded, keeping at 1

### Lutra lutra ####
# Loy, A., Kranz, A., Oleynikov, A., Roos, A., Savage, M. & Duplaix, N. 2022. Lutra lutra (amended version of 2021 assessment). The IUCN Red List of Threatened Species 2022: e.T12419A218069689. https://dx.doi.org/10.2305/IUCN.UK.2022-2.RLTS.T12419A218069689.en. Accessed on 11 March 2024.
"The ongoing recovery in Europe is filling the gap in central Europe, as the 
species is expanding from Austria, Slovenia, Denmark, Netherlands, east Germany, 
and west France. Following this expansion, it returned to Switzerland in 2016 
and North Italy in 2011 (Pavanello et al. 2015). A gap still exists in Northwest 
and central Italy, as the southern Italian population is expanding but is still 
highly isolated from other European populations (Giovacchini et al. 2018)."
# Extending to 5 because all samples are contiguous

### Lycalopex gymnocercus ####
# Lucherini, M. 2016. Lycalopex gymnocercus. The IUCN Red List of Threatened Species 2016: e.T6928A85371194. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T6928A85371194.en. Accessed on 11 March 2024.
# Lucherini, M. 2016. Lycalopex culpaeus. The IUCN Red List of Threatened Species 2016: e.T6929A85324366. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T6929A85324366.en. Accessed on 11 March 2024.
# Outliers are likely other species of South American fox, especially L. culpaeus in Chile
# Sticking with 1 

### Lycaon pictus ####
# Woodroffe, R. & Sillero-Zubiri, C. 2020. Lycaon pictus (amended version of 2012 assessment). The IUCN Red List of Threatened Species 2020: e.T12436A166502262. https://dx.doi.org/10.2305/IUCN.UK.2020-1.RLTS.T12436A166502262.en. Accessed on 11 March 2024.
# "Land where residence was not confirmed (e.g., possible range, unknown range) was excluded."
# Extending to 2 to capture possible missing samples while avoiding outliers

### Lynx lynx ####
# Breitenmoser, U., Breitenmoser-Würsten, C., Lanz, T., von Arx, M., Antonevich, A., Bao, W. & Avgan, B. 2015. Lynx lynx (errata version published in 2017). The IUCN Red List of Threatened Species 2015: e.T12519A121707666. Accessed on 11 March 2024.
# Rodríguez, A. & Calzada, J. 2015. Lynx pardinus (errata version published in 2020). The IUCN Red List of Threatened Species 2015: e.T12520A174111773. https://dx.doi.org/10.2305/IUCN.UK.2015-2.RLTS.T12520A174111773.en. Accessed on 11 March 2024.
# "Was previously considered conspecific with Lynx lynx by some authorities,"
# "It was also absent from the Iberian Peninsula, where the smaller Iberian Lynx Lynx pardinus occurs."
# Extending to 2 to capture contiguous samples without catching Pyrenees samples which are probably historic Iberian lynx (L. pardinus)

### Lynx pardinus ####
# Rodríguez, A. & Calzada, J. 2015. Lynx pardinus (errata version published in 2020). The IUCN Red List of Threatened Species 2015: e.T12520A174111773. https://dx.doi.org/10.2305/IUCN.UK.2015-2.RLTS.T12520A174111773.en. Accessed on 11 March 2024.
# http://www.iberlince.eu/index.php/eng/iberian-lynx/distribution
# IUCN range only includes present range after historic decline
# Expanding to 3 to capture all samples because there are no 

### Lynx rufus ####
# Extending to 2 to catch contiguous samples in South

### Martes americana ####
# Helgen, K. & Reid, F. 2016. Martes americana. The IUCN Red List of Threatened Species 2016: e.T41648A45212861. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41648A45212861.en. Accessed on 11 March 2024.
# "Reintroduction projects in northern Michigan and Wisconsin have, apparently, restored a self-sustaining population in that region (Slough 1994). Reintroduction also has been attempted in New Hampshire and in various other parts of the north-western United States and south-western Canada (Nowak 2005)."
# Outliers seem to be reintroductions, extending to 3

### Martes foina ####
# Abramov, A.V., Kranz, A., Herrero, J., Choudhury, A. & Maran, T. 2016. Martes foina. The IUCN Red List of Threatened Species 2016: e.T29672A45202514. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T29672A45202514.en. Accessed on 11 March 2024.
# "The species was introduced to Ibiza, Balearic Islands (Spain) but it failed"
# Keeping at 1 to avoid Ibiza sample

### Martes martes ####
# Herrero, J., Kranz, A., Skumatov, D., Abramov, A.V., Maran, T. & Monakhov, V.G. 2016. Martes martes. The IUCN Red List of Threatened Species 2016: e.T12848A45199169. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T12848A45199169.en. Accessed on 11 March 2024.
# Extending to 2 to catch contiguous samples

### Martes melampus ####
# Reducing both to 0.3 to stop catching valid samples

### Martes pennanti ####
# Extending to 1.5 to catch contiguous samples without getting possible vagrants

### Meles meles ####
# Kranz, A., Abramov, A.V., Herrero, J. & Maran, T. 2016. Meles meles. The IUCN Red List of Threatened Species 2016: e.T29673A45203002. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T29673A45203002.en. Accessed on 12 March 2024.
# Extending to 7 to catch samples further East without catching possible M. leucurus

### Mephitis mephitis ####
# Extending to 1.5 to catch contiguous samples

### Mustela erminea ####
# Keeping at 1 to avoid possible other weasel species, especially in Caucasus mountains

### Mustela lutreola ####
# Keeping at 1 to avoid outliers

### Mustela putorius ####
# Skumatov, D., Abramov, A.V., Herrero, J., Kitchener, A., Maran, T., Kranz, A., Sándor, A., Saveljev, A., Savour�-Soubelet, A., Guinot-Ghestem, M., Zuberogoitia, I., Birks, J.D.S., Weber, A., Melisch, R. & Ruette, S. 2016. Mustela putorius. The IUCN Red List of Threatened Species 2016: e.T41658A45214384. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41658A45214384.en. Accessed on 12 March 2024.
# "Much information published under the name M. putorius refers specifically to M. furo"
# Reducing to 0.5 to avoid M. furo (feral or domesticated ferrets), especially in UK

### Neovison vison ####
# Keeping at 1 to avoid outliers

### Nyctereutes procyonoides ####
# Reducing both non-native and native buffer to 0.5 to avoid overlap 
# And because majority of non-native is far and excluded from native polygons

### Odocoileus hemionus ####
# https://www.adfg.alaska.gov/index.cfm?adfg=deer.printerfriendly
# Introduced to Kodiak
# Catomeris, C. (2018). The effects of introduced Sitka black-tailed deer (Odocoileus hemionus sitchensis) on plant and soil microbial communities on Haida Gwaii and silvicultural tools to improve western redcedar survival (T). University of British Columbia. Retrieved from https://open.library.ubc.ca/collections/ubctheses/24/items/1.0366156
# Introduced to Haida Gwaii
# Keeping at 1 to avoid introduced island populations

### Odocoileus virginianus ####
# increasing to 3 to catch contiguous samples without getting possible vagrants

### Ourebia ourebi ####
# increasing to 1.5 to catch a contiguous sample

### Ovibos moschatus ####
# Keeping at 1 because all valid samples appear to be captured

### Ovis canadensis ####
# Reducing non-native to 0.5 to reduce number of valid samples removed

### Ovis dalli ####
# leaving at 1 to avoid vagrant and outlier

### Ozotoceros bezoarticus ####
# Increasing to 3 to catch all because of scarcity of samples

### Panthera leo ####
# All valid samples are captured by 1

### Panthera onca ####
# Increasing to 2 to catch contiguous samples

### Panthera pardus ####
# All valid samples caught by 1

### Pecari tajacu ####
# Gongora, J., Reyna-Hurtado, R., Beck, H., Taber, A., Altrichter, M. & Keuroghlian, A. 2011. Pecari tajacu. The IUCN Red List of Threatened Species 2011: e.T41777A10562361. https://dx.doi.org/10.2305/IUCN.UK.2011-2.RLTS.T41777A10562361.en. Accessed on 12 March 2024.
# Their range has recently expanded northward in the southwestern United States (Albert et al. 2004) including into Oklahoma adjacent to Texas (Stangl and Dalquest 1990).
# Extending to 3 to catch range expansion

### Pelea capreolus ####
# Increased to 2 to catch contiguous samples

### Philantomba monticola ####
# IUCN SSC Antelope Specialist Group. 2016. Philantomba walteri. The IUCN Red List of Threatened Species 2016: e.T88418111A88418148. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T88418111A88418148.en. Accessed on 12 March 2024.
# Increasing to 2 to catch more samples while avoiding possible overlap with P. walteri

### Procyon lotor ####
# Timm, R., Cuarón, A.D., Reid, F., Helgen, K. & González-Maya, J.F. 2016. Procyon lotor. The IUCN Red List of Threatened Species 2016: e.T41686A45216638. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T41686A45216638.en. Accessed on 12 March 2024.
# Southern limits for the species are still not clear with the current most comprehensive revision considering Panama as the limit
# Extending to 5 to catch Southern samples

### Puma concolor ####
# Increasing to 2 to catch contiguous samples while avoiding vagrants

### Rangifer tarandus ####
# Gunn, A. 2016. Rangifer tarandus. The IUCN Red List of Threatened Species 2016: e.T29742A22167140. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T29742A22167140.en. Accessed on 12 March 2024.
# avoiding samples outside of iucn range in Scandinavia because of semi-domestication
# Keeping at 1

### Raphicerus campestris ####
# Keeping at 1 to avoid possible vagrant

### Redunca arundinum ####
# Keeping at 1 to avoid likely hunting populations and outliers

### Redunca fulvorufula ####
# Increasing to 2 to catch contiguous samples

### Rupicapra pyrenaica ####
# Keeping at 1 because of outliers

### Rupicapra rupicapra ####
# Most non-native are far away so buffer on non-native doesn't matter so much
# Except in Carpathian mountains
# Reducing non-native to 0.2 to avoid overlap 
# Reducing native to 0.5 to avoid overlap while still capturing spread around French native range

### Spilogale gracilis ####
# Extending to 2 to catch contiguous samples

### Sus scrofa ####
# Keeping it close with 1 to avoid feral, domestic and introduced hunting populations
# Breadth of samples means I'm not too worried about missing ones in Norway
# Maybe I should be because it's their northern range limit

### Syncerus caffer ####
# Reducing non-native to 0.5 to avoid catching valid samples
# Keeping native at 1 to avoid introduced or vagrant

### Taxidea taxus ####
# Increasing to 3 to catch contiguous samples 

### Tragelaphus angasii ####
# IUCN SSC Antelope Specialist Group. 2020. Tragelaphus strepsiceros (amended version of 2016 assessment). The IUCN Red List of Threatened Species 2020: e.T22054A166487759. https://dx.doi.org/10.2305/IUCN.UK.2020-1.RLTS.T22054A166487759.en. Accessed on 12 March 2024.
# Could be mistaken with T. strepsiceros which has overlapping range
# There appear to be clusters of samples just south of the IUCN defined range limits in Pretoria and Durban
# Extending to 3 to capture Pretoria and Durban samples

### Tragelaphus oryx ####
# Extending to 2 to catch contiguous samples 

### Tragelaphus scriptus ####
# Already captures apparently valid samples 

### Tragelaphus strepsiceros ####
# Keeping at 1 to avoid likely hunting introductions

### Urocyon cinereoargenteus ####
# Extending to 1.5 to catch coniguous samples 

### Ursus americanus ####
# Keeping at 1 to avoid likely vagrant

### Ursus arctos ####
# Extending to 1.5 to catch contiguous samples while avoiding vagrants/ outliers

### Ursus maritimus ####
# Keeping at 1 to avoid likely vagrant

### Vulpes lagopus ####
# Keeping at 1 to avoid outliers and vagrants

### Vulpes macrotis ####
# Extending to 2 to catch contiguous samples

### Vulpes velox ####
# Cypher, B. & List, R. 2014. Vulpes macrotis. The IUCN Red List of Threatened Species 2014: e.T41587A62259374. https://dx.doi.org/10.2305/IUCN.UK.2014-3.RLTS.T41587A62259374.en. Accessed on 12 March 2024.
# Roemer, G., Cypher, B. & List, R. 2016. Urocyon cinereoargenteus. The IUCN Red List of Threatened Species 2016: e.T22780A46178068. https://dx.doi.org/10.2305/IUCN.UK.2016-1.RLTS.T22780A46178068.en. Accessed on 12 March 2024.
# Keeping at 1 to avoid potential confusion with V. macrotis or Urocyon cinereoargenteus

### Vulpes vulpes ####
# Extending native to 3 to catch contiguous samples
# Overlapping non-native are already removed within 1

Native_DF <- Native_DF %>% 
  mutate(gbif_buffer = case_when(sci_name %in% c("Martes melampus")         & keep  ~ 0.3,
                                 
                                 sci_name %in% c("Antidorcas marsupialis",
                                                 "Conepatus chinga",
                                                 "Leopardus tigrinus",
                                                 "Mustela putorius",
                                                 "Nyctereutes procyonoides",
                                                 "Rupicapra rupicapra")     & keep  ~ 0.5,
                                 
                                 sci_name %in% c("Acinonyx jubatus",
                                                 "Aepyceros melampus",
                                                 "Antilocapra americana", 
                                                 "Bison bison",
                                                 "Bison bonasus",
                                                 "Blastocerus dichotomus",
                                                 "Canis aureus",
                                                 "Canis lupus",
                                                 "Canis mesomelas",
                                                 "Capra ibex",
                                                 "Capra pyrenaica",
                                                 "Cephalophus natalensis",
                                                 "Cervus nippon",
                                                 "Connochaetes gnou",
                                                 "Connochaetes taurinus",
                                                 "Crocuta crocuta",
                                                 "Cynictis penicillata",
                                                 "Damaliscus lunatus",
                                                 "Equus quagga",
                                                 "Equus zebra",
                                                 "Giraffa camelopardalis",
                                                 "Hippotragus niger",
                                                 "Kobus ellipsiprymnus",
                                                 "Leopardus geoffroyi",
                                                 "Lontra canadensis",
                                                 "Lycalopex gymnocercus",
                                                 "Lynx canadensis",
                                                 "Lynx lynx",
                                                 "Martes foina", 
                                                 "Mustela erminea",
                                                 "Mustela lutreola",
                                                 "Neovison vison",
                                                 "Odocoileus hemionus",
                                                 "Ovibos moschatus",
                                                 "Ovis canadensis",
                                                 "Ovis dalli",
                                                 "Panthera leo",
                                                 "Panthera pardus",
                                                 "Rangifer tarandus",
                                                 "Raphicerus campestris",
                                                 "Redunca arundinum",
                                                 "Rupicapra pyrenaica",
                                                 "Sus scrofa",
                                                 "Sylvicapra grimmia",
                                                 "Tragelaphus scriptus",
                                                 "Tragelaphus strepsiceros",
                                                 "Urocyon littoralis",
                                                 "Ursus americanus",
                                                 "Ursus maritimus",
                                                 "Vulpes lagopus",
                                                 "Vulpes velox")            & keep  ~ 1,
                                 
                                 sci_name %in% c("Canis adustus",
                                                 "Felis silvestris",
                                                 "Martes pennanti",
                                                 "Mephitis mephitis",
                                                 "Ourebia ourebi",
                                                 "Syncerus caffer",
                                                 "Urocyon cinereoargenteus",
                                                 "Ursus arctos")            & keep  ~ 1.5,
                                 
                                 sci_name %in% c("Alces alces",
                                                 "Cerdocyon thous",
                                                 "Cervus elaphus",
                                                 "Chrysocyon brachyurus",
                                                 "Leopardus pardalis",
                                                 "Lycaon pictus",
                                                 "Lynx rufus",
                                                 "Martes martes",
                                                 "Panthera onca",
                                                 "Pelea capreolus",
                                                 "Philantomba monticola",
                                                 "Puma concolor",
                                                 "Redunca fulvorufula",
                                                 "Spilogale gracilis",
                                                 "Tragelaphus oryx",
                                                 "Vulpes macrotis")         & keep  ~ 2, 
                                 
                                 sci_name %in% c("Capreolus capreolus",
                                                 "Lynx pardinus",
                                                 "Martes americana",
                                                 "Odocoileus virginianus",
                                                 "Ozotoceros bezoarticus",
                                                 "Pecari tajacu",
                                                 "Taxidea taxus",
                                                 "Tragelaphus angasii",
                                                 "Vulpes vulpes")           & keep  ~ 3, 
                                 
                                 sci_name %in% c("Canis latrans", 
                                                 "Cervus canadensis",
                                                 "Lutra lutra",
                                                 "Procyon lotor")           & keep  ~ 5,
                                 
                                 sci_name %in% c("Meles meles")             & keep  ~ 7, 
                                 
                                 #______ Non-native _____
                                 sci_name %in% c("Rupicapra rupicapra",
                                                 "Ursus arctos")           & !keep ~ 0.1,
                                 
                                 sci_name %in% c("Martes melampus")         & !keep ~ 0.3,
                                 
                                 sci_name %in% c("Capra ibex",
                                                 "Giraffa camelopardalis",
                                                 "Nyctereutes procyonoides",
                                                 "Ovis canadensis",
                                                 "Syncerus caffer")         & !keep ~ 0.5,
                                 
                                 sci_name %in% c("Aepyceros melampus",
                                                 "Antilocapra americana",
                                                 "Martes martes",
                                                 "Mustela erminea",
                                                 "Rangifer tarandus",
                                                 "Vulpes vulpes")           & !keep ~ 1,
                                 
                                 sci_name %in% c("Ovibos moschatus")        & !keep ~ 2,
                                 
                                 TRUE ~ NA_real_))

GBIF_Data_list <- lapply(Hostlist, function(host, point.data = GBIF_Data, 
                                            range.polygon = IUCN_Data){
  print(host)
  
  if(host != "Canis latrans"){
    # Filter Native DF to species 
    Native_DF <- Native_DF %>% 
      filter(sci_name == host) %>% 
      select(sci_name, gbif_buffer, keep) %>% 
      distinct()
    
    # filter points to species
    point.data <- point.data[point.data$species == host, ]
    
    # First remove points
    buffer.remove <- Native_DF %>% filter(!keep) %>% pull(gbif_buffer)
    
    range.remove <- range.polygon %>% 
      dplyr::filter(sci_name == host & !keep) %>% 
      sf::st_buffer(dist = buffer.remove)
    
    point.data <- point.data[lengths(st_intersects(point.data, range.remove)) == 0,]
    
    rm(range.remove)
    
    # then retain points
    buffer.retain <- Native_DF %>% filter(keep) %>% pull(gbif_buffer)
    
    range.retain <- range.polygon %>% 
      dplyr::filter(sci_name == host & keep) %>% 
      sf::st_buffer(dist = buffer.retain)
    
    point.data <- point.data[lengths(st_intersects(point.data, range.retain)) != 0,]
    
    point.data <- cbind(sf::st_drop_geometry(point.data), 
                        data.frame(sf::st_coordinates(point.data))) %>% 
      dplyr::rename(decimalLongitude = X, decimalLatitude = Y)
    
    rm(range.retain)
    
    # Write data
    write.csv(point.data, 
              file = here::here("GBIF cleaning", host, paste0(host, " GBIF data.csv")), 
              row.names = FALSE)
    
    gc()
    return(host)
    } else {
      # Canis latrans crashes everything on my laptop and 
      # the purpose of the 5 buffer was to catch all but European
      point.data <- point.data %>% 
        filter(species == host & countryCode != "DEU")
      
      point.data <- cbind(sf::st_drop_geometry(point.data), 
                          data.frame(sf::st_coordinates(point.data))) %>% 
        dplyr::rename(decimalLongitude = X, decimalLatitude = Y)
      
      write.csv(point.data, 
                file = here::here("GBIF cleaning", host, paste0(host, " GBIF data.csv")), 
                row.names = FALSE)
      
    }
})

GBIF_Data_list <- lapply(Hostlist, function(host){
  out <- read.csv(here::here("GBIF cleaning", host, paste0(host, " GBIF data.csv")))
  return(out)
})

GBIF_Data <- do.call(rbind, GBIF_Data_list)


write.csv(Native_DF, file = here::here("Data/Data back ups/Native_DF_02.csv"), row.names = FALSE)

GBIF_Data_test2 <- read.csv(here::here("Data/Data back ups/GBIF_Data_test_02.csv"), header = TRUE)




# GBIF_Base_Plots_02 <- apply(Host_Synonyms, MARGIN = 1, FUN = gbif_plotter, dat = GBIF_Data, data_type = "base", range.polygon = IUCN_Native_Data)
# names(GBIF_Base_Plots_02) <- Host_Synonyms$IUCNName
# 
# pdf(file = here::here('GBIF Cleaning/GBIF_Base_Plots_02.pdf'), width = 10, height = 7)
# GBIF_Base_Plots_02
# dev.off()

## Subgroups ####
# TODO: May move prep steps to 02.1 if I don't make other sub-scripts
# GBIF_Data <- GBIF_Data %>%
#   mutate(infraspecificEpithet = case_when(infraspecificEpithet == "" ~ NA_character_,
#                                           TRUE ~ infraspecificEpithet)) %>%
#   mutate(subgroup = infraspecificEpithet)
# 
# GBIF_Data <- as.data.frame(GBIF_Data)
# 
# GBIF_Spatial <- SpatialPointsDataFrame(coords      = GBIF_Data[, c("decimalLongitude", "decimalLatitude")],
#                                        data        = GBIF_Data[, names(GBIF_Data)[!names(GBIF_Data) %in% c("decimalLongitude", "decimalLatitude")]], 
#                                        proj4string = CRS(proj4string(IUCN_Native_Data)))
# 
# source(here::here("02.1Subgrouping GBIF.R"))
# 
# nrow(GBIF_Subgroups) #2285035


# Plotting ####
# prep
# GMPD_Spatial <- SpatialPointsDataFrame(coords      = GMPD_Data[, c("Longitude", "Latitude")],
#                                        data        = GMPD_Data[, names(GMPD_Data)[!names(GMPD_Data) %in% c("Longitude", "Latitude")]], 
#                                        proj4string = CRS(proj4string(IUCN_Native_Data)))
# 
# # polygon legend groups:
# poly_fill <- data.frame(legend = unique(IUCN_Native_Data@data$legend),
#            fill_group = NA)
# poly_fill <- poly_fill %>%
#   mutate(fill_group = case_when(legend == "Extinct" ~ "Extinct",
#                                 legend == "Extinct & Reintroduced" ~ "Extinct",
#                                 str_detect(legend, "Extinct|Presence") ~ "Absence Likely",
#                                 legend == "Extant (resident)" ~ "Extant",
#                                 str_detect(legend, "Extant & R|Extant & I") ~ "Extant",
#                                 str_detect(legend, "Extant") ~ "Presence Likely"))
# 
# # borders and scale
# sub.colours <- c('#ea3c67', '#3cb44b', '#ffbb19', # red, green, yellow
#                  '#4388d8', '#f58231', '#42d4f4', # blue, orange, cyan
#                  '#f032e6', '#469990', # magenta, teal
#                  '#b372ff', '#c37e2e', '#f577a5') # lavendar, brown, pink
# 
# Host_Synonyms <- Host_Synonyms %>% # TODO: is this still relevant?
#   mutate(minlong = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
#                                                                            bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Longitude", "min"]}),
#          maxlong = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
#                                                                            bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Longitude", "max"]}),
#          minlat = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
#                                                                           bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Latitude", "min"]}),
#          maxlat = apply(Host_Synonyms, MARGIN = 1, function(host) {bboxer(bbox(IUCN_Native_Data[IUCN_Native_Data$binomial == host["IUCNName"],]),
#                                                                           bbox(GBIF_Spatial[GBIF_Spatial$species == host["GBIFName"],]))["Latitude", "max"]})) %>%
#   mutate(scaling = case_when(IUCNName == "Vulpes vulpes" ~ "Eurasia",
#                              maxlong - minlong > 300 ~ "Global",
#                              maxlat < 38 & maxlong < 58 & minlong > -20 ~ "Africa",
#                              str_detect(IUCNName, "Herpestes ichneumon|Hyaena hyaena|Panthera leo") ~ "Africa",
#                              str_detect(IUCNName, "Puma concolor|Odocoileus virginianus") ~ "Global", 
#                              minlat > 6 & maxlong < -12  ~ "North",
#                              IUCNName == "Leopardus tigrinus" ~ "Central",
#                              maxlat < 13 & maxlong < -32  ~ "South",
#                              maxlong < -30 ~ "Central",
#                              minlong > 91 ~ "Asia",
#                              minlong > -11  & maxlong < 117 ~ "Europe",
#                              TRUE ~ "Eurasia")) %>%
#   dplyr::select(IUCNName, GBIFName, scaling)
# 
# bboxes <- list("Global" = matrix(c(-180, -90, 180, 90), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "Eurasia" = matrix(c(-180, -90, 180, 90), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "Africa" = matrix(c(-60, -41, 120, 49), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "South" = matrix(c(-180, -67, 0, 23), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "North" = matrix(c(-180, -8, 0, 82), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "Central" = matrix(c(-180, -38, 0, 52), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "Asia" = matrix(c(0, -8, 180, 82), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))),
#                "Europe" = matrix(c(-37, -8, 143, 82), 
#                                  ncol = 2, dimnames = list(c("x", "y"), c("min", "max"))))
# 
# apply(Host_Synonyms[3,], MARGIN = 1, FUN = complete_plot, 
#       dat = GMPD_Spatial, range.dat = GBIF_Spatial, range.polygon = IUCN_Native_Data)

# try plotting in pdf/document
# position legends
# check if hatching shows on small polygons

## extend x/ylim in one direction
# all data at once, symbols differ
# put iucn and gbif next to each other


# look for mathematica trainign courses
# ask bram for mathematica code
# bes, bob o hara



