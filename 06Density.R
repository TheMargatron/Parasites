# Exploring density data

# libraries and data
library(here)
library(tidyverse)
library(tmap)
library(sf)

# GMPD_Analysis_Data <- read.csv(here::here("Data/Data back ups/GMPD_Climate_Data_04.csv"), header = TRUE)
MCDB_Species <- read.csv(here::here("Data/Density/Thibault MCDB/MCDB_species.csv"), header = TRUE)
MCDB_Communities <- read.csv(here::here("Data/Density/Thibault MCDB/MCDB_communities.csv"), header = TRUE)
MCDB_Sites <- read.csv(here::here("Data/Density/Thibault MCDB/MCDB_sites.csv"), header = TRUE)

Hostlist <- unique(GMPD_IUCN_Species$HostCorrectedName)

TetraDENSITY_v1 <- read.csv(here::here("Data/Density/TetraDENSITY/TetraDENSITY_v.1.csv"), header = TRUE)
TetraDENSITY_v2_Data <- read.csv(here::here("Data/Density/TetraDENSITY/Dataset.csv"), header = TRUE)

# First check if they have the species we're interested in
# MCDB ####
MCDB_Species <- MCDB_Species %>% 
  mutate(FullSpecies = paste(Genus, Species, sep = " ")) %>% 
  filter(FullSpecies %in% Hostlist) 

nrow(MCDB_Species)

MCDB_Communities <- MCDB_Communities %>%
  right_join(MCDB_Species, by = join_by(Species_ID)) %>% 
  left_join(MCDB_Sites, by = join_by(Site_ID))

MCDB_Communities %>% 
  filter(Presence_only == 0) %>% 
  group_by(FullSpecies) %>% 
  summarise(n = n())

# They don't have many we're interested in so lets try tetraDENSITY
# TetraDENSITY v1 ####
# V1 from 
# https://doi.org/10.1111/geb.12756
# https://figshare.com/articles/TetraDENSITY_Population_Density_dataset/5371633
Hostlist %in% paste(TetraDENSITY_v1$Genus, TetraDENSITY_v1$Species, sep = " ") %>% summary()

TetraDENSITY_v1 %>% 
  mutate(FullSpecies = paste(Genus, Species, sep = " ")) %>% 
  filter(FullSpecies %in% Hostlist) %>% 
  nrow()

# V2 from
# https://doi.org/10.1111/geb.13476
# https://figshare.com/articles/dataset/Population_density_estimates_for_terrestrial_mammal_species/19087037/2
Hostlist %in% TetraDENSITY_v2_Data$AcceptedName_COL %>% summary()

TetraDENSITY_v2_Data <- TetraDENSITY_v2_Data %>% 
  filter(!is.na(X) & !is.na(Y))

TetraDENSITY_v2_Data %>% 
  filter(AcceptedName_COL %in% Hostlist) %>% 
  nrow()

# V2 has the most species so let's check if they have enough data for our analysis

TetraDENSITY_v2_Data %>% 
  dplyr::filter(AcceptedName_COL %in% Hostlist) %>% 
  group_by(AcceptedName_COL) %>% 
  summarise(n = n()) %>% 
  print(n = nrow(.)) %>% 
  filter(n > 150)

# plot them all to have a proper look
tmap_mode("plot")

Hostlist_Tetra2 <- Hostlist[Hostlist %in% TetraDENSITY_v2_Data$AcceptedName_COL]

temp_map <- function(current_species){
  tm_shape(World) +
    tm_polygons() +
  TetraDENSITY_v2_Data %>% 
    filter(AcceptedName_COL == current_species) %>% 
    sf::st_as_sf(coords = c("X", "Y"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
    tm_shape() + 
    tm_dots(col = "Density_km", size = 0.5, palette = "viridis") +
    tm_shape(sf::st_as_sf(GMPD_IUCN_Species[GMPD_IUCN_Species$HostCorrectedName == current_species,],
                          coords = c("Longitude", "Latitude"),
                          crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) +
    tm_dots(size = 0.25, col = "deeppink") +
    tm_layout(title = current_species,
              title.position = c('right', 'bottom'))
}

pdf(here::here("Figures/TetraDENSITY.pdf"), width = 10, height = 6)

lapply(Hostlist_Tetra2, temp_map)

dev.off()
