####################################################################################################
##### CALCULATE POSITION OF PARASITE SAMPLE LOCATIONS IN ENVIRONMENTAL NICHE  ######################
##### Written by: Margaret Bolton mb804(at)exeter.ac.uk and Regan Early r.early@exeter.ac.uk #######
##### Written on: ... ##############################################################################
##### Modified on: 22nd April 2022   ################################################################
####################################################################################################

##### Outputs #####
# x: The locations of the grid-cell breaks on PCA axis 1
# y: The locations of the grid-cell breaks on PCA axis 2
# hostdens.uncor: Probability density of the host occurring in each PCA environmental grid-cell, uncorrected for environmental availability
# hostdens.cor: Probability density of the host occurring in each PCA environmental grid-cell, corrected for environmental availability
# env.dens: Prevalence of the environmental conditions in each PCA environmental grid-cell
# scores.env: The positions on pca axes 1 and 2 for all env locations
# scores.occ: The positions on pca axes 1 and 2 for all host locations
# pardist.uncor: The distance and angle in PCA environmental space of parasite locations from the PCA grid-cell where the uncorrected host species density is highest, i.e. 1, the centroid of the host environmental niche.
# pardist.cor: The distance and angle in PCA environmental space of parasite locations from the PCA grid-cell where the environmentally-corrected host species density is highest, i.e. 1, the centroid of the host environmental niche.
# hostocc: The grid-cells in environmental PCA space where the species has a non-trivial probability density
# hostcent.uncor: The coordincates of the PCA grid-cell where the uncorrected host species density is highest, i.e. 1, the centroid of the host environmental niche
# hostcent.cor: The coordincates of the PCA grid-cell where the environmentally-corrected host species density is highest, i.e. 1, the centroid of the host environmental niche

##### Libraries and data #####
library(ade4)
library(adehabitatMA)
library(adehabitatHR)
library(raster)
library(tidyverse)
source("E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/CODE_GITHUB_REGAN/NICHE_POSITION/niche_functions.R")
# library(here)

### Load data
dat.wd <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/DATA/22ndApril2022"
GBIF_Data <- read.csv(paste0(dat.wd, "/GBIF_Subgroups_02.csv"), header = TRUE, stringsAsFactors = FALSE) ## Distribution data
GMPD_Data <- read.csv(paste0(dat.wd, "/GMPD_Data_01.csv"), header = TRUE, stringsAsFactors = FALSE) ## Parasite load data

### Set thresholds for removing infinitesimally small densities resulting from the kernel smoother approach. Thresholds are applied in the grid.clim function.
th.sp <- 0
th.env <- 0

### IUCN and GBIF names do not always match so make a new column containing the names that match IUCN
## case_when: LHS argument determines which values match this case. RHS provides the replacement value. TRUE indicates all the other values
## However in the version of the data on 22nd April 2022 only Pekania pennanti is present.
GBIF_Data <- GBIF_Data %>%
  mutate(sp.iucn = case_when(species == "Mustela vison"  ~ "Neovison vison",
                             species == "Taurotragus oryx" ~ "Tragelaphus oryx",
                             species == "Pekania pennanti"  ~ "Martes pennanti",
                             TRUE                           ~ species))

### Climate data at 10 arc-minutes resolution
## bioclim variables described here: https://worldclim.org/data/bioclim.html
bio.dat <- getData('worldclim', var = 'bio', res = 10) ## Great function in raster library for easily accessing geographic data - elevation, administrative boundaries, current and future climate under an array of scenarios

## Isothermality, max temp, min temp, precip wettest quarter, precip driest quarter. 
## See Martins in submission for rationale on isothermality
bio1 <- subset(bio.dat, c(3, 5, 6, 13, 14)) 
bio1$bio5 <- bio1$bio5 / 10 ## Temp variables must be divided by 10. https://worldclim.org/data/v1.4/formats.html
bio1$bio6 <- bio1$bio6 / 10
bio1$bio3 <- bio1$bio3 / 100 ## Isothermality is multiplied by 100 https://worldclim.org/data/bioclim.html

## max temp, min temp, total precip. As used in Estrada papers, Early & Sax....
bio2 <- subset(bio.dat, c(5, 6, 12))
bio2$bio5 <- bio2$bio5 / 10 ## Temp variables must be divided by 10. https://worldclim.org/data/v1.4/formats.html
bio2$bio6 <- bio2$bio6 / 10

## Eight variables used in Petitpierre, Early & Sax 2014.
## Growing Degree Days above 5C
r <- getData('worldclim', var = 'tmean', res = 10)
r <- r/10 ## Temp variables must be divided by 10. https://worldclim.org/data/v1.4/formats.html
gdd5 <- calc(r, fun=function(x){ifelse(x > 5, x*30, 0)})
gdd5 <- sum(gdd5)

## Ratio of actual to potential evapotranspiration. In Early & Sax 2014 when doing the niche calculations I seemed to use data from E:\GIS_DATA\CLIMATE\1961_1990\WORLD\50km_GRID\AFE_AND_WORLD, where there is a variable aetpet, which is at 10 arc-minutes.
## This came from the folder E:\GIS_DATA\CLIMATE\1961_1990\AET_PET.
## See metadata in word doc here: E:\GIS_DATA\CLIMATE\1961_1990\AET_PET\AET
## Original was at 5 arcminutes
aetpet <- raster("E:/GIS_DATA/CLIMATE/1961_1990/AET_PET/aetpet")

## Potential evapotranspiration. In Early & Sax 2014 I seem to have used data from E:\GIS_DATA\CLIMATE\1961_1990\AET_PET.
## See evapo.txt files there for more info.
## The pet data seem to be modified from 50km resolution to 10 minute.
pet <- raster("E:/GIS_DATA/CLIMATE/1961_1990/AET_PET/pet")

bio3 <- subset(bio.dat, c(1, 4, 5, 6, 12))
bio3$bio5 <- bio3$bio5 / 10 ## Temp variables must be divided by 10. https://worldclim.org/data/v1.4/formats.html
bio3$bio6 <- bio3$bio6 / 10
bio1$bio3 <- bio1$bio3 / 100 ## Isothermality is multiplied by 100 https://worldclim.org/data/bioclim.html
bio3$gdd5 <- gdd5
bio3$aetpet <- resample(aetpet, bio3)
bio3$pet <- resample(pet, bio3)
  
### Human impact data
## Global Human Settlement Layer (already in wgs84 projection)
ghsl <- raster("E:/GIS_DATA/WORLD/HUMAN_FOOTPRINT/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2000-30ss-v1-geotiff/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2000-30ss-v1/ghsl-population-built-up-estimates-degree-urban-smod_built-30ss-2000.tif")

## Human Influence Index
# hii <- raster("E:/GIS_DATA/WORLD/HUMAN_FOOTPRINT/hii_global_geo_grid/hii_v2geo")
# test <- projectRaster(hii, ghsl) ## Project coordinate system to WGS84. This makes hii and ghsl align. Another approach is crs=CRS("+init=epsg:4326")
# writeRaster(test, "E:/GIS_DATA/WORLD/HUMAN_FOOTPRINT/hii_global_geo_grid/hii_v2geo_projghsl.tif")
hii <- raster("E:/GIS_DATA/WORLD/HUMAN_FOOTPRINT/hii_global_geo_grid/hii_v2geo_projghsl.tif")

anthro <- stack(ghsl, hii)
names(anthro) <- c("ghsl", "hii")

##### PCA and environmental niche #####
gmpd.cols <- c("ParasiteCorrectedName", 
               "HostsSampled", 
               "Longitude", 
               "Latitude", 
               "Prevalence")#, 
               # "EquatorwardsProp",
               # "EquatorwardsDist",
               # "MedianDist",
               # "MedianProp")

i <- 2 ## The name of the environmental variable set
# for (i in 1:3) {
raw <- get(paste0("bio",i))

### extract env data and calculate PCA space
env.xy <- xyFromCell(raw, 1:ncell(raw))
env <- as.data.frame(na.omit(cbind(env.xy, raster::extract(raw, env.xy)))) 
xVar <- c(3:ncol(env))
rm(env.xy)

pca <- dudi.pca(env[xVar], center = T, scale = T, scannf = F, nf = 2)

### Calculate position in env niche using kernel smoothers applied to the distribution and PCA data
## Example for test runs
# species.name <- "Martes martes"
# spat.dat <- GBIF_Data
# env.dat <- raw
# par.dat <- GMPD_Data

ks.test <- function(species.name, spat.dat, env.dat, par.dat) {
  
  pres <- spat.dat %>%
    filter(sp.iucn == species.name) %>%
    dplyr::select(c("decimalLongitude", "decimalLatitude"))
  
  ## Filter records to one per 10 minute grid-cell. (partially) account for sampling bias by gridding species occurence data and extracting filled cells
  r10 <- raster(resolution = 10/60)
  r.sp <- rasterize(pres, r10, fun = 'count')
  occ <- Which(r.sp, cells = TRUE) ## The cell numbers of the grid-cells that contain records
  
  if(length(occ)>=5) {
    ## extract env variables for points where species occurs
    occ.xy <- xyFromCell(r.sp, occ)
    occ.xy <- as.data.frame(na.omit(cbind(occ.xy, raster::extract(x = env.dat, y = occ.xy))))
    
    nvar <- length(xVar)     
    R <- 100 ## Number of rows and columns in the kernel smoother matrix
    
    #rbind occ data to glob data and add binary weighting column/vector
    #row.w.env <- c(rep(1, 1:nrow(clim)), rep(0, (nrow(clim) + 1):(nrow(clim) + nrow(occ.sp))))   #1 for background clim, zero for occ clim
    #Is this necessary? 
    
    ### Identify the positions of mammal locations and parasite sample points in env space
    scores.env <- suprow(pca, env[, xVar])$lisup ## The pca scores for all env grid-cells
    scores.occ <- suprow(pca, occ.xy[, xVar])$lisup ## The pca scores for grid-cells containing the species
    
    gmpd.sp <- par.dat %>%
      filter(HostCorrectedName == species.name) %>%
      dplyr::select(gmpd.cols)
    gmpd.xy <- as.data.frame(na.omit(cbind(gmpd.sp, raster::extract(env.dat, gmpd.sp[, c("Longitude", "Latitude")])))) ## Locations of parasite sampling
    scores.gmpd <- cbind(suprow(pca, gmpd.xy[, names(env[, xVar])])$lisup, gmpd.xy) ## pca scores for parasite sampling locations
    
    ### Calculate the environmental prevalence, the host species occurrence densities, the distance of the parasite sampling locations from the centre of the species occurrence.
    out <- grid.clim(scores.env, scores.occ, scores.gmpd, R, th.sp, th.env) ## grid.clim is written in niche_functions.R, sourced above.
    
    ## Obtain HII and GHSL for the host points. One value per 30arc-second raster cell only. 
    host <- pres
    coordinates(host) <- ~ decimalLongitude + decimalLatitude ## Make geographic occurrences into spatial points dataframe
    proj4string(host) <- CRS("+init=epsg:4326") ## Set coordinate system to be WGS84
    host.anthro <- unique(as.data.frame(na.omit(cbind(coordinates(host), raster::extract(x=anthro, y=host))))) ## Extract the anthro data at the host presence locations, retaining each raster grid-cell only once
    colnames(host.anthro)[3:4] <- names(anthro)
    
    ## Obtain HII and GHSL for the parasite points. One value per 30arc-second raster cell only. 
    pste <- gmpd.xy
    coordinates(pste) <- ~ Longitude + Latitude ## Make spatial points dataframe
    proj4string(pste) <- CRS("+init=epsg:4326") ## Set coordinate system to be WGS84
    pste.anthro <- as.data.frame(cbind(pste$ParasiteCorrectedName, coordinates(pste), raster::extract(x=anthro, y=pste))) ## Extract the anthro data at the host presence locations, retaining each raster grid-cell only once
    pste.anthro <- unique(na.omit(pste.anthro))
    colnames(pste.anthro)[1] <- "ParasiteCorrectedName"
    
    ## Amalgamate and return data
    out$host.anthro <- host.anthro; out$pste.anthro <- pste.anthro
  } else {out <- "insufficient host distribution data"}
  return(out)
}

### Apply the niche position function to all host species
Hostlist <- unique(GMPD_Data$HostCorrectedName)
Hostlist <- Hostlist[Hostlist %in% GBIF_Data$sp.iucn] ## currently 109 species, should be 110. "Melogale subaurantiaca" missing from GBIF_Data
options(error=recover) ## Can identify the host name where the function failed.
niche.pos <- lapply(Hostlist, ks.test, spat.dat = GBIF_Data, env.dat = raw, par.dat = GMPD_Data)
names(niche.pos) <- Hostlist

saveRDS(niche.pos, file=paste0(dat.wd, "/niche_pos_bio",i))

# }

### Example of error recovery:                                       
# Selection: 2 ## Then use ls() to see the objects in that frame. 
# Browse[7]> species.name
# [1] "Martes pennanti" ## probably naming - see above