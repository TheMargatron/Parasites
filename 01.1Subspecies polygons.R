# Figuring things out ####

# going by subgroup rather than subspecies because of geographic groupings
IUCN_Native_Data@data$subgroup <- IUCN_Native_Data@data$subspecies

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(str_detect(subgroup, " ssp. ") ~ str_split(subgroup, "ssp. ", simplify = TRUE)[,2],
                              TRUE ~ subgroup))

GMPD_Data$subgroup <- GMPD_Data$HostReportedSubspecies

GMPD_Spatial <- SpatialPointsDataFrame(coords      = GMPD_Data[, c("Longitude", "Latitude")],
                                       data        = GMPD_Data[, names(GMPD_Data)[!names(GMPD_Data) %in% c("Longitude", "Latitude")]], 
                                       proj4string = CRS(proj4string(IUCN_Native_Data)))

## Acinonyx jubatus ####
# There are multiple subspecies that are quite geographically separated, 
# but GMPD data only occupies the range of A. j. jubatus
# therefore want to remove all polygons representing other subspecies. 
# I've removed polygons based on wiki distribution map
# https://en.wikipedia.org/wiki/Cheetah#/media/File:Acinonyx_jubatus_subspecies_range_IUCN_2015.png

curr.species <- "Acinonyx jubatus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# upper
clip_poly <- matrix(c(33.0, 0,
                  -10, 0,
                  -10, 40,
                  60.0, 40,
                  60.0, 3.7,
                  45.0, 3.7,
                  35.5, 5,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- gUnion(clip_poly, sprp[13, ])

sprp_jubatus <- sprp - clip_poly
sprp_jubatus@data$subgroup <- "jubatus"
sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_try <- raster::bind(sprp_jubatus, sprp_otherssp)

tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "jubatus"

## Aepyceros melampus ####
# Two subspecies, only the common impala (subsp. melampus) is well represented in GMPD 
# Also black-faced impala seems geographically distinct 
# https://en.wikipedia.org/wiki/Impala#/media/File:Aepyceros_melampus.svg
# IUCN: "In Namibia, the Black-faced Impala is naturally confined to the Kaokoland in the north-west, and neighbouring south-western Angola"

curr.species <- "Aepyceros melampus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

## Alcelaphus buselaphus ####
# 8 subgroup, GMPD data appears to represent major and cokii, and they are geographically distinct
# https://en.wikipedia.org/wiki/Hartebeest#/media/File:Alcelaphus_recent.png
# There is only one sample location for each subspecies so I'll have to drop this one :(

curr.species <- "Alcelaphus buselaphus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

unique(sp.gmpd.points@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Alces alces ####
# There are multiple subspecies, GMPD represents shirasi, gigas, andersoni, and americana in North America
# and alces in Europe/Western Russia, but not buturlini, cameloides, or pfizenmayeri (east of Yenisei river)
# Splitting European polygon around Yenisei should do it. 
# Yenisei data source: https://doi.org/10.1016/j.dib.2018.09.016

curr.species <- "Alces alces"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

# Used combination of Yenisei and Angara because it matched where holes were in polygon, split generously, and crossed the full polygon width
YA_Temp <- River_Data50[(River_Data50$name == "Yenisey" | River_Data50$name == "Angara"),]

YA_Temp <- disaggregate(YA_Temp)
YA_Temp$ID <- LETTERS[1:nrow(YA_Temp)]
tm_shape(YA_Temp) + tm_lines("ID", lwd = 2)
drop_ID <- "A|G|I|K"
tm_shape(YA_Temp[str_detect(YA_Temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)

YA_Temp <- YA_Temp[-grep(drop_ID, YA_Temp$ID),] # to avoid self intersections
YA_coords <- unlist(coordinates(YA_Temp), recursive = FALSE)
names(YA_coords) <- YA_Temp$ID
YA_coords[["C"]] <- YA_coords[["C"]][nrow(YA_coords[["C"]]):1,]
YA_coords[["L"]] <- YA_coords[["L"]][nrow(YA_coords[["L"]]):1,]
YA_coords <- YA_coords[c("D","B","E","C","F","H","J","L")]
YA_coords <- do.call(rbind, YA_coords)

YA_alces_coords <- rbind(YA_coords,
                   matrix(c(YA_coords[nrow(YA_coords), 1], (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), (sprp@bbox[2,1] - 1),
                            YA_coords[1,1], (sprp@bbox[2,1] - 1),
                            YA_coords[1,1], YA_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

clip_poly <- Polygon(YA_alces_coords)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(clip_poly) + tm_polygons()

# alces subspecies
sprp_try <- sprp - clip_poly

alces_poly <- Polygon(matrix(c(0, 75,
                          110, 75,
                          110, 30,
                          0, 30,
                          0, 75), 
                        ncol = 2, byrow = TRUE))
alces_poly <- SpatialPolygons(list(Polygons(list(alces_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_alces <- raster::intersect(sprp_try, alces_poly)
sprp_alces@data$subgroup <- "alces"

# americanus subspecies
sprp_americanus <- sprp_try - alces_poly
sprp_americanus@data$subgroup <- "americana"

# russian americanus subspecies
sprp_otherssp <- raster::intersect(sprp, clip_poly)

sprp_try <- raster::bind(sprp_alces, sprp_americanus, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Antidorcas marsupilis ####
# There are three recognised subspecies (wiki) 
# I only seem to have A. m. marsupialis in GMPD data whose range lies south of Orange river (wiki) 
# In the Eastern part of this range it appears to be restricted by the Vaal rather than the Orange as wiki describes it as extending up to Kimberley which is North of Orange
# "trekbokken" historically occurred across the Orange river
# https://www.ewt.org.za/wp-content/uploads/2019/02/3.-Springbok-Antidorcas-marsupialis_LC.pdf
# However, such large scale migrations no longer occur 
"Anderson C, Schultze E, Codron D, Bissett C, Gaylard A, Child MF. 2016. A conservation
assessment of Antidorcas marsupialis. In Child MF, Roxburgh L, Do Linh San E, Raimondo D, Davies-Mostert HT, editors.
The Red List of Mammals of South Africa, Swaziland and Lesotho. South African National Biodiversity Institute and
Endangered Wildlife Trust, South Africa."
# 
"While the subspecies distinction is debated, morphometric data reveal a difference in size between 
Springbok occurring on either sides of the Orange and Vaal Rivers, which is evidence for maintaining 
subspecies status (Peters & Brink 1992). "

curr.species <- "Antidorcas marsupialis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots()

OV_Temp <- River_Data50[(River_Data50$name == "Orange"| River_Data50$name == "Vaal"),]

tm_shape(sprp) + tm_polygons() + tm_shape(OV_Temp) + tm_lines() +
  tm_shape(sp.gmpd.points) + tm_dots()

OV_Temp <- disaggregate(OV_Temp)
OV_Temp$ID <- LETTERS[1:nrow(OV_Temp)]
tm_shape(OV_Temp) + tm_lines("ID", lwd = 2)
drop_ID <- "G|A|E|B"
tm_shape(OV_Temp[str_detect(OV_Temp$ID, drop_ID, negate = TRUE),]) + tm_lines("ID", lwd = 2)
OV_Temp <- OV_Temp[-grep(drop_ID, OV_Temp$ID),] # to avoid self intersections

OV_edit <- coordinates(OV_Temp[OV_Temp$ID == "F",])[[1]][[1]]
split_point <- nearestPointOnLine(OV_edit, 
                                  tail(coordinates(OV_Temp[OV_Temp$ID == "H",])[[1]][[1]], 1L))

split_row <- which(apply(OV_edit, 1, function(x) all(x == split_point)))
OV_edit <- OV_edit[split_row:nrow(OV_edit),]
tm_shape(OV_Temp) +tm_lines() + 
  tm_shape(SpatialLines(list(Lines(Line(OV_edit), ID = "6"))))  + tm_lines( col = "red")

OV_coords <- unlist(coordinates(OV_Temp), recursive = FALSE)
names(OV_coords) <- OV_Temp$ID
OV_coords[["F"]] <- OV_edit

# Turning them round
OV_coords <- OV_coords[c("J","D","I","C","H","F")]
OV_coords <- do.call(rbind, OV_coords)

OV_coords <- rbind(OV_coords,
                   matrix(c((sprp@bbox[1,1] - 1), OV_coords[nrow(OV_coords), 2],
                            (sprp@bbox[1,1] - 1), (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), OV_coords[1,2], 
                            OV_coords[1,1], OV_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

clip_poly <- Polygon(OV_coords)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(clip_poly) + tm_polygons()

sprp_marsupialis <- sprp - clip_poly
sprp_marsupialis@data$subgroup <- "marsupialis"
sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_try <-  raster::bind(sprp_marsupialis, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "marsupialis"

## Antilocapra americana ####
# three subspecies, but not all represented by GMPD
# "Mitochondrial DNA analyses since the early 1990s support the idea of clines within a wide-ranging species rather than separate subspecies (O’Gara and Yoakum 2004)."
# but peninsularis is quite geographically distinct from other subpopulations
curr.species <- "Antilocapra americana"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") +  tm_shape(sp.gmpd.points) + tm_dots()

sprp_try <- raster::disaggregate(sprp)

sprp_try@data$subgroup <- "americana x sonoriensis"
sprp_try@data[17,"subgroup"] <- "peninsularis"

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "americana x sonoriensis"

## Axis axis ####
# Only have one sample point :(
curr.species <- "Axis axis"
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]

unique(sp.gmpd.points@coords)

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Bison bison ####
# "There are two recognized subspecies in North America: Bison bison bison and B. b. athabascae."

"A herd of hybrid plains bison x wood bison lived wild in the Yukon, Canada. The Wood bison is a 
distinct subspecies that almost became extinct in the 20th century. In an attempt to save the Plains 
bison subspecies, between 1925 and 1928, thousands of Plains bison were released into Wood Buffalo 
Park (a preserve for the Wood bison subspecies). They readily interbred and produced a 12,000 strong 
herd by 1934. The Wood bison was nearly hybridized into extinction. A small genetically pure herd was 
recovered from an isolated area in 1959 and is now being kept isolated from introduced Plains bison. 
Unfortunately, recent genetic testing seems to indicate that these supposedly pure wood bison are 
also hybridized with the plains subspecies, though the majority of the genetic makeup is that of the 
'Wood Buffalo.'" #https://en.wikipedia.org/wiki/Bovid_hybrid

curr.species <- "Bison bison"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data@data[IUCN_Native_Data$binomial == curr.species, "subgroup"] <- "bison x athabascae"
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "bison x athabascae"

## Bison bonasus ####
# "Three subspecies of European bison existed in the recent past, but only one, 
# the nominate subspecies (B. b. bonasus), survives today" (wiki)
# two subspecies are widely recognized as the Lowland Bison (Bison bonasus bonasus) and 
# the Caucasian Bison (Bison bonasus caucasicus) (Kowalczyk and Plumb 2020)" (IUCN)

# "European bison herds, scattered across Central and Eastern Europe, represent two genetic lines"
# the lowland line (Poland, Belarus, and Lithuania) and the lowland-Caucasian line (southern Poland, Russia, Ukraine and Slovakia).
# https://animaldiversity.org/accounts/Bison_bonasus/

# The distribution of the lineages is artificial because it's determined by reintroduction programmes, rather than natural dispersal
# That means that the different habitat patches they occupy will not relate to their lineage, just to their overall species distribution

curr.species <- "Bison bonasus"
IUCN_Native_Data@data[IUCN_Native_Data$binomial == curr.species, "subgroup"] <- "bonasus x caucasicus"
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "bonasus x caucasicus"

## Blastocerus dichotomus ####
# No taxonomic notes on IUCN, no description of subspecies on wiki or msotw

## Canis adustus ####
# There are seven recognized subspecies of the side-striped jackal:[2]

# L. a. adusta (West Africa to most of Angola) – Sundevall's side-striped jackal
## L. a. bweha (East Africa; Kisumu, Kenya) – Elgon side-striped jackal[15]
# L. a. centralis (Central Africa; Cameroon, near the Uham River)
# L. a. grayi (North Africa; Morocco and Tunisia)
# L. a. kaffensis (Kaffa, southwestern Ethiopia) – Kaffa side-striped jackal
## L. a. lateralis (East Africa; Kenya, Uasin Gishu Plateau, south of Gabon)
## L. a. notatus (East Africa; Kenya, Loita Plains, Rift Valley Province) – Loita side-striped jackal[15]

# Leaving this one be because I can't accurately split the polygons and I'm not 100% on the subspecies
# Kenya points could be bweha, lateralis, or notatus
# unclear what Zimbabwe point is as none are described as being in Southern African
# Also range is continuous across Africa with no isolated subpopulations

# I think the species where I can't distinguish between subspecies will still be informative for latitudinal 
# gradient, but less informative for within-range analysis.
# Might be interesting to compare result from uncleaned polygons with cleaned. 
# Split the streams again after adding this step?

curr.species <- "Canis adustus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

### Canis aureus ####
"Recent studies based on mtDNA and morphology have shown that 'Golden Jackals' in Africa are larger in size than those from 
Eurasia and are actually more closely related to the Grey Wolf Canis lupus. African animals hence represent a previously 
overlooked distinct species, the African Wolf, Canis lupaster (see Rueness et al. 2011, Gaubert et al. 2012, Koepfli et al. 
2015, Viranta et al. 2017). However, the putative presence of Golden Jackal in the Sinai Peninsula of Egypt remains unclear 
(see Gaubert et al. 2012, Viranta et al. 2017)." #IUCN

# wiki s 7 subspecies
# aureus: Middle East, Iran, Turkmenistan, Afghanistan, Pakistan and Western India
# cruesemanni: Thailand
# ecsedensis: Pannonian Basin, Central Europe
# indicus: India, Nepal, Bangladesh, Bhutan
# moreotica: Southeastern Europe, Moldova, Asia Minor and the Caucasus
# naria: Coastal South West India, Sri Lanka
# syriacus: Israel, Syria,[38] Lebanon,[62] and Jordan

# Have aureus in Iran and probably moreotica in Greece
# Could subdivide according to : https://doi.org/10.1093/mspecies/sey002
# However: 
"In Greece and Dalmatia, C. aureus is documented from the Holocene (Sommer and Benecke 2005; 
Malez 1984 in Rutkowski et al. 2015) and the ancient Mediterranean populations have persisted 
and merged with jackals coming from Asia (Fabbri et al. 2014; Rutkowski et al. 2015)."
# Indicates good dispersal ability and admixture of subspecies populations so little reason to expect they are isolated

# Following IUCN range description I'm not including the algirensis subspecies (sample in Tunisia)
curr.species <- "Canis aureus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

sp.gmpd.points[countries50[countries50$name == "Tunisia",], ]

## Canis latrans ####
# although subspecies ranges are well described in wiki, there's overlap between them 
# and I don't know exactly where to draw lines. 
# Also range is pretty contiguous (apart from tiburon island, which isn't in polygons anyway). 
# The only subspecies I'm lacking are central american ones which don't appear geographically isolated
curr.species <- "Canis latrans"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

## Canis lupus ####
# for Canis lupus rufus / Canis rufus:
"See Chambers et al. (2012) for a brief review of recent literature concerning the status of this species, 
which they considered a full species, as does this assessment."

# Researching individual subspecies from IUCN taxonomic notes 
## other general sources: 
# https://en.wikipedia.org/wiki/Subspecies_of_Canis_lupus and 
# https://en.wikipedia.org/wiki/List_of_gray_wolf_populations_by_country

curr.species <- "Canis lupus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots()

# North America
# mainland population is generally continuous apart from baileyi in Arizona/New Mexico
# arctos is geographically isolated and not represented
## remaining subspecies are lycaon, occidentalis, and nubilis
clip_poly <- countries50[countries50$continent == "North America",]
clip_poly <- terra::buffer(clip_poly, 0.3) 

sprp_america <- raster::intersect(sprp, clip_poly)
sprp_america@data$subgroup <- "lycaon x occidentalis x nubilis"

sprp_america@data <- sprp_america@data %>%
  mutate(subgroup = case_when(str_detect(island, "Greenland|Ellesmere|Banks|Melville") ~ "arctos",
                              island == "Baffin" ~ "manningi",
                              island == "Vancouver" ~ "crassodon",
                              TRUE ~ subgroup))

# baileyi
sprp_baileyi <- sprp_america[states50[states50$name %in% c("Arizona", "New Mexico"),],]
sprp_baileyi@data$subgroup <- "baileyi"  

# adding america together
sprp_america <- sprp_america - states50[states50$name %in% c("Arizona", "New Mexico"),]
sprp_america <- raster::bind(sprp_america, sprp_baileyi)

tm_shape(sprp_america) + tm_polygons("subgroup")

# Europe:
# signatus is isolated (https://doi.org/10.1111/mec.14824) (Iberian population)
# italicus is isolated (https://doi.org/10.1016/j.mambio.2017.01.005) and does not hybridise much with domestic dogs (https://doi.org/10.1007/BF03194151)

sprp_eurasian <- sprp - clip_poly
sprp_eurasian@data$subgroup <- "lupus x pallipes x chanco"

# signatus
clip_poly <- countries50[countries50$name %in% c("Spain", "Portugal"),]
clip_poly <- terra::buffer(clip_poly, 0.2)
sprp_signatus <- raster::intersect(sprp_eurasian, clip_poly)

clip_poly <- countries50[countries50$name == "Andorra",]
clip_poly <- terra::buffer(clip_poly, 1)
sprp_signatus <- sprp_signatus - clip_poly
sprp_signatus@data$subgroup <- "signatus"

# italicus
clip_poly <- countries50[countries50$name %in% c("France", "Italy", "Andorra"),]
clip_poly <- terra::buffer(clip_poly, 0.5)
clip_poly <- clip_poly - countries50[countries50$name %in% c("Switzerland", "Belgium", "Luxembourg", "Germany", "Liechtenstein", "Austria"),]
sprp_italicus <- raster::intersect(sprp_eurasian, clip_poly)
sprp_italicus <- raster::disaggregate(sprp_italicus)[1,]
sprp_italicus <- raster::bind(sprp_italicus, countries50[countries50$name == "Switzerland",])
sprp_italicus@data$subgroup <- "italicus"
sprp_italicus@data <- sprp_italicus@data[, 1:29]

# Asia:
## very few asian subspecies represented within gmpd, and those that are do not appear to be isolated
## Subpopulations are also less well documented generally 
# isolated arabs population in Arabian peninsula (Yemen, Oman, Southern Saudi Arabia)
# isolated arabs and pallipes population in Sinai peninsula, Israel, Jordan, Lebanon, Southern Syria
# isolated pallipes population in Northwest India

# arabs
clip_poly <- countries50[countries50$name %in% c("Saudi Arabia", "Bahrain", "Qatar", "United Arab Emirates", "Oman", "Yemen"),]
clip_poly <- terra::buffer(clip_poly, 0.5)
clip_poly <- clip_poly - countries50[countries50$name %in% c("Egypt", "Israel", "Jordan", "Iraq"),]
clip_poly <- clip_poly - terra::buffer(countries50[countries50$name == "Kuwait",], 0.05)
sprp_arabs <- raster::intersect(sprp_eurasian, clip_poly)
sprp_arabs@data$subgroup <- "arabs"

# adding all up 
sprp_eurasian <- sprp_eurasian - sprp_signatus - sprp_italicus - sprp_arabs
sprp_try <- raster::bind(sprp_america, sprp_eurasian, sprp_signatus, sprp_italicus, sprp_arabs)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Canis mesomelas ####
# two geographically distinct subspecies, already done :)
curr.species <- "Canis mesomelas"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

# For when I'm doing gmpd subspecies
# sp.over <- over(sp.gmpd.points, gBuffer(sprp, byid = TRUE))
# sp.gmpd$subspecies <- sp.over$subspecies
# sp.gmpd$HostCorrectedName <- paste(sp.gmpd$HostCorrectedName, sp.over$subspecies, sep = " ")

## Canis simensis ####
# Two geographically distinct subspecies and I only seem to have citernii
# https://en.wikipedia.org/wiki/Ethiopian_wolf#/media/File:Canis_simensis_subspecies_range.png
curr.species <- "Canis simensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

clip_poly = matrix(c(37, 9,
                  40, 9,
                  40, 14,
                  37, 14,
                  37, 9), 
                ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- sprp - clip_poly
sprp_try@data$subgroup <- "citernii"

sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_otherssp@data$subgroup <- "simensis"

sprp_try <- raster::bind(sprp_try, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "citernii"

## Capra ibex ####
# no reported subspecies in iucn, wiki, or msotw

## Capra pyrenaica ####
# four subspecies, two extinct. I probably only have hispanica which is geographically distinct from victoriae
# subdivided based on IUCN description of ranges
curr.species <- "Capra pyrenaica"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$temp <- letters[1:nrow(sprp_try@data)]
tm_shape(sprp_try) + tm_polygons("temp")

split_hispanica <- c("a", "b", "c", "d","e", "f", "g", "i", "m", "n")
sprp_try@data <- sprp_try@data %>%
  mutate(subgroup = case_when(temp %in% split_hispanica ~ "hispanica",
                              TRUE ~ "victoriae")) %>% 
  dplyr::select(-temp)

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "hispanica"

## Capreolus capreolus ####
# The IUCN notes the following five confirmed subspecies: italicus, garganta, capreolus, caucasicus, and coxi
# Definitely have italicus and garganta
# uncertain about capreolus but I expect so, plus it's not geographically isolated
# Don't appear to have coxi and caucasicus and they are geographically separated

curr.species <- "Capreolus capreolus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$temp <- 1:nrow(sprp_try)

split_coxi <- c(4, 11, 16, 50, 53)
split_caucasicus <- c(6)

sprp_try@data <- sprp_try@data %>%
  mutate(subgroup = case_when(temp %in% split_coxi ~ "coxi",
                              temp %in% split_caucasicus ~ "caucasicus",
                              TRUE ~ "italicus x garganta x capreolus")) %>% 
  dplyr::select(-temp)

tm_shape(sprp_try) + tm_polygons("subgroup")
sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Capricornis crispus ####
# No subspecies noted by IUCN, wiki, msotw

## Cephalophus natalensis ####
# Two subspecies have been named: C. n. natalensis and C. n. robertsi (north of the Limpopo river)
# I only seem to have natalensis

curr.species <- "Cephalophus natalensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

L_Temp <- River_Data50[(River_Data50$name == "Limpopo"),]
L_coords <- unlist(coordinates(L_Temp), recursive = FALSE)[[1]]
L_coords <- rbind(L_coords,
                matrix(c(33.53, -25.2,
                         33.53, (sprp@bbox[2,1] - 1),
                         (sprp@bbox[1,1] - 1), (sprp@bbox[2,1] - 1),
                         (sprp@bbox[1,1] - 1), L_coords[1,2],
                         L_coords[1,]), 
                       ncol = 2, byrow = TRUE))

clip_poly <- Polygon(L_coords)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_natalensis <- raster::intersect(sprp, clip_poly)
sprp_natalensis@data$subgroup <- "natalensis"

sprp_robertsi <- sprp - clip_poly
sprp_robertsi@data$subgroup <- "robertsi"

sprp_try <- raster::bind(sprp_natalensis, sprp_robertsi)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "natalensis"

## Cerdocyon thous ####
# wiki lists 5 subspecies: thous, azarae, entrerianus, aquilus, germanus
# I seem to have azarae and entrerianus and the others are geographically separate

curr.species <- "Cerdocyon thous"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)
sprp_try@data[1, "subgroup"] <- "azarae x entrerianus"

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "azarae x entrerianus"

## Cervus canadensis ####
# I only have north american subspecies: canadensis, nannodes, roosevelti
# Not possible to separate each subspecies even though there are subpopulations so just going by continent
"Here we recognise the subspecies and their distributions as follows:

C. c. canadensis – N America
C. c. alashanicus – N China
C. c. nannodes - California
C. c. roosevelti - Vancouver Island, Washington state and Oregon
C. c. sibiricus – NE Kazakhstan and N Xinjiang to S Siberia and N Mongolia
C. c. xanthopygus – SE Siberia, Russian Far East, Ussuriland, Manchuria
C. c. macneilli – Lydekker 1909 (Central and SW China (N Qinghai, Gansu, Shaanxi, W Sichuan and E Xizang)
C. c. wallichii - G. Cuvier 1823 (SW China (SE Xizang), Bhutan)"

curr.species <- "Cervus canadensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup", alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# Some samples recorded as Cervus elaphus 
elaphus_temp <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == "Cervus elaphus", ]
tm_shape(sprp) + tm_polygons("subgroup", alpha = 0) + tm_shape(elaphus_temp) + tm_dots()

clip_poly <- matrix(c(-125, 61,
                      -77,61,
                      -77,25,
                      -125,25,
                      -125,61),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_americas <- raster::intersect(sprp, clip_poly)
sprp_americas@data$subgroup <- "canadensis  x nannodes x roosevelti"

sprp_otherssp <- sprp - clip_poly
sprp_try <- raster::bind(sprp_americas, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "canadensis  x nannodes x roosevelti"

## Cervus elaphus ####
# I seem to have: elaphus, italicus, and montanus
# but not: barbarus, corsicanus, maral, brauneri
# Other sources (msotw and wiki) also record scoticus in British Isles, and atlanticus in Norway/South Sweden
"Several subspecies of Western Red Deer have been recognized with their ranges as follows:
C. e. elaphus: Ireland, Great Britain, continental Europe
C. e. barbarus: Atlas Mountains (Algeria, Tunisia)
C. e. corsicanus: Corsica (extinct, reintroduced in 1985), Sardinia
C. e. maral: Anatolia,
C. e. italicus: Italy (Ferrara)
C. e. brauneri Crimea (Russia)
C. e. montanus (syn. Carpathicus) Carpathian mountains" #IUCN

# Not sure which subspecies are around Georgia/East Turkey/Iran, or south of Finland (in Russia) from IUCN description
# But based on wiki they are maral and elaphus, respectively

"However, recent analyses call into question the veracity of the morphology-based subspecific taxonomy of C. elaphus 
more generally. Rather, three genetic lineages were clearly differentiated corresponding approximately with geographical 
factors. One lineage was distributed in central-western Europe (and also Ukraine), one from eastern Europe to the 
Middle East, and a third which corresponded to C. e. barbarus and C. e. corsicanus from North Africa and Sardinia, 
the latter two nomials are considered synonyms and were recommended considered as a subspecies"

# Therefore will group hispanicus in with elaphus as it is contiguous, 
# but keep scoticus and atlanticus separate as they're isolated.
# Though there have been European introductions to England

curr.species <- "Cervus elaphus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# British Isles (scoticus)
buffer_temp <- countries50[countries50$name %in% c("United Kingdom", "Ireland"),]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
sprp_scoticus <- raster::intersect(sprp, buffer_temp)
sprp_scoticus@data$subgroup <- "scoticus"

# North Africa (barbarus)
buffer_temp <- countries50[countries50$name == "Tunisia",]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
buffer_temp <- buffer_temp + countries50[countries50$name == "Morocco",]
sprp_barbarus <- raster::intersect(sprp, buffer_temp)
sprp_barbarus@data$subgroup <- "barbarus"

# Corsica and Sardinia (corsicanus)
sprp_corsicanus <- sprp[sprp$island %in% c("Corsica", "Sardinia"),]
sprp_corsicanus@data$subgroup <- "corsicanus"

# Norway and Sweden (atlanticus)
buffer_temp <- countries50[countries50$name %in% c("Norway", "Sweden"),]
buffer_temp <- terra::buffer(buffer_temp, 0.3)
sprp_atlanticus <- raster::intersect(sprp, buffer_temp)
sprp_atlanticus@data$subgroup <- "atlanticus"

# Italy (italicus)
clip_poly <- matrix(c(18, 40,
                      16, 39,
                      7, 45,
                      13, 45,
                      18, 40),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_italicus <- raster::intersect(sprp, clip_poly)
sprp_italicus@data$subgroup <- "italicus"

# Crimea (brauneri)
clip_poly <- matrix(c(32, 46,
                      36, 46,
                      36, 44,
                      32, 44,
                      32, 46),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_brauneri <- raster::intersect(sprp, clip_poly)
sprp_brauneri@data$subgroup <- "brauneri"

# Anatolia, Caucasus, Iran (maral)
clip_poly <- matrix(c(28.6, 40.7,
                      37, 46,
                      57, 46,
                      57, 35,
                      30, 35,
                      25.5, 39.5,
                      27.5, 40.7,
                      28.6, 40.7),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_maral <- raster::intersect(sprp, clip_poly)
sprp_maral@data$subgroup <- "maral"

# Most of Europe (elaphus)
sprp_try <- raster::bind(sprp_scoticus, sprp_barbarus, sprp_corsicanus, sprp_atlanticus, sprp_italicus, sprp_brauneri, sprp_maral)
sprp_try <- sprp_try[, 1:29]
sprp_elaphus <- sprp - sprp_try
sprp_elaphus@data$subgroup <- "elaphus"

sprp_try <- raster::bind(sprp_try, sprp_elaphus)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Cervus nippon ####
# Sources:
# https://doi.org/10.1007/s10344-005-0011-5            Groves 2005
# https://doi.org/10.1046/j.1365-294X.2001.01277.x     Goodman et al 2001
# IUCN describes species status of yesoensis and hortulorum as uncertain so maintain as one species under nippon

# splitting at the Hyogo prefecture Eastern boundary based on descriptions in sources and 
# locations of gmpd samples which have subspecies name recorded

curr.species <- "Cervus nippon"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- matrix(c(134.87, 35.8,
                      134.87, 35.66,
                      134.85, 35.59,
                      134.92, 35.54,
                      134.93, 35.51,
                      135.04, 35.54,
                      135.05, 35.41,
                      134.99, 35.39,
                      134.94, 35.41,
                      134.92, 35.31,
                      135.39, 35.1,
                      135.35, 34.96,
                      135.47, 34.92,
                      135.42, 34.82,
                      135.46, 34.74,
                      135.41, 34.69,
                      134.8, 33.8,
                      134.8, 30.9,
                      129.5, 30.9,
                      129.5, 33.3,
                      130.8, 34.5,
                      132.8, 35.7,
                      134.87, 35.8),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
tm_shape(clip_poly) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(alpha = 0)

sprp_nippon <- raster::intersect(sprp, clip_poly)
sprp_nippon@data$subgroup <- "nippon"

sprp_aplodontus <- sprp[sprp$island == "Honshu" & !is.na(sprp$island),] - clip_poly - countries50[countries50$name %in% c("Russia", "China"),]
sprp_aplodontus@data$subgroup <- "aplodontus"

sprp_yesoensis <- sprp[sprp$island == "Hokkaido" & !is.na(sprp$island),]
sprp_yesoensis@data$subgroup <- "yesoensis"

sprp_try <- raster::bind(sprp_nippon, sprp_aplodontus, sprp_yesoensis)
sprp_otherssp <- sprp - sprp_try
sprp_try <- raster::bind(sprp_otherssp, sprp_try)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Chrysocyon brachyurus ####
# no subspecies described by IUCN, wiki, or msotw

## Civettictis civetta #### 
# IUCN doesn't describe any subspecies but they are recognised on wiki and msotw
# Subspecific range descriptions lack detail, are merely type specimen locations
# Range appears continuous and hard to subdivide based on biogeographic barriers

curr.species <- "Civettictis civetta"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Conepatus chinga ####
# No subspecies described by IUCN or wiki, and range descriptions not available on msotw
# Also range appears continuous 
curr.species <- "Conepatus chinga"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Connochaetes gnou ####
# No subspecies described by IUCN, wiki, or msotw

## Connochaetes taurinus ####
# already included in species range polygons <3
curr.species <- "Connochaetes taurinus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Crocuta crocuta ####
# neither IUCN nor msotw list any subspecies and wiki states:
# "all the variation seen in the then recognised subspecies could also be found in a single population"

## Cynictis penicillata ####
# no subspecies listed by IUCN, no range descriptions on wiki, and the 12 listed on msotw have no range descriptions
# Also range appears continuous
curr.species <- "Cynictis penicillata"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Damaliscus lunatus ####
# already included in species range polygons <3
# Not sure where Zambian sample falls

curr.species <- "Damaliscus lunatus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Damaliscus pygargus ####
# One subspecies already included, pygargus is the other described by IUCN
curr.species <- "Damaliscus pygargus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & is.na(subspecies) ~ "pygargus",
                              TRUE ~ subgroup))
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]) + tm_polygons("subgroup")

## Equus grevyi ####
# No subspecies on msotw
# "However, Groves and Bell (2004) concluded that the species is indeed monotypic."

## Equus quagga ####
# No subspecies on msotw
# "The molecular data represented a genetic cline" so it's monotypic

## Equus zebra ####
# "We continue to recognize Mountain Zebra as a single species comprising two subspecies."
# polygons already labelled <3
curr.species <- "Equus zebra"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Felis silvestris and lybica ####
# "A revised taxonomy of the Felidae. The final report of the Cat Classification Task Force of the IUCN/SSC Cat Specialist Group. Cat News Special Issue 11, 80 pp."
# two subspecies (silvestris and caucasica) 
# and possibly a third which I will include because it's an island population (grampia)
# I seem to also have Felis libyca
curr.species <- "Felis silvestris"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Felis silvestris
# grampia
clip_poly <- countries50[countries50$name == "United Kingdom",]
clip_poly <- raster::buffer(clip_poly, 0.2)
sprp_grampia <- raster::intersect(sprp, clip_poly)
sprp_grampia@data$subgroup <- "grampia"

# silvestris
turkey_rangepol <- countries50[countries50$name == "Turkey",]
turkey_rangepol <- raster::disaggregate(turkey_rangepol)
clip_poly <- countries50[countries50$continent == "Europe" & !countries50$name %in% c("United Kingdom", "Ireland", "Russia"),]
clip_poly <- raster::bind(clip_poly, turkey_rangepol[3,])
clip_poly <- terra::buffer(clip_poly, 0.1)
clip_poly <- gUnion(clip_poly, terra::buffer(countries50[countries50$name %in% c("Ukraine", "Croatia"),], 0.3))

sprp_silvestris <- raster::intersect(sprp, clip_poly)
sprp_silvestris@data$subgroup <- "silvestris"

# caucasica
clip_poly <- matrix(c(36.7, 44.6,
                      37.4, 46.6,
                      49.7, 43.6,
                      49, 40.5,
                      46.5, 38.9,
                      46.2, 38.9,
                      44.7, 39.8,
                      42.7, 37.3,
                      36, 35.7,
                      26.2, 36.3,
                      26.3, 40.05,
                      36.7, 44.6),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_caucasica <- raster::intersect(sprp, clip_poly)
sprp_caucasica@data$subgroup <- "caucasica"

sprp_try <- raster::bind(sprp_grampia, sprp_silvestris, sprp_caucasica)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("Prevalence")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Felis lybica
clip_poly <- matrix(c(7.9, 43.2,
                      10, 43.2,
                      10, 38.5,
                      7.9, 38.5,
                      7.9, 43.2),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp <- sprp - sprp_try
sprp <- raster::bind(sprp, sprp_silvestris[clip_poly,])
sprp@data$binomial <- "Felis lybica"

# ornata
clip_poly <- matrix(c(42.6, 48.6,
                      107.6, 48.6,
                      107.6, 17.5,
                      62.5, 14.8, 
                      56.65, 26.6,
                      50.2, 24.4,
                      47.2, 32.1,
                      42.6, 35.5,
                      42.6, 48.6),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_ornata <- raster::intersect(sprp, clip_poly)
sprp_ornata@data$subgroup <- "ornata"

# cafra
clip_poly <- matrix(c(10,2,
                      29.5, -6,
                      40.4, -10.5,
                      40.8, -10.2,
                      42, -35,
                      10, -35,
                      10,2),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - countries50[countries50$name == "Tanzania",]
clip_poly <- terra::buffer(clip_poly, 0.01)

sprp_cafra <- raster::intersect(sprp, clip_poly)
sprp_cafra@data$subgroup <- "cafra"

# lybica
sprp_lybica <- sprp - sprp_ornata - sprp_cafra
sprp_lybica@data$subgroup <- "lybica"

sprp_try <- raster::bind(sprp_ornata, sprp_cafra, sprp_lybica)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("Prevalence")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data, sprp_try)

## Galictis cuja ####
# from wiki:
"Four subspecies are recognised:

Galictis cuja cuja – southwestern Bolivia, western Argentina, central Chile
Galictis cuja furax – southern Brazil, northeastern Argentina, Uruguay, and Paraguay
Galictis cuja huronax – south-central Bolivia, eastern Argentina
Galictis cuja luteola – extreme southern Peru, western Bolivia and northern Chile"

# However, range appears continuous and I have no knowledge of barriers between subpopulations
curr.species <- "Galictis cuja"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Genetta genetta ####
# wiki:
"G. g. genetta (Linnaeus), 1758 — Spain, Portugal and France
G. g. afra (Cuvier), 1825 — North Africa[32]
G. g. senegalensis (Fischer), 1829 — sub-Saharan Africa[33]
G. g. dongolana (Hemprich and Ehrenberg), 1832 — Arabia[34]"

# IUCN: 
"This assessment includes the South African Small-spotted Genet (Genetta felina (Thunberg, 1811))"

# I have genetta in Spain, but excluding because introduced
# I don't have afra or dongolana
# Southern African population is G. felina

curr.species <- "Genetta genetta"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia", "Libya"),]
clip_poly <- terra::buffer(clip_poly, 0.1)
sprp_afra <- raster::intersect(sprp, clip_poly)
sprp_afra@data$subgroup <- "afra"

clip_poly <- countries50[countries50$name %in% c("South Africa", "Botswana", "Namibia", "Angola", "Zimbabwe"),]
clip_poly <- terra::buffer(clip_poly, 1.7)
sprp_felina <- raster::intersect(sprp, clip_poly)
sprp_felina@data$subgroup <- "felina"

clip_poly <- countries50[countries50$name %in% c("Saudi Arabia", "Yemen", "Oman"),]
clip_poly <- terra::buffer(clip_poly, 0.2)
sprp_dongolana <- raster::intersect(sprp, clip_poly)
sprp_dongolana@data$subgroup <- "dongolana"

sprp_try <- raster::bind(sprp_afra, sprp_felina, sprp_dongolana)
sprp_senegalensis <- sprp - sprp_try
sprp_senegalensis@data$subgroup <- "senegalensis"

sprp_try <- raster::bind(sprp_try, sprp_senegalensis)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "senegalensis"
## Genetta thierryi ####
# No subspecies listed by IUCN, wiki, or msotw
# plus it has a small continuous range

## Giraffa camelopardalis ####
# already recorded <3
curr.species <- "Giraffa camelopardalis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Herpestes ichneumon ####
# none listed on IUCN, but 11 listed by msotw (location notes from wiki)

# SUBSPECIES ichneumon      from the area of the Nile River in Egypt
# SUBSPECIES angolensis     a male specimen from Quissange in Angola
# SUBSPECIES cafra          based on a description of a specimen from the Cape of Good Hope
# SUBSPECIES centralis      two specimens from Beni, Democratic Republic of the Congo
# SUBSPECIES funestus       a specimen from Naivasha in British East Africa (Kenya)
# SUBSPECIES mababiensis    a specimen from Mababe in northern Bechuanaland (Botswana)
## SUBSPECIES numidicus      two individuals from Algiers in Algeria
# SUBSPECIES parvidens      three specimens collected near the lower Congo River in Congo Free State (roughly DRC)
# SUBSPECIES sabiensis      specimen from Sabi Sand Game Reserve in Southern Africa
## SUBSPECIES sangronizi     specimen from Mogador in Morocco
## SUBSPECIES widdringtonii  specimen from Sierra Morena in Spain

# I might have nominate subspecies, but hard to tell
# Much of range is continuous, apart from Iberian Peninsula and Northwest Africa which I will exclude
# Not possible to justify any further subdivision of range based on available information

curr.species <- "Herpestes ichneumon"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

# Northwest Africa (sangronizi and numidicus) and Iberian peninsula (widdringtonii)
clip_poly <- matrix(c(-16, 44,
                      18, 44,
                      18, 19,
                      -16, 19,
                      -16, 44),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_ichneumon <- sprp - clip_poly
sprp_ichneumon@data$subgroup <- "ichneumon"

sprp_try <- raster::bind(sprp_otherssp, sprp_ichneumon)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "ichneumon"

## Hippopotamus amphibius ####
# three subspecies in msotw, 
# though it's possible tschadensis and constrictus are missing because of lack of info rather than anything else
"Genetic analyses have tested the existence of three of these putative subspecies. A study examining 
mitochondrial DNA from skin biopsies taken from 13 sampling locations, considered genetic diversity and 
structure among hippo populations across the continent. The authors found low, but significant, genetic 
differentiation among H. a. amphibius, H. a. capensis, and H. a. kiboko. Neither H. a. tschadensis nor 
H. a. constrictus has been tested.[9][10]" #(wiki)
# I likely only have capensis, described in wiki as from Zambia to South Africa.
# However, it's really hard to tell if there's admixture because their population is so subdivided
# Also subspecies ranges described in wiki have a lot of overlap

curr.species <- "Hippopotamus amphibius"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Hippotragus niger
# Only variani is geographically separate, and descriptions of ranges do not specify biogeographical boundaries
"Four subspecies are usually recognized: H. n. niger, H. n. kirkii, H. n. roosevelti and the isolated Giant Sable 
(H. n. variani) from Angola. As for many other antelope species, the validity and precise distribution of most of the 
described subspecies are uncertain. An extensive study of the geographical genetic structure of Hippotragus niger 
identified three genetic subdivisions representing a Kenya and east Tanzania clade (H. n. roosevelti), a west Tanzania 
clade (H. n. kirkii), and a southern African clade (H. n. niger) (Pitra et al. 2002)." #(IUCN)

curr.species <- "Hippotragus niger"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data[IUCN_Native_Data$binomial == curr.species & is.na(IUCN_Native_Data$subgroup), "subgroup"] <- "niger x kirkii x roosevelti"
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "niger x kirkii x roosevelti"

## Hyaena hyaena ####
# As of 2005,[3] no subspecies are recognised. (wiki, from msotw)

## Kobus ellipsiprymnus ####
# already labelled <3
curr.species <- "Kobus ellipsiprymnus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Kobus kob ####
# already done <3
curr.species <- "Kobus kob"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "kob"

## Kobus leche ####
# dooone <£
curr.species <- "Kobus leche"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Leopardus geoffroyi ####
# according to cat specialist group there's only one ssp
curr.species <- "Leopardus geoffroyi"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Leopardus pardalis ####
# cat specialist group say there's two: pardalis (north) and mitis (south)
# "borders between subspecies are speculative" and I have one right on the border
"The morphological differentiation between Central and South American forms is clear and supported partly 
by molecular data as well as a clear biogeographical barrier, the Andes."
# I'll follow their apparent division at the border between Panama and Costa Rica

curr.species <- "Leopardus pardalis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- countries50[countries50$continent == "North America" & countries50$name != "Panama",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with
tm_shape(sprp) + tm_polygons(col = "red") + tm_shape(buffer_temp) + tm_polygons(alpha = 0)

sprp_pardalis <- raster::intersect(sprp, buffer_temp)
sprp_pardalis@data$subgroup <- "pardalis"
sprp_mitis <- sprp - buffer_temp
sprp_mitis@data$subgroup <- "mitis"

sprp_try <- raster::bind(sprp_pardalis, sprp_mitis)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Leopardus tigrinus ####
# "Until then L. tigrinus is recognised as having two subspecies: tigrinis (south) and oncilla" (IUCN)
# They are geographically distinct and I only have one
curr.species <- "Leopardus tigrinus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$subgroup <- c("tigrinus", "oncilla")
sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "tigrinus"

## Lontra canadensis ####
# (wiki)
"L. c. canadensis (Schreber, 1777) – (eastern Canada, U.S., Newfoundland)
L. c. kodiacensis (Goldman, 1935) – (Kodiak Island, Alaska)
L. c. lataxina (Cuvier, 1823) – (U.S.)
L. c. mira (Goldman, 1935) – (Alaska, British Columbia)
L. c. pacifica (J. A. Allen, 1898) – (Alaska, Canada, northern U.S., south to central California, northern Nevada, and northeastern Utah)
L. c. periclyzomae (Elliot, 1905) – (Queen Charlotte Islands, British Columbia)
L. c. sonora (Rhoads, 1898) – (U.S., Mexico)"

# I seem to have canadensis, pacifica, sonora, and lataxina
# Not kodiacensis, mira, periclyzomae
# Only kodiacensis and periclyzomae appear to be isolated
# Also, GBIF distribution of points and wiki distribution map show a more limited range, excluding Hawaiian islands and most of Nunavut
# I'm going to follow suit by removing Hawaii and islands North of Hudson Bay

curr.species <- "Lontra canadensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- matrix(c(-78.5, 62.7,
                      -72.7, 62.9,
                      -63.4, 60.8,
                      -48.5, 48,
                      -79.4, 23.9,
                      -169.3, 23.9,
                      -169.3, 75.3,
                      -114.35, 68.97,
                      -113.2, 68.34,
                      -109.37, 68.37,
                      -106.75, 69.16,
                      -101.9, 68.24,
                      -95.79, 71.94,
                      -96.32, 74.27,
                      -88.5, 74.22,
                      -91.93, 71.83,
                      -87.09, 69.35,
                      -85.4, 69.95,
                      -83.1, 69.81,
                      -82.56, 69.71,
                      -81.34, 69.75, 
                      -78.5, 62.7),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- raster::intersect(sprp, clip_poly)

sprp_kodiacensis <- sprp_try[sprp_try$island %in% c("Kodiak", "Afognak", "Shuyak"),]
sprp_kodiacensis@data$subgroup <- "kodiacensis"

sprp_periclyzomae <- sprp_try[sprp_try$island %in% c("Langara", "Graham", "Moresby", "Louise"),]
sprp_periclyzomae@data$subgroup <- "periclyzomae"

sprp_canadensis <- sprp_try[!(sprp_try$island %in% c("Kodiak", "Afognak", "Shuyak", "Langara", "Graham", "Moresby", "Louise")),]
sprp_canadensis@data$subgroup <- "canadensis and other"

sprp_try <- raster::bind(sprp_kodiacensis, sprp_periclyzomae, sprp_canadensis)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "canadensis and other"

## Lutra lutra ####
# following:
# https://doi.org/10.1093/mspecies/sew011
# but too hard to split aurobrunnea, kutab, and monticolus from each other (doesn't matter anyway though)

curr.species <- "Lutra lutra"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# angustifrons
clip_poly <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly <- terra::buffer(clip_poly, 0.1)
sprp_angustifrons <- raster::intersect(sprp, clip_poly)
sprp_angustifrons@data$subgroup <- "angustifrons"

# nair
clip_poly <- countries50[countries50$name == "Sri Lanka",]
clip_poly <- terra::buffer(clip_poly, 5)
sprp_nair <- raster::intersect(sprp, clip_poly)
sprp_nair@data$subgroup <- "nair"
#tm_shape(clip_poly) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(alpha = 0)

# barang
clip_poly <- matrix(c(97.4, 23.8,
                      98.2, 24.3, 
                      108, 24.6,
                      108.05, 21.55,
                      108.05, 18.8,
                      110.2, 13.2, 
                      110.2, -7.9,
                      94.3, -7.9,
                      94.8, 23.1,
                      97.4, 23.8),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - countries50[countries50$name == "China",]

sprp_barang <- raster::intersect(sprp, clip_poly)
sprp_barang@data$subgroup <- "barang"

# hainana
sprp_hainana <- sprp[sprp$island == "Hainan" & !is.na(sprp$island),]
sprp_hainana@data$subgroup <- "hainana"

# chinensis
clip_poly <- matrix(c(97.4, 23.8,
                      97.55, 28.5, 
                      97.55, 37.8,
                      123.1, 37.8,
                      123.1, 20.2,
                      97.4, 20.2, 
                      97.4, 23.8),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - sprp_barang - countries50[countries50$name == "Myanmar",]

sprp_chinensis <- raster::intersect(sprp, clip_poly)
sprp_chinensis@data$subgroup <- "chinensis"

# aurobrunnea x kutab x monticolus
clip_poly <- matrix(c(101, 38.2,
                      101, 21.2,
                      71.8, 21.2,
                      71.8, 38.2,
                      101, 38.2),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- gUnion(clip_poly, countries50[countries50$name == "Tajikistan",]) - sprp_chinensis - sprp_barang

poly_temp <- matrix(c(69.45, 39.7,
                      68.3, 38.1,
                      66.1, 40.1,
                      70.4, 41.5, 
                      71.6, 39.7,
                      69.45, 39.7),
                    ncol = 2, byrow = TRUE)
poly_temp <- Polygon(poly_temp)
poly_temp <- SpatialPolygons(list(Polygons(list(poly_temp), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - poly_temp

sprp_akm <- raster::intersect(sprp, clip_poly)
sprp_akm@data$subgroup <- "aurobrunnea x kutab x monticolus"

# seistanica
clip_poly <- matrix(c(59.6, 28.7,
                      56.6, 38.5,
                      60.7, 44.6, 
                      69.2, 44.6,
                      69.3, 40,
                      70, 38.9,
                      72.3, 38.5, 
                      76.5, 28.7,
                      59.6, 28.7),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - sprp_akm - countries50[countries50$name %in% c("Kyrgyzstan", "Iran"),]

sprp_seistanica <- raster::intersect(sprp, clip_poly)
sprp_seistanica@data$subgroup <- "seistanica"

# meridionalis
clip_poly <- matrix(c(44.1, 39.31,
                      44.05, 39.34, 
                      41.2, 39.3,
                      41.2, 43.5,
                      49.8, 43.5, 
                      51.7, 38.1,
                      60.2, 38.6, 
                      66.9, 32.4,
                      56.1, 26.3, 
                      40.5, 27.7,
                      46.3, 36.6,
                      44.1, 39.31),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - countries50[countries50$name %in% c("Turkey", "Turkmenistan", "Afghanistan"),]

sprp_meridionalis <- raster::intersect(sprp, clip_poly)
sprp_meridionalis@data$subgroup <- "meridionalis"

# lutra
sprp_try <- raster::bind(sprp_angustifrons, sprp_meridionalis, sprp_seistanica, sprp_nair, sprp_akm, sprp_barang, sprp_chinensis, sprp_hainana)
sprp_lutra <- sprp - sprp_try
sprp_lutra@data$subgroup <- "lutra"

sprp_try <- raster::bind(sprp_lutra, sprp_try)
tm_shape(sprp_try) +tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName == "Lutra lutra", "subgroup"] <- "lutra"

## Lycalopex culpaeus ####
# only one sample point :o
curr.species <- "Lycalopex culpaeus"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Lycalopex fulvipes ####
# only one sample
curr.species <- "Lycalopex fulvipes"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Lycalopex gymnocercus ####
# (wiki)
"Five subspecies are currently recognised, although the geographic range of each is unclear, and the type 
localities of three of them lie outside the present-day range of the species:[1][4]"
# Also has an apparently continuous range and subspecies range descriptions are not 

curr.species <- "Lycalopex gymnocercus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Lycaon pictus ####
# five subspecies in wiki, though:
"these subspecific designations are not universally accepted"
"extensive intermixing has occurred between East African and Southern African populations in the past"
"Southern African and northeastern African populations, with a transition zone"
"West African wild dog population may possess a unique haplotype, thus possibly constituting a truly distinct subspecies"
"The West African wild dog used to be widespread from western to central Africa, from Senegal to Nigeria. "
# Seems as though there has been lots of mixing between the subspecies, apart from manguensis
# I at least have pictus and lupinus based on gbif, probably also somalicus
# Too hard to know what's what and where there is or has been contact between subpopulations

curr.species <- "Lycaon pictus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("legend") + tm_shape(sp.gmpd.points) + tm_dots() 

## Lynx canadensis ####
# cat group says "Therefore we conclude that Lynx canadensis is a monotypic species"

## Lynx lynx ####
# On the basis of current evidence we propose the following six subspecies" (cat group)
# I only have ssp. lynx is described as:
"Distribution: Scandinavia, Finland, Baltic States, Belarus, European part of Russia E to the Yenissei River."
# Also see their map
curr.species <- "Lynx lynx"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# using Yenisei to try and replicate cat map
# I am assuming that the Yenisei acts as a biogeographic boundary between subspecies
# Although the range does extend South beyond the extent of the river

Yeni_Temp <- River_Data50[River_Data50$name == "Yenisey",]
Yeni_Temp <- disaggregate(Yeni_Temp)

Yeni_coords <- unlist(coordinates(Yeni_Temp), recursive = FALSE)[[1]]
Yeni_coords <- Yeni_coords[nrow(Yeni_coords):1,]

YA_lynx_coords <- rbind(Yeni_coords,
                         matrix(c(Yeni_coords[nrow(Yeni_coords), 1], (sprp@bbox[2,2] + 1),
                                  (sprp@bbox[1,2] + 1), (sprp@bbox[2,2] + 1),
                                  (sprp@bbox[1,2] + 1), (sprp@bbox[2,1] - 1),
                                  1.4, (sprp@bbox[2,1] - 1),
                                  1.4, 56,
                                  38, 47,
                                  80, 49,
                                  86, 47.5,
                                  100, 42,
                                  109, 42,
                                  96, 48,
                                  Yeni_coords[1,1], Yeni_coords[1,2]), 
                                ncol = 2, byrow = TRUE))

clip_poly <- Polygon(YA_lynx_coords)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

tm_shape(clip_poly) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(alpha = 0)

sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_lynx <- sprp - clip_poly
sprp_lynx@data$subgroup <- "lynx"

sprp_try <- raster::bind(sprp_otherssp, sprp_lynx)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "lynx"

## Lynx pardinus ####
# This is a monotypic species - cat group

## Lynx rufus ####
# cat group splits into 2 on either side of the great plains which acts as a barrier
# Nice to see it's even apparent in the locations of parasite samples
# They also group mexican populations into a rough subspecies
curr.species <- "Lynx rufus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# Going for quite a crude split trying to replicate cat map

buffer_temp <- countries50[countries50$name == "Mexico",]
buffer_temp <- terra::buffer(buffer_temp, 0.1)
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(buffer_temp) + tm_polygons(alpha = 0)

sprp_mexican <- raster::intersect(sprp, buffer_temp)
sprp_mexican@data$subgroup <- "mexican"

sprp_north <- sprp - buffer_temp

clip_poly <- matrix(c(-101, 29.5,
                      -102, 32,
                      -101.5, 35.5,
                      -100, 39,
                      -101.5, 45,
                      -106, 49,
                      -110, 51,
                      -58, 51,
                      -58, 23,
                      -101, 23,
                      -101, 29),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

clip_poly <- clip_poly - buffer_temp

sprp_rufus <- raster::intersect(sprp_north, clip_poly)
sprp_rufus@data$subgroup <- "rufus"

sprp_fasciatus <- sprp_north - clip_poly
sprp_fasciatus@data$subgroup <- "fasciatus"

sprp_try <- raster::bind(sprp_rufus, sprp_mexican, sprp_fasciatus)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Martes americana ####
# thirteen subspecies noted by msotw but no range descriptions, and I have points across the range
# "Martes americana, which is now considered to comprise two subspecies-groups (americana and caurina)" - IUCN
# "found to intergrade in Montana and British Columbia" - IUCN
# Also don't know of any island or otherwise isolated subspecies
curr.species <- "Martes americana"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Martes foina ####
"Some of the island populations are morphologically quite distinct, although the taxonomic significance 
of this is not yet clear"  #IUCN
# based on wiki descriptions of eleven subspecies, I have:
# foina, mediterranea
# but not: bosniaca, bunites, intermedia, kozlovi, milleri, nehringi, rosanowi, syriaca, toufoeus

# bosniaca, nehringi, rosanowi, and syriaca are contiguous with foina
# intermedia, bunites, kozlovi, milleri, toufoeus are geographically separate in IUCN polygons
# bunites and milleri are island populations

# bunites: only crete is included in polygon, not skopelos, erimomilos, karpathos, samothrace, seriphos, kythnos, or naxos

curr.species <- "Martes foina"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- c("Mongolia", "China", "Kyrgyzstan", "Tajikistan", "Uzbekistan", "Afghanistan", "Turkmenistan", "Pakistan")
buffer_temp <- countries50[countries50$name %in% buffer_temp,]
buffer_temp <- raster::aggregate(buffer_temp)
buffer_temp <- terra::buffer(buffer_temp, 2.8)
tm_shape(buffer_temp) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(alpha = 0)

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$subgroup <- "western"
sprp_try@data <- sprp_try@data %>%
  mutate(subgroup = case_when(island == "Crete" ~ "bunites",
                             island == "Rhodes" ~ "milleri",
                             TRUE ~ subgroup))

sprp_eastern <- raster::intersect(sprp, buffer_temp)
sprp_eastern <- raster::bind(sprp_eastern, sprp_try[4,])
sprp_eastern@data$subgroup <- "eastern"

sprp_try <- sprp_try - sprp_eastern
sprp_try <- raster::bind(sprp_try, sprp_eastern)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "western"

## Martes martes #####
# has eight subspecies listed by msotw, but no range descriptions available anywhere
curr.species <- "Martes martes"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Martes melampus ####
# I only have melampus
"The two confirmed subspecies of Japanese marten are:
M. m. melampus lives on several of the Japanese islands.
M. m. tsuensis is found on Tsushima Island.[3]" # wiki

curr.species <- "Martes melampus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("island") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & (is.na(island) | island == "Kamishima") ~ "tsuensis",
                              binomial == curr.species ~ "melampus",
                              TRUE ~ subgroup))

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup")

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "melampus"

## Martes pennanti ####
# "in general, the fisher is recognized to be a monotypic species with no extant subspecies.[11]" (wiki)

## Meles meles ####
# based on wiki description I have subspecies meles, marianensis, and possibly milleri, 
# not heptneri, though I can't remove the range of this one as its description overlaps with meles
# Need to remove the range of M. canescens by splitting at the Caucasus mountains as it's now a separate species (msotw)

curr.species <- "Meles meles"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- c("Georgia", "Azerbaijan", "Armenia", "Syria", 
                 "Israel", "Lebanon", "Palestine", "Jordan", 
                 "Iraq", "Iran", "Egypt", "Afghanistan",
                 "Turkmenistan", "Uzbekistan", "Tajikistan")
buffer_temp <- countries50[countries50$name %in% buffer_temp,]

buffer_temp <- raster::bind(buffer_temp, turkey_rangepol[c(1,2),])

buffer_temp <- raster::aggregate(buffer_temp)
buffer_temp <- terra::buffer(buffer_temp, 0.1)
buffer_temp <- buffer_temp - terra::buffer(turkey_rangepol[3,], 0.05)

sprp_try <- sprp - buffer_temp
tm_shape(sprp_try) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(alpha = 0)

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Melogale moschata/subaurantiaca ####
# Melogale subaurantiaca was previously considered subspecies of M. moschata
# now considered a full species occupying Taiwan
# All my samples appear to be M subaurantiaca so will adjust both polygons and gmpd

curr.species <- "Melogale moschata"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

GMPD_Data <- GMPD_Data %>% 
  mutate(HostCorrectedName = case_when(HostCorrectedName == curr.species ~ "Melogale subaurantiaca",
                                       TRUE ~ HostCorrectedName)) %>%
  mutate(HostReportedSubspecies = case_when(HostCorrectedName == "Melogale subaurantiaca" ~ NA_character_,
                                            TRUE ~ HostReportedSubspecies))

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(binomial = case_when(binomial == curr.species & island == "Taiwan" ~ "Melogale subaurantiaca",
                              TRUE ~ binomial))

curr.species <- "Melogale subaurantiaca"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Mephitis mephitis ####
# I have good sample coverage over most subspecies and they have a continuous range over North America
# Not aware of any island or otherwise isolated subspecies
curr.species <- "Mephitis mephitis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Mungos mungo ####
# only one sample location means this one would be dropped later anyway
curr.species <- "Mungos mungo"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Mustela erminea ####
# Following wiki description of 21 subspecies ranges, where available
# North America: kadiacensis, arctica, polaris
# range description of arctica does not include US or Eastern Canada, 
# but no subspecies attributable to this area and range is contiguous with arctica anyway

# Europe: hibernica, ricinae, stabilis, minima, erminea, aestiva
# IUCN polygon does not include Caucasus where teberdina is described
# same goes for ricina in Hebrides
# treating minima and erminea as part of aestiva as not sure where to split
# splitting at the Urals

# Asia: tobolica, ferghanae, mongolica, lymani, nippon, karaginensis, kaneii
# Uncertain how to split south of the Yenisei, but none of these are represented anyway
curr.species <- "Mustela erminea"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# North America
# kadaiacensis
sprp_kadiacensis <- sprp[sprp$island == "Kodiak" & !is.na(sprp$island),]
sprp_kadiacensis@data$subgroup <- "kadiacensis"

# polaris
sprp_polaris <- sprp[sprp$island == "Greenland" & !is.na(sprp$island),]
sprp_polaris@data$subgroup <- "polaris"

# arctica
clip_poly <- countries50[countries50$name %in% c("Canada", "United States"),]
clip_poly <- terra::buffer(clip_poly, 0.3)
sprp_arctica <- raster::intersect(sprp, clip_poly)
sprp_arctica <- sprp_arctica[!sprp_arctica$island %in% c("Kodiak", "Greenland"),]
sprp_arctica@data$subgroup <- "arctica"

# Europe
# hibernica
poly_temp <- countries50[countries50$name == "United Kingdom",]
poly_temp <- raster::disaggregate(poly_temp)
clip_poly <- countries50[countries50$name %in% c("Ireland", "Isle of Man"),]
clip_poly <- raster::bind(clip_poly, poly_temp[3,])
clip_poly <- terra::buffer(clip_poly, 0.15)

sprp_hibernica <- raster::intersect(sprp, clip_poly)
sprp_hibernica@data$subgroup <- "hibernica"

# stabilis
clip_poly <- poly_temp - clip_poly
sprp_stabilis <- raster::intersect(sprp, clip_poly)
sprp_stabilis@data$subgroup <- "stabilis"

# aestiva
clip_poly <- matrix(c(67.4, 68.8,
                      66.1, 68,
                      65.7, 67.2,
                      63.5, 66.5,
                      62.8, 65.8, 
                      61.1, 65,
                      59.9, 65.1,
                      58.9, 59.4,
                      59.7, 55.3, 
                      57.8, 54.5, 
                      57.1, 52.5,
                      57.2, 50.6,
                      51.3, 45.8,
                      -13.5, 39,
                      -13.5, 71.3, 
                      64.1, 71.3,
                      67.4, 68.8),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - terra::buffer(countries50[countries50$name %in% c("Kazakhstan", "United Kingdom", "Ireland"),], 0.15)

sprp_aestiva <- raster::intersect(sprp, clip_poly)
sprp_aestiva@data$subgroup <- "aestiva"

# remaining subspecies
sprp_try <- raster::bind(sprp_kadiacensis, sprp_polaris, sprp_arctica, sprp_stabilis, sprp_hibernica, sprp_aestiva)
sprp_otherssp <- sprp - sprp_try
sprp_otherssp@data$subgroup <- "tobolica x kaneii"
sprp_try <- raster::bind(sprp_try, sprp_otherssp)

tm_shape(sprp_try) +tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots("Prevalence")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Mustela lutreola ####
# seven subspecies, I only seem to have the french mink (biedermanni)
# uncertain if any subspecies are isolated
# but either way it's not possible to remove any effectively because of continuous range
curr.species <- "Mustela lutreola"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Mustela nivalis ####
# Basing division on geography, rather than morphological size based categories (from wiki)
# North American subspecies: allegheniensis, campestris, eskimo, rixosa
# European/Asian: nivalis, boccamela, caucasica, heptneri, namiyei, numidica, pallida, pygmaea, russelliana, vulgaris
# from range descriptions, island populations seem to be treated as part of mainland subspecies, e.g. pygmaea or numidica
curr.species <- "Mustela nivalis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- countries50[countries50$continent == "North America",]
clip_poly <- raster::aggregate(clip_poly)
clip_poly <- terra::buffer(clip_poly, 0.8)

sprp_americas <- raster::intersect(sprp, clip_poly)
sprp_americas@data$subgroup <- "american"

sprp_eurasian <- sprp - clip_poly
sprp_eurasian@data$subgroup <- "eurasian"

sprp_try <- raster::bind(sprp_americas, sprp_eurasian)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "eurasian"

## Mustela putorius ####
# msotw has seven subspecies, and based on wiki ranges I have:
# I defs have putorius and probably furo
# not aureola, mosquensis, or rothschildi but these are contiguous with putorius so I won't remove
# not anglia or caledoniae and I can remove as they are isolated

curr.species <- "Mustela putorius"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)

buffer_temp <- countries50[countries50$name == "United Kingdom",]
buffer_temp <- terra::buffer(buffer_temp, 0.3)

sprp_ukssp <- raster::intersect(sprp_try, buffer_temp)
sprp_ukssp@data$subgroup <- "anglia x caledoniae"

sprp_putorius <- sprp_try - buffer_temp
sprp_putorius@data$subgroup <- "putorius and european"

sprp_try <- raster::bind(sprp_putorius, sprp_ukssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "putorius and european"

## Nasua nasua ####
# only one sample location
curr.species <- "Nasua nasua"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Neovison vison ####
# 15 subspecies in msotw, with range descriptions from wiki
# I have multiple subspecies represented: vison, energumenos, evergladensis, and others 
# only one that would be reasonable to exclude is nesolestes but Admiralty Island is not included in the range anyway
curr.species <- "Neovison vison"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Nyctereutes procyonoides ####
# Going to follow current IUCN treatment of japanese populations as subspecies 
# I don't have: koreensis, orestes, procyonoides, ussuriensis
"There are six recognized subspecies: albus, koreensis, orestes, procyonoides, ussuriensis and viverrinus. 
It has been suggested recently that Japanese Raccoon Dogs should be classified as a distinct species Nyctereutes 
viverrinus with two subspecies N. v. viverrinus and N. v. albus (Sang-In et al. 2015). The classification is based 
on chromosomal, molecular and morphological differences between Japanese and mainland populations (see also Kauhala 
and Saeki 2004a)." #IUCN

"The raccoon dogs from Hokkaido are sometimes recognized as a different subspecies from the mainland tanuki as 
Nyctereutes procyonoides albus (Hornaday, 1904) (or N. viverrinus albus if recognized as a distinct species). 
This taxon is synonymized with N. p. viverrinus in Mammal Species of the World,[13][16] but comparative morphometric 
analysis supports recognizing the Hokkaido population as a distinct subspecific unit.[13]" #wiki

curr.species <- "Nyctereutes procyonoides"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("island") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & island == "Hokkaido" ~ "albus",
                              binomial == curr.species & str_detect(island, "Honshu|Kyushu|Sado|Shikoku") ~ "viverrinus",
                              TRUE ~ subgroup))

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Odocoileus hemionus ####
# 10 subspecies according to msotw, also recognised by IUCN
# Following wiki description of ranges: 
# I have californicus, columbianus, crooki, hemionus, inyoensis
# Don't have cerrosensis, fuliginatus, peninsulae, sheldoni, sitkensis
# Don't need to exclude cerrosensis or sheldoni as Cedros Island and Tiburon island not included in polygon anyway

# the IUCN range polygon differs from the ranges reported on wiki
# if wiki is correct, the remaining subspecies not represented in gmpd are contiguous with those that are 
# if the IUCN range is correct, peninsulae subspecies is isolated

# Going to leave it be as removal of peninsulae wouldn't affect range extent anyway

curr.species <- "Odocoileus hemionus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Odocoileus virginianus ####
# following wiki image of subspecies distributions
# Most of North American range is continuous and generally well sampled
# North american island subspecies: hiltonensis, mcilhennyi, rothschildi, taurinsulae, venatorius
# Some of these small islands are not represented in IUCN range polygons, but are inconsequential to analysis anyway

# Gap in wiki polygons is at panama border so I will split there
# https://commons.wikimedia.org/wiki/File:Odocoileus_virginianus_SA_map.svg

curr.species <- "Odocoileus virginianus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- countries50[countries50$continent == "North America",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with

sprp_north <- raster::intersect(sprp, buffer_temp)
sprp_north@data$subgroup <- "northern"
sprp_south <- sprp - buffer_temp
sprp_south@data$subgroup <- "southern"

sprp_try <- raster::bind(sprp_north, sprp_south)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "northern"

## Oreamnos americanus ####
# no subspecies on iucn, wiki, or msotw

## Otocolobus manul ####
# only one sample point :(
curr.species <- "Otocolobus manul"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]

GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Otocyon megalotis ####
# already got the subspecies <3
curr.species <- "Otocyon megalotis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Ourebia ourebi ####
# only extant isolated subspecies is already labelled, so referring to main population as nominate subspecies:
"Numerous subspecies of the Oribi have been described but most of these reflect individual variation and have little 
or no validity. Haggard's Oribi (O. o. haggardi) of eastern coastal Kenya and adjacent Somalia is a geographically 
isolated subspecies which is well differentiated in size and colour from other Oribi. Another distinctive subspecies 
from East Africa, the Kenya Oribi (O. o. keniae) from the lower slopes of Mount Kenya, is now apparently extinct."

curr.species <- "Ourebia ourebi"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & is.na(subgroup) ~ "oribi",
                              TRUE ~ subgroup)) 
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "oribi"

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup")

## Ovibos moschatus ####
# no subspecies on iucn, wiki, msotw

## Ovis ammon ####
# Sample locations don't correspond at all to Ovis ammon range
# see "Subspecies info" as well 
curr.species <- "Ovis ammon"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Ovis canadensis ####
# Not following the IUCN usage of 7 subspecies. following wiki instead as it's more recent
# I have: canadensis, sierrae, and nelsoni

"The 2016 genetics study suggested more modest divergence of this desert bighorn sheep into three lineages 
consistent with the earlier work of Cowan: Nelson's (O. c. nelsoni), Mexican (O. c. mexicana), and Peninsular 
(O. c. cremnobates). These three lineages occupy desert biomes that vary significantly in climate, suggesting 
exposure to different selection regimens."

# unable to split according to lineages or subspecies as their range descriptions are not specific enough

curr.species <- "Ovis canadensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Ovis dalli ####
# IUCN states there are two subspecies, but they are known to admix
# see Fannin sheep (O. d. fannini)

## Ozotoceros bezoarticus ####
# Not entirely sure which subspecies I have (only leucogaster, or bezoarticus as well) but can
"O. b. bezoarticus (Linnaeus, 1758), in the Cerrado ranging from eastern and central Brazil, south of the Amazon river 
between the plateau of Mato Grosso and the upper San Francisco river. While there is a population further north on the 
island of Marajo, at the mouth of the Amazon River (Rossetti and de Toledo 2006).

O. b. celer Cabrera, 1943, inhabiting the entire Argentinean pampas from the Atlantic coast to Sub Andean foothills and 
southward to the Rio Negro basin.

O. b. leucogaster (Goldfüss, 1817), located in the seasonally flooded grasslands of southwestern Brazil in southern 
Mato Grosso, southeastern Bolivia, Paraguay and in the Chaco savannas in northern Argentina (northern Santiago del 
Estero, Santa Fe, Formosa and Corrientes).

O. b. arerunguaensis González, Álvares-Valin and Maldonado, 2002. Type locality - Uruguay; northwestern Salto 
Department, Arerunguá, El Tapado, 31º41’51”S, 56º43’31”W.

O. b. uruguayensis González, Álvares-Valin and Maldonado, 2002. Type locality - Uruguay, grasslands of eastern Rocha 
Department, Sierra de Los Ajos, 33º50’01”S; 54º01’34”W."

curr.species <- "Ozotoceros bezoarticus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# celer
clip_points <- matrix(c(-62.8, -38.8,
                        -66.3, -34.8,
                        -57.2, -36.7),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- raster::disaggregate(sprp)
sprp_celer <- sprp_try[clip_points,]
sprp_celer@data$subgroup <- "celer"

# uruguayensis
clip_points <- matrix(c(-54, -33.5),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_uruguayensis <- sprp_try[clip_points,]
sprp_uruguayensis@data$subgroup <- "uruguayensis"

# arerunguaensis
clip_points <- matrix(c(-56.7, -31.3),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_arerunguaensis <- sprp_try[clip_points,]
sprp_arerunguaensis@data$subgroup <- "arerunguaensis"

# leucogaster
clip_points <- matrix(c(-60.8, -28.2,
                        -56.2, -28,
                        -51, -27.8,
                        -57.4, -22.2,
                        -58.5, -13.4,
                        -62, -16.1,
                        -60.9, -14.5,
                        -63.2, -13,
                        -68.3, -12.6,
                        -58.5, -17),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_leucogaster <- sprp_try[clip_points,]
sprp_leucogaster@data$subgroup <- "leucogaster"

# bezoarticus
clip_points <- matrix(c(-51.8, -22.3,
                        -45.9, -21,
                        -52.2, -16.4,
                        -45.5, -12.9,
                        -52, -12.1, 
                        -49.6, -9.5,
                        -50.8, -6.5),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_bezoarticus <- sprp_try[clip_points,]
sprp_bezoarticus@data$subgroup <- "bezoarticus"

sprp_try <- raster::bind(sprp_celer, sprp_uruguayensis, sprp_arerunguaensis, sprp_leucogaster, sprp_bezoarticus)
tm_shape(sprp_try) + tm_polygons("subgroup")

## Paguma larvata ####
# I only have samples in their introduced range
curr.species <- "Paguma larvata"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Panthera leo ####
# following the now outdated subspecies taxonomy of Asia (persica) and Africa (leo) as Asian population
# is so geographically isolated
curr.species <- "Panthera leo"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & is.na(dist_comm) ~ "leo",
                              binomial == curr.species ~ "persica",
                              TRUE ~ subgroup))

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "leo"

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Panthera onca ####
# Cat group describes it as a monotypic species, 
# but there are sufficiently isolated subpopulations that I'm going to split
# Can't split South American subgroups because of lack of specificity over which Amazon tributaries
# Also does not appear to be as isolated as Northern and Central American subgroups
"The status of the subspecies is unclear. Although eight subspecies have been recognized (Seymour 1989), 
morphological and genetic analyses do not support the existence of discrete subspecies (Larson 1997, Eizirik 
et al. 2001, Ruiz-Garcia et al. 2006). While not elevating the regional differences to the subspecies level, 
Eizirk et al. (2001) found evidence for four incompletely isolated phylogeographic groups: Mexico and Guatemala,
southern Central America, northern South America, and South America south of the Amazon River. Similarly, 
Ruiz-Garcia et al. (2006) found that the Andes Mountains incompletely isolates Jaguar populations in Colombia."

curr.species <- "Panthera onca"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mexico and Guatemala
clip_poly <- countries50[countries50$name %in% c("Mexico", "Guatemala", "Belize"), ]
clip_poly <- terra::buffer(clip_poly, 0.3)
sprp_northern <- raster::intersect(sprp, clip_poly)
sprp_northern@data$subgroup <- "northern"

# southern Central America
clip_poly <- matrix(c(-71.3, 9.3,
                      -72.1, 8.9,
                      -72.1, 8.1,
                      -72.7, 7.5,
                      -72.7, 6,
                      -78.5, 0,
                      -80, -3,
                      -88, -3,
                      -88, 17,
                      -71.3, 17,
                      -71.3, 9.3),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
sprp_central <- raster::intersect(sprp, clip_poly)
sprp_central@data$subgroup <- "central"

# Southern
sprp_try <- terra::buffer(raster::bind(sprp_central, sprp_northcentral), 0)
sprp_south <- sprp - sprp_try
sprp_south@data$subgroup <- "southern"

sprp_try <- raster::bind(sprp_south, sprp_central, sprp_northcentral)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "southern"

## Panthera pardus ####
# already has subspecies designations <3
curr.species <- "Panthera pardus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "pardus"

## Pecari tajacu ####
# splitting at Panama border
"DNA studies suggest that P. tajacu may consist of at least two major clades or lineages comprising specimens 
from North/Central and South America (Gongora et al. 2006, 2011) with structural chromosomal differences 
(Gongora et al. 2000, Adega et al. 2006). Additional studies are needed to clarify whether intra and/or 
inter-specific genetic and chromosomal variation is occurring within this species."

curr.species <- "Pecari tajacu"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- countries50[countries50$continent == "North America" & countries50$name != "Trinidad and Tobago",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with

sprp_north <- raster::intersect(sprp, buffer_temp)
sprp_north@data$subgroup <- "northern"
sprp_south <- sprp - buffer_temp
sprp_south@data$subgroup <- "southern"

sprp_try <- raster::bind(sprp_north, sprp_south)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "northern"

## Pelea capreolus ####
# continuous range with samples across it, also no subspecies described by iucn, wiki, msotw

### Phacochoerus aethiopicus ####
# range only includes extant subspecies, delamerei
"the Cape Warthog, Phacochoerus aethiopicus aethiopicus, endemic to South Africa, became extinct in the 1870s"
"the Somali Warthog, P. a. delamerei, occurs in Kenya and the Horn of Africa."

curr.species <- "Phacochoerus aethiopicus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp, bbox = bbox(sp.gmpd.points)) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# Some suspicious GMPD samples
# odd ones are probably reintroduced for hunting

View(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species,])

" See Grubb and d’Huart (2010) for a detailed historic overview of the classification of Phacochoerus."
# https://doi.org/10.2982/028.099.0204

## Philantomba monticola ####
# lots of subspecies, most of which have contiguous ranges
# I seem to have bicolor which has a range description as follows:
# "The range extends from Zanzibar to the KwaZulu Natal region in South Africa."
# Will therefore leave the range as is for now
curr.species <- "Philantomba monticola"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Procapra gutturosa ####
"Procapra gutturosa has not given rise to distinct geographic races. Specimens from the Mongolian Altai are
indistinguishable from those of eastern Mongolia and hence the subspecies P. g. altaica described from 
Bayan-Tsagan-Gobi is not accepted (Sokolov and Lushchekina 1997)." #IUCN

## Procyon lotor ####
# Many subspecies. Based on wiki I appear to have:
# lotor, elucus, fuscipes, hirtus, litoreus, marinus, megalodous, pacificus, psora, simus
# All of these are mainland subspecies except:
# litoreus: "Coastal strip and islands of Georgia."
# marinus: "Keys of the Ten Thousand Islands Group, and adjoining mainland of southwestern Florida from Cape Sable north through the Everglades to Lake Okeechobee."
## Both of these are hard to separate and not completely isolated from mainland subspecies, so will include with mainland

# Don't appear to have
# excelsus (contiguous with others), auspicatus (polygon not included), gloveralleni (extinct and not included),
# grinnelli (isolated), hernandezii (contiguous), incautus (polygon not included), inesperatus (not fully isolated), 
# insularis (isolated), maynardi (polygon not included), minor (polygon not included), pallidus (continguous), 
# pumilus (contiguous), vancouverensis (isolated but polygon not included)
## only separating fully isolated subspecies (grinnelli and insularis)

curr.species <- "Procyon lotor"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# grinnelli
clip_poly <- matrix(c(-112.1, 27.9,
                      -107.7, 23,
                      -110.6, 21.5,
                      -115.3, 27.2,
                      -112.1, 27.9),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_grinnelli <- raster::intersect(sprp, clip_poly)
sprp_grinnelli@data$subgroup <- "grinnelli"

# insularis
clip_poly <- matrix(c(-106.3, 22.5,
                      -105.9, 21.5,
                      -106.3, 20.5,
                      -107.3, 21.5,
                      -106.3, 22.5),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_insularis <- raster::intersect(sprp, clip_poly)
sprp_insularis@data$subgroup <- "insularis"

sprp_try <- raster::disaggregate(sprp)
sprp_mainland <- sprp_try - sprp_grinnelli - sprp_insularis
sprp_mainland@data$subgroup <- "mainland"

sprp_try <- raster::bind(sprp_mainland, sprp_grinnelli, sprp_insularis)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "mainland"

## Procyon pygmaeus ####
# monotypic

## Puma concolor ####
# cat group recognises two subspecies: 
# concolor "South America, possibly excluding W of Andes in north."
# couguar "North and Central America, possibly N South America W of Andes."
# Also going to treat Florida as its own population as it's so isolated

curr.species <- "Puma concolor"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- countries50[countries50$continent == "North America",]
buffer_temp <- terra::buffer(buffer_temp, 0.08) # smallest buffer I could get away with
florida_temp <- terra::buffer(states50[states50$name == "Florida",], 0.1)

tm_shape(sprp) + tm_polygons(col = "red") + 
  tm_shape(buffer_temp) + tm_polygons(alpha = 0) +
  tm_shape(florida_temp) + tm_polygons(alpha = 0)

sprp_couguar <- raster::intersect(sprp, buffer_temp)
sprp_couguar <- sprp_couguar - florida_temp
sprp_couguar@data$subgroup <- "couguar"

sprp_florida <- raster::intersect(sprp, florida_temp)
sprp_florida@data$subgroup <- "florida"

sprp_concolor <- sprp - buffer_temp - florida_temp
sprp_concolor@data$subgroup <- "concolor"

sprp_try <- raster::bind(sprp_couguar, sprp_florida, sprp_concolor)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Rangifer tarandus ####
# 14 subspecies listed by msotw
# North America:
# I may have: caboti, caribou, granti, groenlandicus, pearyi, terraenovae
# I don't have: osborni (contiguous with others)

# Russia & Finland:
# I don't have: buskensis, fennicus, pearsoni, phylarchus, sibiricus, valentinae

# Europe:
# Have: platyrhynchus in Svalbard, tarandus in Norway

# Based on geographic overlap between subspecies, I'm grouping them as follows:
# all north american and Greenland subspecies = american
# Russia and Finland 
# Norway = tarandus
# Svalbard = platyrhynchus

# Excluding outlying points in North Norway because of overlap with semi-domestic herds (see Sami wiki page)
# Although caribou in North America are hunted, they remain undomesticated

curr.species <- "Rangifer tarandus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

buffer_temp <- countries50[countries50$continent == "North America",]
buffer_temp <- terra::buffer(buffer_temp, 0.3) 
tm_shape(buffer_temp) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(col = "red")

sprp_americas <- raster::intersect(sprp, buffer_temp)
sprp_americas@data$subgroup <- "american"

norway_temp <- countries50[countries50$name == "Norway",]
norway_temp <- terra::buffer(norway_temp, 0.3) 
tm_shape(norway_temp) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons(col = "red")

sprp_tara_platy <- raster::intersect(sprp, norway_temp)
sprp_tara_platy@data$subgroup <- "platyrhynchus"

sprp_tara_platy <- raster::disaggregate(sprp_tara_platy)
sprp_tara_platy@data[2, "subgroup"] <- "tarandus"

sprp_otherssp <- sprp - buffer_temp - norway_temp
sprp_try <- raster::bind(sprp_americas, sprp_tara_platy, sprp_otherssp)

tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Raphicerus campestris ####
# already done <3
curr.species <- "Raphicerus campestris"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "campestris"

## Redunca arundinum ####
# iucn, wiki, msotw list no subspecies

## Redunca fulvorufula ####
# already done <£
curr.species <- "Redunca fulvorufula"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "fulvorufula"

## Rupicapra pyrenaica ####
# partially done, and wiki states:
"R. p. pyrenaica (Pyrenean chamois): France and Spain
R. p. parva (Cantabrian chamois): Spain
R. p. ornata (Abruzzo chamois): Central and southern Italy"
# this agrees with IUCN geographic range description 

curr.species <- "Rupicapra pyrenaica"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_ornata <- raster::intersect(sprp, countries50[countries50$name == "Italy",])
sprp_ornata@data$subgroup <- "ornata"

sprp_parva <- raster::intersect(sprp[is.na(sprp$subspecies),], countries50[countries50$name == "Spain",])
sprp_parva@data$subgroup <- "parva"

sprp_try <- raster::bind(sprp_ornata,
                         sprp_parva,
                         sprp[!is.na(sprp$subspecies),])

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
sprp_try@data <- sprp_try@data[,1:29]
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Rupicapra rupicapra ####
"The species R. rupicapra is categorized into seven subspecies:

R. r. asiatica (Anatolian chamois or Turkish chamois): Turkey
R. r. balcanica (Balkan chamois): Albania, Bosnia and Herzegovina, Bulgaria, Croatia, northern Greece 
(the Pindus Mountains),[4] North Macedonia, Serbia, Montenegro, and Slovenia (isolated populations)
R. r. carpatica (Carpathian chamois): Romania
R. r. cartusiana (Chartreuse chamois): France
R. r. caucasica (Caucasian chamois): Azerbaijan, Georgia, Russia
R. r. rupicapra (Alpine chamois): Austria, France, Germany, Italy, Switzerland, Slovenia,[5] Slovakia (Veľká Fatra, Slovak Paradise)
R. r. tatrica (Tatra chamois): Slovakia (Tatras and Low Tatras) and Poland (Tatras)"

# asiatica, caucasica, carpatica, and tatrica subspecies are isolated and not represented
# Don't have cartusiana based on range description here: http://www.wilddocu.de/chartreuse-chamois-rupicapra-rupicapra-cartusiana/
# balcanica is a bit hard to separate, but need to do it anyway

curr.species <- "Rupicapra rupicapra"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# tatrica
tatrica_temp <- terra::buffer(countries50[countries50$name == "Slovakia",], 0.1)
sprp_tatrica <- raster::intersect(sprp, tatrica_temp)
sprp_tatrica@data$subgroup <- "tatrica"

# carpatica
carpatica_temp <- countries50[countries50$name == "Romania",]
sprp_carpatica <- raster::intersect(sprp, carpatica_temp)
sprp_carpatica@data$subgroup <- "carpatica"

# asiatica
asiatica_temp <- terra::buffer(turkey_rangepol[2,], 1)
sprp_asiatica <- raster::intersect(sprp, asiatica_temp)
sprp_asiatica@data$subgroup <- "asiatica"

# caucasica
caucasica_temp <- countries50[countries50$name == "Georgia",]
caucasica_temp <- terra::buffer(caucasica_temp, 1.5)
caucasica_temp <- caucasica_temp - asiatica_temp
sprp_caucasica <- raster::intersect(sprp, caucasica_temp)
sprp_caucasica@data$subgroup <- "caucasica"

# balcanica
balcanica_temp <- countries50[countries50$name %in% c("Greece", "Bulgaria", "Macedonia", "Albania", 
                                                      "Montenegro", "Serbia", "Bosnia and Herz.", "Kosovo"),]
balcanica_temp <- terra::buffer(raster::aggregate(balcanica_temp), 0.25)

clip_poly <- matrix(c(17, 44,
                      14.9, 44,
                      14.9, 45.3,
                      15.4, 45.5,
                      14.9, 45.8,
                      16.5, 46.5,
                      17, 44),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

balcanica_temp <- gUnion(balcanica_temp, clip_poly)
sprp_balcanica <- raster::intersect(sprp, balcanica_temp) 
sprp_balcanica@data$subgroup <- "balcanica"

# rupicapra and cartusiana
clip_poly <- raster::bind(tatrica_temp, carpatica_temp, asiatica_temp, caucasica_temp, balcanica_temp)
sprp_rupicapra <- sprp - clip_poly
sprp_rupicapra@data$subgroup <- "rupicapra"
sprp_rupicapra <- raster::disaggregate(sprp_rupicapra)
sprp_rupicapra@data[15, "subgroup"] <- "cartusiana"
sprp_rupicapra <- raster::aggregate(sprp_rupicapra, by = names(sprp_rupicapra))

sprp_try <- raster::bind(sprp_tatrica, sprp_carpatica, sprp_asiatica, sprp_caucasica, sprp_balcanica, sprp_rupicapra)
sprp_try@data <- sprp_try@data[,1:29]
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Rusa unicolor ####
# only one sample
curr.species <- "Rusa unicolor"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Saiga tatarica ####
# only have one sample point
curr.species <- "Saiga tatarica"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Spilogale gracilis ####
# I only seem to have phenax, but apart from amphialus there is no apparent barrier to admixture
# amphialus not sampled and not recognised in range polyons so don't need to adjust
"S. g. amphialus Dickey, 1929 — Channel Islands spotted skunk (Channel Islands of California)
S. g. gracilis Merriam, 1890 — from south-eastern Washington to the extreme west of Oklahoma
S. g. latifrons Merriam, 1890 — southwestern British Columbia to western Oregon
S. g. leucoparia Merriam, 1890 — southern Arizona, New Mexico, and Texas, and northern Mexico
S. g. lucasana Merriam, 1890 — southern Baja California
S. g. martirensis Elliot, 1903 — northern and central Baja California
S. g. phenax Merriam, 1890 — California"

curr.species <- "Spilogale gracilis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Spilogale putorius ####
# only have one sample location
curr.species <- "Spilogale putorius"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Suricata suricatta ####
# only one sample
curr.species <- "Suricata suricatta"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Sus scrofa ####
# 16 subspecies, with four regional groupings:
# Indonesian: vittatus (not represented)
# Eastern: sibiricus, ussuricus, leucomystax, riukiuanus, taivanus, moupinensis (none represented)
# Indian: davidi, cristatus (none represented)
# Western: scrofa, meridionalis, algira, attila, lybicus, nigripes
# only scrofa represented, but not possible to distinguish range boundaries in Iran and Iraq
# also lots of overlap in range descriptions of attila and lybicus
# can at least split into regional, algira, meridionalis, and majori
curr.species <- "Sus scrofa"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# Indonesian
clip_poly <- countries50[countries50$name %in% c("Malaysia", "Indonesia"),]
clip_poly <- terra::buffer(clip_poly, 0.1)
sprp_vittatus <- raster::intersect(sprp, clip_poly)
sprp_vittatus@data$subgroup <- "vittatus"

# Eastern
clip_poly <- matrix(c(95.3, 48.6,
                      95.3, 55.3,
                      147, 55.3,
                      147, 7.1,
                      105.1, 7.1,
                      100.4, 13.3,
                      99.9, 20.9,
                      96.4, 24.3,
                      97.6, 28.5,
                      97.2, 30.1,
                      95.3, 48.6),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - countries50[countries50$name == "Myanmar",]

sprp_eastern <- raster::intersect(sprp, clip_poly)
sprp_eastern@data$subgroup <- "eastern"

# Indian
clip_poly <- matrix(c(64.6, 25.1,
                      69.8, 31.9,
                      73.9, 35.2,
                      98.3, 30.9,
                      104.9, 4.6,
                      64.6, 4.6,
                      64.6, 25.1),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly <- clip_poly - sprp_eastern - sprp_vittatus

sprp_indian <- raster::intersect(sprp, clip_poly)
sprp_indian@data$subgroup <- "indian"

# algira
clip_poly <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly <- terra::buffer(clip_poly, 0.1)
sprp_algira <- raster::intersect(sprp, clip_poly)
sprp_algira@data$subgroup <- "algira"

# meridionalis
clip_poly <- matrix(c(-8.9, 36,
                      -1.5, 38.8,
                      0.5, 36.4,
                      -8.9, 36),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_meridionalis <- raster::intersect(sprp, clip_poly)
sprp_meridionalis <- raster::bind(sprp_meridionalis, sprp[sprp$island %in% c("Sardinia", "Corsica"),])
sprp_meridionalis@data$subgroup <- "meridionalis"

# majori
clip_poly <- matrix(c(8.51, 44.3,
                      8.52, 44.9,
                      13.2, 44.9,
                      19.6, 39,
                      16.4, 37.1,
                      8.51, 44.3),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_majori <- raster::intersect(sprp, clip_poly)
sprp_majori@data$subgroup <- "majori"

# Western
sprp_try <- raster::bind(sprp_vittatus, sprp_eastern, sprp_indian, sprp_algira, sprp_majori, sprp_meridionalis)
sprp_western <- sprp - sprp_try
sprp_western@data$subgroup <- "western"

# finishing off
sprp_try <- raster::bind(sprp_try, sprp_western)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName ==curr.species, "subgroup"] <- "western"

## Sylvicapra grimmia ####
# can't really do much with the available information, and range is continuously connected
"Fourteen subspecies were recognized by Grubb and Groves (2001) and Wilson (2013). Distribution is continuous, 
there are many cases of intergradation and geographical boundaries between forms have not been delineated 
accurately"
curr.species <- "Sylvicapra grimmia"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Syncerus caffer ####
# already done <3
curr.species <- "Syncerus caffer"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Taxidea taxus ####
# wiki:
"Four subspecies have been recognized based on differences in skull size and pelage colour (Long 1972): 
Taxidea taxus berlandieri, in the southern United States; 
T. t. jacksoni, in the north-central United States and southern Ontario in Canada; 
T. t. taxus, in the Great Plains ecosystem from the United States into the prairie provinces of Canada; 
and T. t. jeffersonii, in the western United States and southern British Columbia."
# Likely to be taxus, but range is continuous

curr.species <- "Taxidea taxus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Tragelaphus angasii ####
# no subspecies described by IUCN, wiki, or msotw
curr.species <- "Tragelaphus angasii"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Tragelaphus eurycerus ####
# only one sample location
curr.species <- "Tragelaphus eurycerus"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Tragelaphus oryx ####
# "Three subspecies of Common Eland have been recognized, although their validity requires investigation" - IUCN
# oryx, livingstonii, and pattersonianus
# according to wiki I have all three and their ranges are contiguous
curr.species <- "Tragelaphus oryx"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Tragelaphus scriptus ####
# very complex taxonomic history (see iucn taxonomic notes and wiki discussion)
# As a result (and because range is largely continuous) it's not possible to divide into subspecies
# although this would be desirable considering there are two divergent lineages that may be distinct species
curr.species <- "Tragelaphus scriptus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Tragelaphus spekii ####
# unable to split based on basins because of lack of specificity in descriptions
"The species might even be monotypic,[6] however, based on different drainage systems, three distinct subspecies are currently recognised:[14][15]

T. s. spekii (Speke, 1863): Nile sitatunga or East African sitatunga. Found in the Nile watershed.
T. s. gratus (Sclater, 1880): Congo sitatunga or forest sitatunga. Found in western and central Africa.
T. s. selousi (W. Rothschild, 1898): Southern sitatunga or Zambezi sitatunga. Found in southern Africa."

curr.species <- "Tragelaphus spekii"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Tragelaphus strepsiceros ####
# seems as though I have strepsiceros, which is connected to chora by rift valley
# https://doi.org/10.1046/j.1365-294x.2001.01205.x
# Based on above paper's finding that their genotypes are quite divergent, even in geographically close samples
# I'm going to split the polygon in Kenya, though the exact line will be a bit arbitrary
"T. s. strepsiceros – southern parts of the range from southern Kenya to Namibia, Botswana, and South Africa
T. s. chora – northeastern Africa from northern Kenya through Ethiopia to eastern Sudan, Somalia, and Eritrea
T. s. cottoni – Chad and western Sudan"
curr.species <- "Tragelaphus strepsiceros"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- matrix(c(26, 8,
                      17, 8,
                      17, 15,
                      26, 15,
                      26, 8),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_cottoni <- raster::intersect(sprp, clip_poly)
sprp_cottoni@data$subgroup <- "cottoni"

clip_poly <- matrix(c(37, -1.2,
                      35, -0.9,
                      32, -0.9,
                      32, 21,
                      48, 21,
                      48, -1.2,
                      37, -1.2),
                    ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_chora <- raster::intersect(sprp, clip_poly)
sprp_chora@data$subgroup <- "chora"

sprp_strepsiceros <- sprp - clip_poly - sprp_cottoni
sprp_strepsiceros@data$subgroup <- "strepsiceros"

sprp_try <- raster::bind(sprp_cottoni, sprp_chora, sprp_strepsiceros)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

GMPD_Data[GMPD_Data$HostCorrectedName ==curr.species, "subgroup"] <- "strepsiceros"

## Urocyon cinereoargenteus ####
# Many many subspecies, but only venezuelae (Colombia and Venezuela) is geographically isolated
curr.species <- "Urocyon cinereoargenteus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)
sprp_try@data[2, "subgroup"] <- "venezuelae"
sprp_try@data[c(1,3), "subgroup"] <- "northern"

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Urocyon littoralis ####
# Only have one sample location for each island so it's impossible to find any useful info
"Six distinct subspecies are recognized, one on each of the islands where they occur:

San Miguel Island Fox (Urocyon littoralis littoralis (Baird, 1858)), San Miguel Island,
Santa Rosa Island Fox (U. l. santarosae Grinnell & Linsdale, 1930), Santa Rosa Island,
Santa Cruz Island Fox (U. l. santacruzae Merriam, 1903), Santa Cruz Island,
Santa Catalina Island Fox (U. l. catalinae Merriam, 1903), Santa Caralina Island,
San Nicolas Island Fox (U. l. dickeyi Grinnell & Linsdale, 1930), San Nicolas Island, and
San Clemente Island Fox (U. l. clementae Merriam, 1903), San Clemente Island."

# I have Santa Catalina (catalinae), San Nicolas (dickeyi), and San Clemente (clementae)
curr.species <- "Urocyon littoralis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Ursus americanus ####
# subspecific taxonomy is too complex and mainland populations are too admixed to separate reasonably
# coupled with the fact that I have samples across the full range
# can therefore only really exclude island subspecies populations
# I don't have Haida Gwaii (carlottae), Dall Island (pugnax), Vancouver Island (vancouveri), Kenai (perniger)
# afaik these are the only geographically isolated subpops
# https://doi.org/10.1093/molbev/msv114
curr.species <- "Ursus americanus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# Dall island and Haida Gwaii
haida_temp <- matrix(c(-134, 54.5,
                      -132.7, 55,
                      -132, 54.4,
                      -131, 54.2,
                      -130.5, 51.5,
                      -133.5, 53,
                      -134, 54.5),
                    ncol = 2, byrow = TRUE)
haida_temp <- Polygon(haida_temp)
haida_temp <- SpatialPolygons(list(Polygons(list(haida_temp), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- sprp - haida_temp
haida_temp <- raster::intersect(sprp, haida_temp)

# Kenai peninsula
kenai_temp <- matrix(c(-150.5, 61.1,
                      -149.36, 60.91,
                      -149.156, 60.91,
                      -149.01, 60.83,
                      -148.7, 60.78,
                      -148.57, 60.81,
                      -148.46, 60.81,
                      -148.23, 60.77,
                      -147.4, 60.4,
                      -148, 59.9,
                      -151.5, 58.9,
                      -152.5, 59.4,
                      -151.6, 60.7, 
                      -150.5, 61.1),
                    ncol = 2, byrow = TRUE)
kenai_temp <- Polygon(kenai_temp)
kenai_temp <- SpatialPolygons(list(Polygons(list(kenai_temp), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- sprp_try - kenai_temp
kenai_temp <- raster::intersect(sprp, kenai_temp)

# Vancouver Island
clip_points <- matrix(c(-125.8, 49.8,
                        -126.7, 49.7,
                        -126.1, 49.3,
                        -125.8, 49.2,
                        -123.5, 48.8,
                        -123.4, 48.9,
                        -126.9, 50.65),
                      ncol = 2, byrow = TRUE)
clip_points <- SpatialPoints(clip_points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- raster::disaggregate(sprp_try)
vancouver_temp <- raster::intersect(sprp_try, clip_points)
sprp_try <- sprp_try - vancouver_temp
sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
sprp_try@data$subgroup <- "mainland"

sprp_try <- raster::bind(sprp_try, haida_temp, kenai_temp, vancouver_temp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)
GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, "subgroup"] <- "mainland"

### Ursus arctos ####
# Brown bear taxonomy and subspecies classification has been described as "formidable and confusing,"
# North America:
"DNA analysis shows that, apart from recent human-caused population fragmentation,[39] brown bears 
in North America are generally part of a single interconnected population system, with the exception 
of the population (or subspecies) in the Kodiak Archipelago, which has probably been isolated since 
the end of the last Ice Age.[40][41] "

curr.species <- "Ursus arctos"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# North America
# middendorffi
sprp_middendorffi <- sprp[sprp$island == "Kodiak" & !is.na(sprp$island),]
sprp_middendorffi@data$subgroup <- "middendorffi"

# ungavaensis 
clip_poly <- states50[states50$name %in% c("Québec", "Newfoundland and Labrador"),]
clip_poly <- terra::buffer(clip_poly)
sprp_ungavaensis <- raster::intersect(sprp, clip_poly)
sprp_ungavaensis@data$subgroup <- "ungavaensis"

# sitkensis
sprp_sitkensis <- sprp[sprp$island %in% c("Admiralty", "Chichago & Baranof"),]
sprp_sitkensis@data$subgroup <- "sitkensis"

# horribilis
clip_poly <- countries50[countries50$continent == "North America", ]
clip_poly <- terra::buffer(clip_poly, 0.6)
sprp_try <- raster::intersect(sprp, clip_poly)
sprp_horribilis <- sprp_try - sprp_sitkensis - sprp_ungavaensis - sprp_middendorffi
sprp_horribilis@data$subgroup <- "horribilis"

sprp_america <- raster::bind(sprp_middendorffi, sprp_ungavaensis, sprp_sitkensis, sprp_horribilis)
tm_shape(sprp_america) + tm_polygons("subgroup")

# Eurasia and North Africa
# crowtheri
clip_poly <- countries50[countries50$name %in% c("Morocco", "Algeria", "Tunisia"),]
clip_poly <- terra::buffer(clip_poly, 0.1)
sprp_crowtheri <- raster::intersect(sprp, clip_poly)
sprp_crowtheri@data$subgroup <- "crowtheri"

# pyrenaicus
clip_poly <- countries50[countries50$name == "Spain",]
clip_poly <- terra::buffer(clip_poly, 0.8) - sprp_crowtheri
sprp_pyrenaicus <- raster::intersect(sprp, clip_poly)
sprp_pyrenaicus@data$subgroup <- "pyrenaicus"

# marsicanus
clip_poly <- matrix(c(13, 44,
                      16, 42,
                      13, 40, 
                      9, 43,
                      13, 44),
                     ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
sprp_marsicanus <- raster::intersect(sprp, clip_poly)
sprp_marsicanus@data$subgroup <- "marsicanus"

# syriacus
clip_poly <-  matrix(c(23.9, 39, 
                       35, 43,
                       38, 45.5,
                       62, 43.6,
                       64, 37,
                       77, 27,
                       25, 27, 
                       23.9, 39),
                     ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- raster::disaggregate(sprp)
sprp_syriacus <- sprp_try[clip_poly,]
sprp_syriacus@data$subgroup <- "syriacus"

# pruinosus, isabellinus, gobiensis
clip_poly <-  matrix(c(67, 48,
                       91, 45, 
                       106, 46,
                       106, 27,
                       77, 27,
                       75.5, 32,
                       74.8, 33,
                       74, 34.8,
                       72, 35.5,
                       70, 37,
                       65, 37),
                     ncol = 2, byrow = TRUE)
clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
tm_shape(clip_poly) + tm_polygons(alpha = 0) +tm_shape(sprp) + tm_polygons(alpha = 0)

sprp_pruinosus <- sprp_try[clip_poly,]
sprp_pruinosus@data$subgroup <- "pruinosus x isabellinus x gobiensis"

#

tm_shape(sprp_pruinosus) + tm_polygons()

## Ursus maritimus ####
# no subpsecies on IUCN or msotw, and wiki states: 
"When the polar bear was originally documented, two subspecies were identified [...] This 
distinction has since been invalidated."
curr.species <- "Ursus maritimus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

## Vulpes corsac ####
# only one sample
curr.species <- "Vulpes corsac"
unique(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]
GMPD_Data <- GMPD_Data[GMPD_Data$HostCorrectedName != curr.species,]

## Vulpes lagopus ####
"Besides the nominate subspecies, the common Arctic fox, V. l. lagopus, four other subspecies of this fox have been described:

Bering Islands Arctic fox, V. l. beringensis
Greenland Arctic fox, V. l. foragoapusis
Iceland Arctic fox, V. l. fuliginosus
Pribilof Islands Arctic fox, V. l. pribilofensis" # (wiki)

# I don't have beringensis or pribilofensis
# I do have foragoapusis, fuliginosus, and lagopus

curr.species <- "Vulpes lagopus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

# foragoapusis (Iceland)
sprp_foragoapusis <- sprp[!is.na(sprp$island) & sprp@data$island == "Iceland",]
sprp_foragoapusis@data$subgroup <- "foragoapusis"

# fuliginosus (Greenland)
sprp_fuliginosus <- sprp[!is.na(sprp$island) & sprp$island == "Greenland",]
sprp_fuliginosus@data$subgroup <- "fuliginosus"

# lagopus
bering_temp <- c("St. Matthew", "Hall", "St. Lawrence")
pribilof_temp <- c("St. Paul", "St. George")
sprp_lagopus <- sprp[!sprp$island %in% c(bering_temp, pribilof_temp, "Greenland", "Iceland"),]
sprp_lagopus@data$subgroup <- "lagopus"

sprp_lagopus <- raster::disaggregate(sprp_lagopus)
sprp_lagopus@data[26, "subgroup"] <- "fuliginosus"
sprp_lagopus <- raster::aggregate(sprp_lagopus, by = names(sprp_lagopus))

# beringensis
sprp_beringensis <- sprp[sprp$island %in% bering_temp, ]
sprp_beringensis@data$subgroup <- "beringensis"

# pribilofensis
sprp_pribilofensis <- sprp[sprp$island %in% pribilof_temp, ]
sprp_pribilofensis@data$subgroup <- "pribilofensis"

sprp_try <- raster::bind(sprp_lagopus, sprp_foragoapusis, sprp_fuliginosus, sprp_beringensis, sprp_pribilofensis)
tm_shape(sprp_try) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots()

IUCN_Native_Data@data <- IUCN_Native_Data@data %>% 
  mutate(subgroup = case_when(binomial == curr.species & is.na(subgroup) ~ "macrotis",
                              TRUE ~ subgroup))

## Vulpes macrotis ####
"most available data suggest that kit foxes in the San Joaquin Valley of California are likely 
to warrant a subspecific designation, V. m. mutica, due to geographical isolation, and that any 
other kit foxes may be included in a second subspecies, V. m. macrotis"
curr.species <- "Vulpes macrotis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>% 
  mutate(subgroup = case_when(binomial == curr.species & is.na(subgroup) ~ "macrotis",
                              TRUE ~ subgroup))

## Vulpes velox ####
# no subspecies in IUCN, wiki, or msotw
curr.species <- "Vulpes velox"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup") + tm_shape(sp.gmpd.points) + tm_dots() 

### Vulpes vulpes ####
# 45 subspecies with patchy range descriptions which do not cover full extent of iucn range polygon
# In Europe I have: crucifera, silacea, vulpes
# In Africa: barbara
# In Eurasia, moving West to East: 
# palaestina/arabica, kurdistanica, flavescens, unknown, japonica, schrencki

curr.species <- "Vulpes vulpes"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

# barbarica x atlantica
# silacea
# ichnusae
# indutus
# niloticus
# crucifera
# arabica?
# vulpes?
# japonica
# schrencki
# splendidissima






# tidying up
rm(list = ls(pattern = "^sprp"))
rm(list = ls(pattern = "_temp$"))
rm(list = ls(pattern = "_Temp$"))







## discard pile ####
curr.species <- "Ovibos moschatus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sprp <- IUCN_Mammals[IUCN_Mammals$binomial == curr.species, ]

sp.gmpd.points <- GMPD_Spatial[GMPD_Spatial$HostCorrectedName == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostReportedName == "Viverra civetta", c("Longitude", "Latitude")])

tm_shape(World, bbox = sprp@bbox) + tm_polygons() + tm_shape(sprp) + tm_polygons("legend") + tm_shape(sp.gmpd.points) + tm_dots()



sp.bbox <- sprp@bbox
sprp$id <- rownames(sprp@data)
sprp$poly <- 1:36
sprp_fort <- merge(fortify(sprp), sprp@data, bi = "id")

# upper
coords = matrix(c(33.0, 0,
                  -10, 0,
                  -10, 40,
                  60.0, 40,
                  60.0, 3.7,
                  45.0, 3.7,
                  35.5, 5,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

# lower?
coords = matrix(c(33.0, 0,
                  35.5, 5,
                  45.0, 3.7,
                  60.0, 3.7,
                  60.0, -30,
                  -10, -30,
                  -10, 0,
                  33.0, 0), 
                ncol = 2, byrow = TRUE)

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = c(30,50), ylim = c(-5,15), expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sprp_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = as.data.frame(coords), aes(x = V1, y = V2), fill = NA, colour = "black")

clip_poly <- Polygon(coords)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
clip_poly_try <- gUnion(clip_poly, sprp[sprp$poly == 13, ])

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = sp.bbox[1,], ylim = sp.bbox[2,], expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sprp_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = fortify(clip_poly_try), aes(x = long, y = lat), fill = NA, colour = "black")

sprp.test <- st_as_sf(sprp)
clip_poly.test <- st_as_sf(clip_poly_try)
sprp_try <- st_difference(sprp.test, clip_poly.test)

ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(colour = "grey65", size = 0.15) + theme_bw() +
  coord_sf(xlim = sp.bbox[1,], ylim = sp.bbox[2,], expand = TRUE) +
  xlab("Longitude") + ylab("Latitude") +
  
  geom_polypath(data = sprp_fort, 
                aes(x = long, y = lat, group = group, colour = poly, fill = poly),
                size = 0.05) +
  geom_polypath(data = (sprp_try), aes(x = long, y = lat), fill = "black", colour = "black")



##
#library(mapview)
#geojson.io
sf::sf_use_s2(FALSE)
library(tmap)
tmap_mode("plot")
tm_shape(sprp_try) + tm_polygons("binomial")