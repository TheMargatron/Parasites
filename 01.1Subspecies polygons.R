# Figuring things out ####

# going by subgroup rather than subspecies because of geographic groupings
IUCN_Native_Data@data$subgroup <- IUCN_Native_Data@data$subspecies

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(str_detect(subgroup, " ssp. ") ~ str_split(subgroup, "ssp. ", simplify = TRUE)[,2],
                              TRUE ~ subgroup))


## Acinonyx jubatus ####
# rationale: There are multiple subspecies that are quite geographically separated, 
# but GMPD data only occupies the range of A. j. jubatus
# therefore want to remove all polygons representing other subspecies. 
# I've removed polygons based on wiki distribution map
# https://en.wikipedia.org/wiki/Cheetah#/media/File:Acinonyx_jubatus_subspecies_range_IUCN_2015.png

curr.species <- "Acinonyx jubatus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup")

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

tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Aepyceros melampus ####
# rationale: Two subspecies, only the common impala (subsp. melampus) is well represented in GMPD based on wiki and iucn maps
# Also black-faced impala seems geographically distinct 
# https://en.wikipedia.org/wiki/Impala#/media/File:Aepyceros_melampus.svg
# IUCN: "In Namibia, the Black-faced Impala is naturally confined to the Kaokoland in the north-west, and neighbouring south-western Angola"

curr.species <- "Aepyceros melampus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sprp) + tm_polygons("subgroup")  

## Alcelaphus buselaphus ####
# rationale: 8 subgroup, GMPD data appears to represent major and cokii, and they are geographically distinct
# https://en.wikipedia.org/wiki/Hartebeest#/media/File:Alcelaphus_recent.png
# There is only one sample location for each subspecies so I'll have to drop this one :(

curr.species <- "Alcelaphus buselaphus"
unique(GMPD_Data[which(GMPD_Data$HostCorrectedName == curr.species),c("Latitude", "Longitude")])

IUCN_Native_Data <- IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,]

## Alces alces ####
# rationale: There are multiple subspecies, GMPD represents shirasi, gigas, andersoni, and americana in North America
# and alces in Europe/Western Russia, but not buturlini, cameloides, or pfizenmayeri (east of Yenisei river)
# Splitting European polygon around Yenisei should do it. 
# Yenisei data source: https://doi.org/10.1016/j.dib.2018.09.016

curr.species <- "Alces alces"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]

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

YA_coords <- rbind(YA_coords,
                   matrix(c(YA_coords[nrow(YA_coords), 1], (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), (sprp@bbox[2,2] + 1),
                            (sprp@bbox[1,2] + 1), (sprp@bbox[2,1] - 1),
                            YA_coords[1,1], (sprp@bbox[2,1] - 1),
                            YA_coords[1,1], YA_coords[1,2]), 
                          ncol = 2, byrow = TRUE))

clip_poly <- Polygon(YA_coords)
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
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Antidorcas marsupilis ####
# rationale: There are three recognised subspecies (wiki) but their individual distributions are not well defined
# They are also not described by the IUCN
# I only seem to have A. m. marsupialis in GMPD data whose range lies south of Orange river (wiki) 
# In the Eastern part of this range it appears to be restricted by the Vaal rather than the Orange as wiki describes it as extending up to Kimberley which is North of Orange

curr.species <- "Antidorcas marsupialis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
OV_Temp <- River_Data50[(River_Data50$name == "Orange"| River_Data50$name == "Vaal"),]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

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

## Antilocapra americana ####
# rationale: three subspecies, but not all represented by GMPD
# "Mitochondrial DNA analyses since the early 1990s support the idea of clines within a wide-ranging species rather than separate subspecies (O’Gara and Yoakum 2004)."

curr.species <- "Antilocapra americana"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])

tm_shape(sprp) + tm_polygons("legend") +  tm_shape(sp.gmpd.points) + tm_dots()

sprp_try <- raster::disaggregate(sprp)
tm_shape(sprp_try[17,]) + tm_polygons()
sprp_try@data$subgroup <- "americana x sonoriensis"
sprp_try@data[17,"subgroup"] <- "peninsularis"
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Axis axis ####
# rationale: No description of subspecies in either IUCN or main wiki page, 
# but there's a page for sri lankan subspecies: https://en.wikipedia.org/wiki/Sri_Lankan_axis_deer
# and it's geographically isolated

curr.species <- "Axis axis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
tm_shape(sprp) + tm_polygons("island")

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & island == "Sri Lanka" ~ "ceylonensis",
                              binomial == curr.species ~ "axis",
                              TRUE ~ subgroup))

tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species,]) + tm_polygons("subgroup")

## Bison bison ####
# rationale: "There are two recognized subspecies in North America: Bison bison bison and B. b. athabascae."

curr.species <- "Bison bison"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

# not sure how to split these for now because gbif data does not have clear delineation
# and subspecies data seems to be the other way round in gbif?
# probs just won't split them as they are not too geographically distinct anyway

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species ~ "bison x athabascae",
                              TRUE ~ subgroup))

## Bison bonasus ####
# rationale: "Three subspecies of European bison existed in the recent past, but only one, 
# the nominate subspecies (B. b. bonasus), survives today" (wiki)
# two subspecies are widely recognized as the Lowland Bison (Bison bonasus bonasus) and 
# the Caucasian Bison (Bison bonasus caucasicus) (Kowalczyk and Plumb 2020)" (IUCN)

# "European bison herds, scattered across Central and Eastern Europe, represent two genetic lines"
# the lowland line (Poland, Belarus, and Lithuania) and the lowland-Caucasian line (southern Poland, Russia, Ukraine and Slovakia).
# https://animaldiversity.org/accounts/Bison_bonasus/

# animaldiversity doesn't describe them as subspecies, only genetic lines
# The distribution of the lineages is artificial because it's determined by reintroduction programmes, rather than natural dispersal
# That means that the different habitat patches they occupy will not relate to their lineage, just to their overall species distribution

curr.species <- "Bison bonasus"
IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species ~ "bonasus x caucasicus",
                              TRUE ~ subgroup))

## Blastocerus dichotomus ####
# rationale: No taxonomic notes on IUCN, no description of subspecies on wiki or animaldiversity.org

## Canis adustus ####
# rationale: There are seven recognized subspecies of the side-striped jackal:[2]

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

# I think the species where I can't distinguish between subspecies will still be informative for latitudinal gradient,
# but less informative for within-range analysis.
# Might be interesting to compare result from uncleaned polygons with cleaned. Split the streams again after adding this step?

curr.species <- "Canis adustus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

## Canis aureus ####
"Recent studies based on mtDNA and morphology have shown that 'Golden Jackals' in Africa are larger in size than those from 
Eurasia and are actually more closely related to the Grey Wolf Canis lupus. African animals hence represent a previously 
overlooked distinct species, the African Wolf, Canis lupaster (see Rueness et al. 2011, Gaubert et al. 2012, Koepfli et al. 
2015, Viranta et al. 2017). However, the putative presence of Golden Jackal in the Sinai Peninsula of Egypt remains unclear 
(see Gaubert et al. 2012, Viranta et al. 2017)."

# rationale: wiki s 7 subspecies
# aureus: Middle East, Iran, Turkmenistan, Afghanistan, Pakistan and Western India
# cruesemanni: Thailand
# ecsedensis: Pannonian Basin, Central Europe
# indicus: India, Nepal, Bangladesh, Bhutan
# moreotica: Southeastern Europe, Moldova, Asia Minor and the Caucasus
# naria: Coastal South West India, Sri Lanka
# syriacus: Israel, Syria,[38] Lebanon,[62] and Jordan

# Have aureus in Iran and probably moreotica in Greece
# subdividing according to : https://doi.org/10.1093/mspecies/sey002
# Following IUCN range description I'm not including the algirensis subspecies
curr.species <- "Canis aureus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

# moreotica
sprp_try <- raster::disaggregate(sprp)
sprp_moreotica <- sprp_try[2,]

# removing ecsedensis (Pannonian basin/Hungary)
sprp_moreotica <- sprp_moreotica - countries50[countries50$name == "Hungary",]
sprp_moreotica <- raster::disaggregate(sprp_moreotica)
sprp_moreotica <- sprp_moreotica[1,]

# adding Turkey
turkey_rangepol <- countries50[countries50$name == "Turkey",]
turkey_rangepol <- raster::intersect(sprp_try[1,], turkey_rangepol)
turkey_rangepol@data <- turkey_rangepol@data[, names(sprp_try)]

sprp_moreotica <- raster::bind(sprp_moreotica, 
                                           turkey_rangepol)

# adding range around azov sea, split at Georgia boundary
azov_rangepol <- Polygon(matrix(c(39.75, 43,
                                    34, 43,
                                    34, 49,
                                    43.75, 49,
                                    39.75, 43), 
                                  ncol = 2, byrow = TRUE))
azov_rangepol <- SpatialPolygons(list(Polygons(list(azov_rangepol), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
sprp_moreotica <- raster::bind(sprp_moreotica, 
                                           raster::intersect(azov_rangepol, sprp_try[1,]))
sprp_moreotica@data$subgroup <- "moreotica"
tm_shape(sprp_moreotica) + tm_polygons("subgroup")

# aureus
clip_poly <- raster::bind(turkey_rangepol, 
                       countries50[countries50$name == "Jordan",],
                       countries50[countries50$name == "Syria",],
                       countries50[countries50$name == "India",],
                       azov_rangepol)
sprp_aureus <- sprp_try[1,] - clip_poly
sprp_aureus <- raster::disaggregate(sprp_aureus)
sprp_aureus <- sprp_aureus[1,]
sprp_aureus@data$subgroup <- "aureus"
tm_shape(sprp_aureus) + tm_polygons("subgroup")

# adding it all up
sprp_try <- raster::bind(sprp_moreotica, sprp_aureus)
sprp_otherssp <- sprp - sprp_try
sprp_try <- raster::bind(sprp_try, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

# Getting rid of crumbs
sprp_try <- raster::disaggregate(sprp_try)
sprp_try@data$SHAPE_Area <- raster::area(sprp_try)
sprp_try <- sprp_try[sprp_try$SHAPE_Area > 100000000,]
sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))

tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Canis latrans ####
# rationale: although subspecies ranges are well described in wiki, there's overlap between them 
# and I don't know exactly where to draw lines. 
# Also range is pretty contiguous (apart from tiburon island). The only subspecies I'm lacking are central american ones

### Canis lupus ####
# for Canis lupus rufus / Canis rufus:
"See Chambers et al. (2012) for a brief review of recent literature concerning the status of this species, 
which they considered a full species, as does this assessment."
# rationale: absolutely loads of subspecies: https://en.wikipedia.org/wiki/Subspecies_of_Canis_lupus
# based on: https://upload.wikimedia.org/wikipedia/commons/1/16/Present_distribution_of_gray_wolf_%28canis_lupus%29_subspecies.png
# Not got: baileyi, arctos, albus, arabs, nubilis
# Got: occidentalis, lycaon, signatus, lupus, italicus, pallipes

# use IUCN taxonomic notes


# removing arctos
sprp_try <- sprp[-grep("Greenland|Ellesmere|Banks|Melville", sprp$island),]

## Canis mesomelas ####
# rationale: two geographically distinct subspecies
curr.species <- "Canis mesomelas"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

# For when I'm doing gmpd subspecies
# sp.over <- over(sp.gmpd.points, gBuffer(sprp, byid = TRUE))
# sp.gmpd$subspecies <- sp.over$subspecies
# sp.gmpd$HostCorrectedName <- paste(sp.gmpd$HostCorrectedName, sp.over$subspecies, sep = " ")

## Canis simensis ####
# rationale: Two geographically distinct subspecies and I only seem to have citernii
# https://en.wikipedia.org/wiki/Ethiopian_wolf#/media/File:Canis_simensis_subspecies_range.png
curr.species <- "Canis simensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

clip_poly = matrix(c(37, 9,
                  40, 9,
                  40, 14,
                  37, 14,
                  37, 9), 
                ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sprp_try <- sprp - clip_poly
tm_shape(clip_poly) + tm_polygons() + tm_shape(sprp_try) + tm_polygons()
sprp_try@data$subgroup <- "citernii"

sprp_otherssp <- raster::intersect(sprp, clip_poly)
sprp_otherssp@data$subgroup <- "simensis"

sprp_try <- raster::bind(sprp_try, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Capra ibex ####
# rationale: no reported subspecies

## Capra pyrenaica ####
# rationale: four subspecies, two extinct. I probably only have hispanica
# subdivided based on IUCN description of ranges
curr.species <- "Capra pyrenaica"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$temp <- letters[1:nrow(sprp_try@data)]
tm_shape(sprp_try) + tm_polygons("temp")

split_hispanica <- c("a", "b", "c", "d","e", "f", "g", "i", "m", "n")
sprp_try@data <- sprp_try@data %>%
  mutate(subgroup = case_when(temp %in% split_hispanica ~ "hispanica",
                              TRUE ~ "victoriae")) %>% 
  dplyr::select(-temp)

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Capreolus capreolus ####
# rationale: The IUCN notes the following five confirmed subspecies: italicus, garganta, capreolus, caucasicus, and coxi
# Definitely have italicus and garganta
# uncertain about capreolus but I expect so
# Don't appear to have coxi and caucasicus

curr.species <- "Capreolus capreolus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots()

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
# rationale: No subspecies noted by IUCN

## Cephalophus natalensis ####
# rationale: Two subspecies have been named: C. n. natalensis and C. n. robertsi (north of the Limpopo river)
# I only seem to have natalensis

curr.species <- "Cephalophus natalensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

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
tm_shape(clip_poly) + tm_polygons(alpha = 0) #+ tm_shape(sprp) + tm_polygons()

sprp_natalensis <- raster::intersect(sprp, clip_poly)
sprp_natalensis@data$subgroup <- "natalensis"

sprp_robertsi <- sprp - clip_poly
sprp_robertsi@data$subgroup <- "robertsi"

sprp_try <- raster::bind(sprp_natalensis, sprp_robertsi)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Cerdocyon thous ####
# rationale: wiki lists 5 subspecies: thous, azarae, entrerianus, aquilus, germanus
# I seem to have azarae and entrerianus

curr.species <- "Cerdocyon thous"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies", alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)
sprp_try@data[1, "subgroup"] <- "azarae x entrerianus"

sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

## Cervus canadensis ####
# rationale: I only have north american subspecies: canadensis, nannodes, roosevelti
# Not possible to separate each one based on subspecies so going by continent
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
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies", alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

clip_poly <- matrix(c(-125, 61,
                      -77,61,
                      -77,25,
                      -125,25,
                      -125,61),
                    ncol = 2, byrow = TRUE)

clip_poly <- Polygon(clip_poly)
clip_poly <- SpatialPolygons(list(Polygons(list(clip_poly), ID = "a")), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
tm_shape(clip_poly) + tm_polygons(alpha = 0) + tm_shape(sprp) + tm_polygons()

sprp_americas <- raster::intersect(sprp, clip_poly)
sprp_americas@data$subgroup <- "canadensis  x nannodes x roosevelti"

sprp_otherssp <- sprp - clip_poly
sprp_try <- raster::bind(sprp_americas, sprp_otherssp)
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

### Cervus elaphus ####
# rationale: I seem to have: elaphus, italicus
# and maybe montanus
# but not: barbarus, corsicanus, maral, brauneri
"Several subspecies of Western Red Deer have been recognized with their ranges as follows:
C. e. elaphus: Ireland, Great Britain,   continental Europe
C. e. barbarus: Atlas Mountains (Algeria, Tunisia)
C. e. corsicanus: Corsica (extinct, reintroduced in 1985), Sardinia
C. e. maral: Anatolia,
C. e. italicus: Italy (Ferrara)
C. e. brauneri Crimea (Russia)
C. e. montanus (syn. Carpathicus) Carpathian mountains"

# Not sure which subspecies are in Morocco or around Georgia/East Turkey, or south of Finland

"However, recent analyses call into question the veracity of the morphology-based subspecific taxonomy of C. elaphus 
more generally. Rather, three genetic lineages were clearly differentiated corresponding approximately with geographical 
factors. One lineage was distributed in central-western Europe (and also Ukraine), one from eastern Europe to the 
Middle East, and a third which corresponded to C. e. barbarus and C. e. corsicanus from North Africa and Sardinia, 
the latter two nomials are considered synonyms and were recommended considered as a subspecies"

# Removing North African polygons

curr.species <- "Cervus elaphus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies", alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 


IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

### Cervus nippon ####
# rationale:
"C. n. aplodontus (North Honshu Sika);
C. n. mantchuricus (Manchurian Sika): former USSR;
C. n. mandarinus (North China Sika) China (considered extinct);
C. n. grassianus, (Shansi Sika): China (considered extinct);
C. n. sichuanicus (Sichuan Sika): China;
C. n. yesoensis (Hokkaido Sika): Japan;
C. n. taiouanus (Formosan or Taiwan Sika): Taiwan;
C. n. pseudaxis (Viet Namese or Tonkin Sika): Viet Nam;
C. n. kopschi, (South China or Kopschi Sika): China;
C. n. keramae (Ryukyu or Kerama Sika): Japan;
C. n. pulchellus: Japan - Tsushima Islands."
# and C. n. "hortulorum of the mainland range."

# I have only Japanese samples from Honshu, Hokkaido, and Kyushu 
# Only Hokkaido sika (yesoensis) matches 
# Honshu samples are likely to be North Honshu sika (aplodontus)

curr.species <- "Cervus nippon"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies", alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

sprp@data[sprp$island == "Hokkaido", "subgroup"] <- "yesoensis"

## Chrysocyon brachyurus ####
# rationale: no subspecies described by IUCN or elsewhere

### Civettictis civetta #### 
# rationale: IUCN doesn't describe any subspecies but they are recognised elsewhere
# Not sure where to split polygon and I have multiple subspecies so will likely leave this one

curr.species <- "Civettictis civetta"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Conepatus chinga ####
# rationale: No subspecies described by IUCN or wiki, and range descriptions not available on msotw

## Connochaetes gnou ####
# rationale: No subspecies described by IUCN

## Connochaetes taurinus ####
# rationale: already included in species range polygons <3

## Crocuta crocuta ####
# rationale: IUCN doesn't list any subspecies and wiki states:
# "all the variation seen in the then recognised subspecies could also be found in a single population"

### Cynictis penicillata ####
# rationale: no subspecies listed by IUCN, and the 12 listed on msotw have no distribution descriptions

## Damaliscus lunatus ####
# rationale: already included in species range polygons <3
# Not sure where Zambian sample falls

## Damaliscus pygargus ####
# rationale: One subspecies already included, pygargus is the other described by IUCN
curr.species <- "Damaliscus pygargus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

IUCN_Native_Data@data <- IUCN_Native_Data@data %>%
  mutate(subgroup = case_when(binomial == curr.species & is.na(subspecies) ~ "pygargus",
                              TRUE ~ subgroup))
tm_shape(IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]) + tm_polygons("subgroup")

## Equus grevyi ####
# rationale: "However, Groves and Bell (2004) concluded that the species is indeed monotypic."

## Equus quagga ####
# rationale: "The molecular data represented a genetic cline" so it's monotypic

## Equus zebra ####
# rationale: "We continue to recognize Mountain Zebra as a single species comprising two subspecies."
# polygons already labelled <3

### Felis silvestris ####
# rationale: "A revised taxonomy of the Felidae. The final report of the Cat Classification Task Force of the IUCN/SSC Cat Specialist Group. Cat News Special Issue 11, 80 pp."
# two subspecies
# I seem to also have Felis libyca
curr.species <- "Felis silvestris"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Galictis cuja ####
# rationale: from wiki:
"Four subspecies are recognised:

Galictis cuja cuja – southwestern Bolivia, western Argentina, central Chile
Galictis cuja furax – southern Brazil, northeastern Argentina, Uruguay, and Paraguay
Galictis cuja huronax – south-central Bolivia, eastern Argentina
Galictis cuja luteola – extreme southern Peru, western Bolivia and northern Chile"
curr.species <- "Galictis cuja"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

### Genetta genetta ####
# rationale: 
"There is a high degree of intraspecific variation in this species, which has resulted in many described subspecies; 
the validity of many of these is unknown, while others might actually represent distinct species"
curr.species <- "Genetta genetta"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

## Genetta thierryi ####
# rationale: No subspecies listed by IUCN, wiki, or msotw
# plus it has a small continuous range

## Giraffa camelopardalis ####
# rationale: already recorded <3

### Herpestes ichneumon ####
# rationale: none listed on IUCN
curr.species <- "Herpestes ichneumon"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Hippopotamus amphibius ####
# rationale:
curr.species <- "Hippopotamus amphibius"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Hippotragus niger
# rationale: may not be able to subdivide
"Four subspecies are usually recognized: H. n. niger, H. n. kirkii, H. n. roosevelti and the isolated Giant Sable 
(H. n. variani) from Angola. As for many other antelope species, the validity and precise distribution of most of the 
described subspecies are uncertain. An extensive study of the geographical genetic structure of Hippotragus niger 
identified three genetic subdivisions representing a Kenya and east Tanzania clade (H. n. roosevelti), a west Tanzania 
clade (H. n. kirkii), and a southern African clade (H. n. niger) (Pitra et al. 2002)."

curr.species <- "Hippotragus niger"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons(alpha = 0) + tm_shape(sp.gmpd.points) + tm_dots() 

## Hyaena hyaena ####
# rationale: As of 2005,[3] no subspecies are recognised. (wiki)

## Kobus ellipsiprymnus ####
# rationale: already labelled <3

## Kobus kob ####
# rationale: already done <3

## Kobus leche ####
# rationale: dooone <

## Leopardus geoffroyi ####
# rationale: according to cat specialist group there's only one ssp

### Leopardus pardalis ####
# rationale: cat specialist group say there's two: pardalis (north) and mitis (south)
# "borders between subspecies are speculative" and I have one right on the border
curr.species <- "Leopardus pardalis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

## Leopardus tigrinus ####
# rationale: Until then L. tigrinus is recognised as having two subspecies: tigrinis (south) and oncilla
curr.species <- "Leopardus tigrinus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

sprp_try <- raster::disaggregate(sprp)
sprp_try@data$subgroup <- c("tigrinus", "oncilla")
sprp_try <- raster::aggregate(sprp_try, by = names(sprp_try))
tm_shape(sprp_try) + tm_polygons("subgroup")

IUCN_Native_Data <- raster::bind(IUCN_Native_Data[IUCN_Native_Data$binomial != curr.species,],
                                 sprp_try)

### Lontra canadensis ####
# rationale: (wiki)
"L. c. canadensis (Schreber, 1777) – (eastern Canada, U.S., Newfoundland)
L. c. kodiacensis (Goldman, 1935) – (Kodiak Island, Alaska)
L. c. lataxina (Cuvier, 1823) – (U.S.)
L. c. mira (Goldman, 1935) – (Alaska, British Columbia)
L. c. pacifica (J. A. Allen, 1898) – (Alaska, Canada, northern U.S., south to central California, northern Nevada, and northeastern Utah)
L. c. periclyzomae (Elliot, 1905) – (Queen Charlotte Islands, British Columbia)
L. c. sonora (Rhoads, 1898) – (U.S., Mexico)"
curr.species <- "Lontra canadensis"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Lutra lutra ####
# rationale:
"Lutra lutra angustifrons in North Africa;
Lutra lutra aurobrunneus in Garhwal Himalayas in northern India and higher altitudes in Nepal;
Lutra lutra barang in southeast Asia (Thailand, Viet Nam, Indonesia and Sumatra);
Lutra lutra chinensis in southern China and Taiwan;
Lutra lutra hainana in Hainan Island, China;
Lutra lutra kutab in northern India (Kashmir);
Lutra lutra lutra is the most widely distributed, spanning from Portugal to South Korea;
Lutra lutra meridionalis in from Georgia through Armenia, Azerbaijan and Iran
Lutra lutra monticolus in northern India (Punjab, Kumaon, Himachal Pradesh, Sikkim and Assam) Nepal, Bhutan and Myanmar;
Lutra lutra nair in southern India and Sri Lanka;
Lutra lutra seistanica in Afghanistan, Eastern Iran, Kazakhstan, Uzbekistan, and Turkmenistan; and
the Japanese subspecies, Lutra lutra whiteleyi, was considered a distinct species (L. nippon) by Suzuki et al. (1996)."

curr.species <- "Lutra lutra"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Lycalopex culpaeus ####
# five subspecies in wiki
curr.species <- "Lycalopex culpaeus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

## Lycalopex fulvipes ####
# rationale: no subspecies recorded by iucn, wiki, or msotw

### Lycalopex gymnocercus ####
# rationale: (wiki)
"Five subspecies are currently recognised, although the geographic range of each is unclear, and the type localities 
of three of them lie outside the present-day range of the species:[1][4]"
curr.species <- "Lycalopex gymnocercus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

### Lycaon pictus ####
# rationale: five subspecies in wiki
curr.species <- "Lycaon pictus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

## Lynx canadensis ####
# rationale: cat group says "Therefore we conclude that Lynx canadensis is a monotypic species:"

### Lynx lynx ####
# rationale: On the basis of current evidence we propose the following six subspecies" (cat group)
# and they have a map 
curr.species <- "Lynx lynx"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Lynx pardinus
curr.species <- "Lynx pardinus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Lynx rufus
curr.species <- "Lynx rufus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Martes americana
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Martes foina
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Martes martes
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Martes melampus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Martes pennanti
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Meles meles
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Melogale moschata
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mephitis mephitis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mungos mungo
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mustela erminea
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mustela lutreola
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mustela nivalis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Mustela putorius
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Nasua nasua
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Neovison vison
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Nyctereutes procyonoides
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Odocoileus hemionus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Odocoileus virginianus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Oreamnos americanus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Otocolobus manul
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Otocyon megalotis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ourebia ourebi
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ovibos moschatus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ovis ammon
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ovis canadensis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ovis dalli
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ozotoceros bezoarticus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Paguma larvata
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Panthera leo
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Panthera onca
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Panthera pardus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Pecari tajacu
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Pelea capreolus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Phacochoerus aethiopicus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Philantomba monticola
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Procapra gutturosa
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Procyon lotor
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Procyon pygmaeus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Puma concolor
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Rangifer tarandus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Raphicerus campestris
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Redunca arundinum
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Redunca fulvorufula
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Rupicapra pyrenaica
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Rupicapra rupicapra
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Rusa unicolor
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Saiga tatarica
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Spilogale gracilis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Spilogale putorius
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Suricata suricatta
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Sus scrofa
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Sylvicapra grimmia
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Syncerus caffer
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Taxidea taxus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus angasii
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus eurycerus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus oryx
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus scriptus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus spekii
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Tragelaphus strepsiceros
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Urocyon cinereoargenteus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Urocyon littoralis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ursus americanus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ursus arctos
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Ursus maritimus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Vulpes corsac
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Vulpes lagopus
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Vulpes macrotis
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Vulpes velox
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 

# Vulpes vulpes
curr.species <- ""
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
tm_shape(sprp) + tm_polygons("subspecies") + tm_shape(sp.gmpd.points) + tm_dots() 









# tidying up
rm(list = ls(pattern = "^sprp"))

## discard pile ####
curr.species <- "Ovibos moschatus"
sprp <- IUCN_Native_Data[IUCN_Native_Data$binomial == curr.species, ]
sprp <- IUCN_Mammals[IUCN_Mammals$binomial == curr.species, ]

sp.gmpd.points <- SpatialPoints(GMPD_Data[GMPD_Data$HostCorrectedName == curr.species, c("Longitude", "Latitude")])
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