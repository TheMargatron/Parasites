# Plotting analysis of distance metrics (without climate)
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(MCMCglmm)
library(here)
library(postMCMCglmm)
library(beepr)
library(broom.mixed)
library(maps)
library(rgdal)
extrafont::loadfonts(device = "win")
source(here::here("Functions.R"))
# source(here::here("03.1Analysing distance metrics"))

GMPD_Distances_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), header = TRUE, stringsAsFactors = FALSE)

# Prepping data ####
GMPD_IUCN_Species <- GMPD_Distances_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected,
         HostSubgroup           = paste(HostCorrectedName, subgroup))

set.seed(75114)
# 7b: Lat + Med random :)* ####
Mc_IUCN_Species_L_M_Hb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random  = ~ us(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName,
                                   data    = GMPD_IUCN_Species,
                                   family  = "multinomial2",
                                   nitt    = 130000,
                                   thin    = 100,
                                   burnin  = 30000,
                                   pr      = TRUE)
M1_summary <- summary(Mc_IUCN_Species_L_M_Hb)
# DIC:        561086.3 
# Intercept:  -0.445467 **
# Latitude:   -0.521344 ***
# MedianProp: 0.205312 .
Sol1 <- autocorr(Mc_IUCN_Species_L_M_Hb$Sol, lags = 1)
VCV1 <- autocorr(Mc_IUCN_Species_L_M_Hb$VCV, lags = 1)
# TODO: iterations
#write_rds(Mc_IUCN_Species_L_M_Hb, here::here("BackupModel.rds"))


## plots ####
plot_dat <- GMPD_IUCN_Species

medIntercept <- summary(Mc_IUCN_Species_L_M_Hb)$solutions["(Intercept)","post.mean"]
medLatitude <- summary(Mc_IUCN_Species_L_M_Hb)$solutions["LatitudeScaled","post.mean"]
medMedian <- summary(Mc_IUCN_Species_L_M_Hb)$solutions["MedianPropScaled","post.mean"]

meanLat <- mean(abs(plot_dat$Latitude))
meanMed <- mean(plot_dat$MedianProp)

# bigline
latrange <- seq(min(abs(plot_dat$Latitude)), max(abs(plot_dat$Latitude)), length.out = 2000)
medrange <- seq(min(plot_dat$MedianProp), max(plot_dat$MedianProp), length.out = 2000)

latrange.sc <- (latrange - meanLat)/(sd(abs(plot_dat$Latitude) - meanLat))
medrange.sc <- (medrange - meanMed)/(sd(plot_dat$MedianProp - meanMed))

bigline <- data.frame("latrange" = latrange, "medrange" = medrange, "latrange.sc" = latrange.sc, "medrange.sc" = medrange.sc)

logitlat <- medIntercept + (mean(plot_dat$MedianPropScaled)*medMedian) + (bigline$latrange.sc*medLatitude)
logitpol <- medIntercept + (mean(plot_dat$LatitudeScaled)*medLatitude) + (bigline$medrange.sc*medMedian)
bigline$linelat <- 1/(1+exp(-logitlat))
bigline$linepol <- 1/(1+exp(-logitpol))

# small lines
plotranef <- broom.mixed::tidy(Mc_IUCN_Species_L_M_Hb, effects="ran_vals")
plotranef <- plotranef %>% 
  select(-std.error) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  as.data.frame()

plot_dat$smallatlines <- NA
plot_dat$smalmedlines <- NA

for(i in 1:nrow(plotranef)){
  p.rows <- which(plot_dat$HostCorrectedName == plotranef[i, "level"])
  
  c.logitlat <- medIntercept + plotranef[i, "(Intercept)"]+ (meanMed*(medMedian + plotranef[i, "MedianPropScaled"])) + (plot_dat[p.rows,"LatitudeScaled"]*(medLatitude + plotranef[i, "LatitudeScaled"]))
  c.logitmed <- medIntercept + plotranef[i, "(Intercept)"]+ (meanLat*(medLatitude + plotranef[i, "LatitudeScaled"])) + (plot_dat[p.rows,"MedianPropScaled"]*(medMedian + plotranef[i, "MedianPropScaled"]))
  
  plot_dat[p.rows,"smallatlines"] <- 1/(1 + exp(-c.logitlat))
  plot_dat[p.rows,"smalmedlines"] <- 1/(1 + exp(-c.logitmed))
}


### Latitude base ####
ggplot(plot_dat, aes(Latitude, Prevalence)) +
  geom_line(data = bigline, aes(y = linelat, x = latrange),
            col = "#940039", linewidth = 2, lineend = "round") +
  theme_void() +
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
        ) +
  labs(x = "Latitude", y = "Parasitism rate")+
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), limits = c(0,90)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) 


### Latitude host ####
ggplot(plot_dat, aes(Latitude, Prevalence)) +
  geom_line(aes(y = smallatlines, x = Latitude, col = HostCorrectedName), size = 1) +
  geom_line(data = bigline, aes(y = linelat, x = latrange),
            col = "#940039", linewidth = 2, lineend = "round") +
  theme_void() +
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  labs(x = "Latitude", y = "Parasitism rate")+
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), limits = c(0,90)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  scale_color_manual(values = colorRampPalette(c("#60A561","#60A561"))(94))


### Median prop base ####
ggplot(plot_dat, aes(MedianProp, Prevalence)) +
  geom_line(data = bigline, aes(y = linepol, x = medrange),
            col = "#940039", linewidth = 2, lineend = "round") +
  theme_void() +
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  labs(x = "Proximity to range margin", y = "Parasitism rate") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) 


### Median prop host ####
ggplot(plot_dat, aes(MedianProp, Prevalence)) +
  geom_line(aes(y = smalmedlines, x = MedianProp, col = HostCorrectedName), size = 1) +
  geom_line(data = bigline, aes(y = linepol, x = medrange),
            col = "#940039", linewidth = 2, lineend = "round") +
  theme_void() +
  theme(axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(angle = 90, 
                                    margin = margin(t=0,r=15,b=0,l=0)),
        axis.line = element_line(linewidth = 2, colour = "#656565", lineend = "round"),
        axis.ticks = element_line(linewidth = 1, colour = "#656565", lineend = "round"),
        axis.ticks.length = unit(0.3, "lines"),
        axis.text.x = element_text(size = 15, margin = margin(t=5,r=0,b=0,l=0)),
        axis.text.y = element_text(size = 15, margin = margin(t=0,r=5,b=0,l=0)),
        
        plot.background = element_rect(fill = "#FFECCD", colour = "#FFECCD"),
        plot.margin = margin(t=20,r=25,b=10,l=20),
        text = element_text(family = "Outfit", size = 30),
        legend.position = "none",
        aspect.ratio = 0.7
  ) +
  labs(x = "Proximity to range margin", y = "Parasitism rate") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  scale_color_manual(values = colorRampPalette(c("#60A561","#60A561"))(94))


#### 9b: subgroup full ####
host_and_subgroup2 <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random = ~ idh(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName + idh(LatitudeScaled + MedianPropScaled + 1):HostSubgroup,
                                   data   = GMPD_IUCN_Species,
                                   family = "multinomial2",
                                   nitt   = 130000,
                                   thin   = 100,
                                   burnin = 30000,
                                   pr     = TRUE
                               )
M2_summary <- summary(host_and_subgroup2)


#### 9d: host idh ####
host_not_subgroup4 <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               random = ~ idh(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               nitt   = 130000,
                               thin   = 100,
                               burnin = 30000,
                               pr     = TRUE
                               )
M3_summary <-summary(host_not_subgroup4)


# Aesthetic plots ####
# read in mammals from gbif
# apply raster to count layers in each cell
# ggplot raster

IUCN_Mammals <- readOGR(here::here("Data/IUCN"), "MAMMALS") # plenty time to make a cup of tea
IUCN_Terrestrial <- IUCN_Mammals[IUCN_Mammals$terrestial == "true",]
write_rds(IUCN_Terrestrial, here::here("Terrestrial_Mammals.rds"))
mammal_count_raster <- raster::rasterize(x = IUCN_Mammals, 
                                         y = raster(),
                                         field = 1,
                                         fun = "count" 
                                         ) 
mammal_count_raster_terr <- raster::rasterize(x = IUCN_Terrestrial, 
                                         y = raster(),
                                         field = 1,
                                         fun = "count" 
) 

world <- map_data("world")

ggplot() +  
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region) # change map id, maybe make blank variable called "id" 
  ) #+
  geom_raster(
    data = fortify(mammal_count_raster), 
    aes(fill = layer, x=x, y=y), # might have to change value to count, deleted x=x and y=y
    alpha = 0.8
    ) + 
  #coord_equal() +
  #theme_map() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank()
        ) +
  scale_fill_gradient(low = "#E4C1F9", high = "#940039")


## using tmap ####
library(tmap)
data("World")

pal <- colorRampPalette(c("#E4C1F9","#940039"))
tm_shape(world) + tm_polygons(col = "#656565") #+

tm_shape(World) + tm_polygons(col = "#656565") + 
tm_shape(mammal_count_raster_terr) + tm_raster(col = "layer", palette = pal(201),
                                          n = 201, legend.show = FALSE) +
tm_layout(bg.color = "#FFECCD")


write_rds(mammal_count_raster_terr, here::here("BackupMapTerr.rds"))


masked <- mask(x = mammal_count_raster_terr, mask = World)
equatorline <- SpatialLines(list(Lines(Line(cbind(c(-180, 180),c(0,0))), ID="a")))


tm_shape(World) + tm_polygons(col = "#656565") + 
  tm_shape(masked) + tm_raster(col = "layer", palette = pal(201),
                                                 n = 201, legend.show = FALSE) +
  tm_shape(equatorline) + tm_lines(col = "#656565", lty = "dashed") +
  tm_layout(bg.color = "#FFECCD") 
