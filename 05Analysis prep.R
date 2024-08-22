# Analysis and plotting of final dataset
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ###########################################################

library(here)               #
library(tidyverse)          #
library(corrplot)
library(tmap)
library(patchwork)
library(cowplot)
library(lme4)
library(DHARMa)
library(grid)
library(gridExtra)
library(cowplot)

source(here::here("Functions.R"))

GMPD_Climate_Data <- read.csv(here::here("Data/Data back ups/GMPD_Climate_Data_04.csv"), 
                              header = TRUE, 
                              stringsAsFactors = FALSE) 

plot_outputs <- FALSE
rerun_models <- FALSE

# Prep and data exploration ####################################################

# https://stats.stackexchange.com/questions/117497/comparing-between-random-effects-structures-in-a-linear-mixed-effects-model
# Can't compare AIC as I usually would for random effects so could use anova instead
# https://stats.stackexchange.com/questions/418593/understanding-anova-to-compare-mixed-model-with-a-gzlm

# Can't use anova because there is no difference in the number of parameters between models
# Could then return to comparing with AIC? The reason it's not used for random effects is 
# because it penalises according to the number of parameters and that's not the point for random
# But then Phylum will always come out on top just because it has more levels
# Turns out I can use anova as long as I'm comparing nested random effects because they do count as terms

GMPD_Analysis_Data %>% 
  dplyr::filter(RestrAll, CleanAll) %>% 
  group_by(ParPhylum, ParType) %>% 
  summarise(n = n()) %>% 
  print() %>% 
  group_by(ParType) %>% 
  summarise(n = sum(n))

# The argument for Phylum is that it has more levels and so captures variation more accurately
# The flipside is true of Type which is coarser
# The argument against Phylum is that the levels are very uneven in their sampling
# That's where I think Type wins out, because the smaller samples coalesce without the
# groupings losing meaning as they are based on (hopefully) informative traits

## Correlations ####
test_data <- GMPD_Analysis_Data %>% 
  filter(RestrAll,
         CleanAll) %>% 
  mutate(LatitudeScaled    = scale(abs(Latitude)),
         MedianPropSquared = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                              TRUE       ~ MedianProp_iucn),
         RowID                   = row.names(.)) %>% 
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp))

if(plot_outputs){
  test_data %>% 
    dplyr::select(Prevalence, SampleSize, Latitude, 
                  MedianProp_iucn, MedianProp_gbif,
                  CorrectedDensity, CorrectedDistance, CorrectedDistProp, CorrectedAngle,
                  Axis1, Axis2, 
                  wc2.1_10m_bio_5, wc2.1_10m_bio_6, wc2.1_10m_bio_12) %>% 
    cor() %>% 
    corrplot::corrplot.mixed(diag = "n", tl.pos = "lt")
  
  corr_plot <- test_data %>% 
    dplyr::rename(`Range position` = MedianProp_iucn,
                  `Niche position` = CorrectedDistProp) %>% 
    dplyr::select(Prevalence, SampleSize, Latitude,
                  `Range position`, `Niche position`) %>% 
    cor() %>% .[,5:1] %>% 
    ggcorrplot(outline.col = "white", 
               colors = c("#222E50", "white", "#D1495B"))  +
    theme(text = element_text(family = "Outfit", size = 15)) + 
    geom_text(aes(label = value), family = "Outfit")
  
  ggsave(here::here("Figures/correlation plot.pdf"), corr_plot, width = 10, height = 9, 
         device = cairo_pdf)
  
}

# https://doi.org/10.1111/j.1600-0587.2012.07348.x
# None of the variables I plan on modelling 
# (Latitude, MedianProp, CorrectedDistance)
# have particularly high collinearity (above 0.7)


## Non-normal distributions ####
# Latitude
hist(test_data$Latitude)
hist(abs(test_data$Latitude))
if(plot_outputs){
  Latitude_hist <- test_data %>% 
    pretty_hist(xvar = Latitude, 
                breaks = seq(-80, 80, length.out = 40)) +
    scale_x_continuous(limits = c(-80,80), 
                       breaks = c(-80,-40,0,40,80)) +
    labs(x = "Absolute latitude", y = "") +
    annotate(geom = "text", x = -80, y = 0.85*600,
             label = "A",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 600) 
  
  abslatitude_hist <- test_data %>% 
    pretty_hist(xvar = abs(Latitude), 
                breaks = seq(0, 80, length.out = 40)) +
    labs(x = "Absolute latitude", y = "") +
    annotate(geom = "text", x = 0, y = 0.85*600,
             label = "B",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 600) 
  
  Latitude_hists <- Latitude_hist +
    abslatitude_hist +
    plot_layout(design = "
              AB
              ") 
  
  ggsave(here::here("Figures/Latitude hists.pdf"), Latitude_hists, width = 10, height = 4,
         device = cairo_pdf)
}

# MedianProp
hist(test_data$MedianProp_iucn)
test_data %>% 
  dplyr::mutate(MedianPropSquared_iucn = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                                     TRUE       ~ MedianProp_iucn)) %>% 
  dplyr::pull(MedianPropSquared_iucn) %>% 
  hist()

hist(test_data$MedianProp_gbif)
test_data %>% 
  dplyr::mutate(MedianPropSquared_gbif = dplyr::case_when(!AboveMedn_gbif ~ MedianProp_gbif * -1,
                                                          TRUE       ~ MedianProp_gbif)) %>% 
  dplyr::pull(MedianPropSquared_gbif) %>% 
  hist()

if(plot_outputs){
  MedianProp_hist <- test_data %>% 
    pretty_hist(xvar = MedianProp_iucn, 
                breaks = seq(0, 1, length.out = 20)) +
    labs(x = "Range position", y = "") +
    annotate(geom = "text", x = 0, y = 0.85*600,
             label = "C",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 600) 
  
  MedianProp_s_hist <- test_data %>% 
    dplyr::mutate(MedianPropSquared_iucn = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                                            TRUE       ~ MedianProp_iucn)) %>% 
    pretty_hist(xvar = MedianPropSquared_iucn, 
                breaks = seq(-1, 1, length.out = 20)) +
    labs(x = "Range position", y = "") +
    annotate(geom = "text", x = -1, y = 0.85*600,
             label = "D",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 600) 
  
  MedianProp_hists <- MedianProp_hist +
    MedianProp_s_hist +
    plot_layout(design = "
              AB
              ") 
  
  ggsave(here::here("Figures/Range position hists.pdf"), MedianProp_hists, width = 10, height = 4,
         device = cairo_pdf)
}

# CorrectedDistance
hist(test_data$CorrectedDistance)
hist(test_data$CorrectedDistProp)
hist(log(test_data$CorrectedDistance+1))
hist(sqrt(test_data$CorrectedDistance))
test_data %>% 
  dplyr::arrange(desc(pmax(CorrectedDistance, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared = case_when(!TopHalf ~ CorrectedDistance * -1,
                                                     TRUE ~ CorrectedDistance)) %>% 
  dplyr::pull(CorrectedDistanceSquared) %>% 
  hist()
test_data %>% 
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp)) %>% 
  dplyr::pull(CorrectedDistPropSquared) %>% 
  hist()

if(plot_outputs){
  DistProp_hist <- test_data %>% 
    pretty_hist(xvar = CorrectedDistProp, 
                breaks = seq(0, 1, length.out = 20)) +
    labs(x = "Niche position", y = "") +
    annotate(geom = "text", x = 0, y = 0.85*900,
             label = "E",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 900)
  
  DistProp_s_hist <- test_data %>% 
    dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>% 
    dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                  CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                       TRUE ~ CorrectedDistProp)) %>% 
    pretty_hist(xvar = CorrectedDistPropSquared, 
                breaks = seq(-1, 1, length.out = 20)) +
    labs(x = "Niche position", y = "") +
    annotate(geom = "text", x = -1, y = 0.85*900,
             label = "F",
             family = "Outfit", size = 5, hjust = 0) +
    ylim(0, 900)
  
  DistProp_hists <- DistProp_hist +
    DistProp_s_hist +
    plot_layout(design = "
              AB
              ") 
  
  ggsave(here::here("Figures/Niche position hists.pdf"), DistProp_hists, width = 10, height = 4,
         device = cairo_pdf)

  all_hists <- Latitude_hist + abslatitude_hist +
    MedianProp_hist + MedianProp_s_hist +
    DistProp_hist + DistProp_s_hist +
    plot_layout(design = "
              AB
              CD
              EF
              ") 

  ggsave(here::here("Figures/all hists.pdf"), all_hists, width = 10, height = 9, 
         device = cairo_pdf)
  
}

# I think the last is the most appropriate given that it's a true half-normal 
# distribution that was calculated from a 2d space
# Or at least, something approximating this approach
# Debatable whether I should use angle in sorting


rm(list = ls(pattern = "test_"))

## Prepping Species data ####
GMPD_Both_Species <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll,
                CleanAll) %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                MedianPropScaled_iucn        = base::scale(MedianProp_iucn),
                MedianPropSquared_iucn       = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                                                TRUE       ~ MedianProp_iucn),
                MedianPropSquScaled_iucn     = base::scale(MedianPropSquared_iucn),
                
                MedianPropScaled_gbif   = base::scale(MedianProp_gbif),
                MedianPropSquared_gbif  = dplyr::case_when(!AboveMedn_gbif ~ MedianProp_gbif * -1,
                                                           TRUE       ~ MedianProp_gbif),
                MedianPropSquScaled_gbif = base::scale(MedianPropSquared_gbif),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                CorrectedDensityScaled_5  = base::scale(CorrectedDensity_5),
                CorrectedDistanceScaled_5 = base::scale(CorrectedDistance_5),
                CorrectedAngleScaled_5    = base::scale(CorrectedAngle_5),
                
                CorrectedDensityScaled_20  = base::scale(CorrectedDensity_20),
                CorrectedDistanceScaled_20 = base::scale(CorrectedDistance_20),
                CorrectedAngleScaled_20    = base::scale(CorrectedAngle_20),
                
                ParasiteDetected        = as.integer(round(SampleSize * Prevalence, 0)),
                ParasiteUndetected      = SampleSize - ParasiteDetected,
                
                # HostSubgroup            = paste(HostCorrectedName, subgroup, sep = "_"),
                HostParasite            = paste(HostCorrectedName, ParasiteCorrectedName, sep = "_"),
                
                RowID                   = row.names(.)) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistance, CorrectedAngle))) %>% 
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared = case_when(!TopHalf ~ CorrectedDistance * -1,
                                                     TRUE ~ CorrectedDistance),
                CorrectedDistanceSquScaled = base::scale(CorrectedDistanceSquared)) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistance_5, CorrectedAngle_5))) %>%
  dplyr::mutate(TopHalf_5 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared_5 = case_when(!TopHalf ~ CorrectedDistance_5 * -1,
                                                       TRUE ~ CorrectedDistance_5),
                CorrectedDistanceSquScaled_5 = base::scale(CorrectedDistanceSquared_5)) %>%
  dplyr::arrange(desc(pmax(CorrectedDistance_20, CorrectedAngle_20))) %>%
  
  dplyr::mutate(TopHalf_20 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistanceSquared_20 = case_when(!TopHalf ~ CorrectedDistance_20 * -1,
                                                        TRUE ~ CorrectedDistance_20),
                CorrectedDistanceSquScaled_20 = base::scale(CorrectedDistanceSquared_20)) %>% 
  
  # Proportional distance
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>%
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp),
                CorrectedDistPropSquScaled = base::scale(CorrectedDistPropSquared)) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistProp_5, CorrectedAngle_5))) %>%
  dplyr::mutate(TopHalf_5 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_5 = case_when(!TopHalf ~ CorrectedDistProp_5 * -1,
                                                       TRUE ~ CorrectedDistProp_5),
                CorrectedDistPropSquScaled_5 = base::scale(CorrectedDistPropSquared_5)) %>%
  dplyr::arrange(desc(pmax(CorrectedDistProp_20, CorrectedAngle_20))) %>%
  
  dplyr::mutate(TopHalf_20 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_20 = case_when(!TopHalf ~ CorrectedDistProp_20 * -1,
                                                        TRUE ~ CorrectedDistProp_20),
                CorrectedDistPropSquScaled_20 = base::scale(CorrectedDistPropSquared_20))


GMPD_Both_Species %>% 
  dplyr::select(Prevalence, SampleSize, Latitude, 
                MedianProp_iucn, MedianPropSquared_iucn,
                MedianProp_gbif, MedianPropSquared_gbif,
                CorrectedDistance, CorrectedDistanceSquared,
                CorrectedDistProp, CorrectedDistPropSquared) %>% 
  mutate(abslatitude = abs(Latitude)) %>% 
  cor() %>% 
  corrplot::corrplot.mixed(diag = "n", tl.pos = "lt")


# Analysis  ####################################################################
## Frequentist ####

### Model naming convention
## Model main lists:
# LM = Frequentist rather than Bayesian
# IUCN/GBIF = method used for restricting data and calculating position within range
LM_IUCN_Species    <- list()
LM_GBIF_Species    <- list()

## Model sub lists
# Fixed         = Fixed effects only
## Fixed_Grp    = Fixed effects only, with host group also a fixed effect
## Fixed_Typ    = Fixed effects only, with parasite type also a fixed effect

# Host          = Host species random effect
## HostGroup    = Host species nested within host group as random effects
## Group        = Host group as a fixed effect and host species as a random effect

# Parasite      = Parasite species as random effect
## Type         = Parasite type as random effect
## ParType      = Parasite species nested within Parasite type as random effects
## TypePar      = Parasite type as a fixed effect and parasite species as a random effect

# Both
## BothType     = Host species and parasite type as random effects
## BothSpecies  = Host species and parasite species as random effects
## CrossType    = Host species crossed with parasite type as random effects
## CrossSpecies = Host species crossed with parasite species as random effects


if(rerun_models){
  
  LM_IUCN_Species <- all_models(GMPD_Both_Species, method = "iucn")
  saveRDS(LM_IUCN_Species, here::here("Data/Model back ups/LM_IUCN_Species_05.rds"))
  
  LM_GBIF_Species <- all_models(GMPD_Both_Species, method = "gbif")
  saveRDS(LM_GBIF_Species, here::here("Data/Model back ups/LM_GBIF_Species_05.rds"))
  
} else {
  
  LM_IUCN_Species <- readRDS(here::here("Data/Model back ups/LM_IUCN_Species_05.rds"))
  LM_GBIF_Species <- readRDS(here::here("Data/Model back ups/LM_GBIF_Species_05.rds"))
  
}
### Random effects quickly ####
# Doing a quick comparison of hosts vs parasites
anova(LM_IUCN_Species$Host$Null, LM_IUCN_Species$Fixed$Null)
# host better than nothing

anova(LM_IUCN_Species$Group$Null, LM_IUCN_Species$Fixed$Null)
anova(LM_IUCN_Species$HostGroup$Null, LM_IUCN_Species$Host$Null)
anova(LM_IUCN_Species$HostGroup$Null, LM_IUCN_Species$Group$Null)
# group better than nothing, but doesn't add anything on top of host
# host adds things on top of group

anova(LM_IUCN_Species$Parasite$Null, LM_IUCN_Species$Fixed$Null)
# parasite better than nothing

anova(LM_IUCN_Species$Type$Null, LM_IUCN_Species$Fixed$Null)
anova(LM_IUCN_Species$ParType$Null, LM_IUCN_Species$Parasite$Null)
anova(LM_IUCN_Species$ParType$Null, LM_IUCN_Species$Type$Null)
# Type better than nothing, and adds a little on top of parasite
# parasite adds things on top of type

anova(LM_IUCN_Species$BothSpecies$Null, LM_IUCN_Species$Parasite$Null)
anova(LM_IUCN_Species$BothSpecies$Null, LM_IUCN_Species$Host$Null)
# Host adds more on top of parasite, parasite adds more on top of host

anova(LM_IUCN_Species$CrossSpecies$Null, LM_IUCN_Species$BothSpecies$Null)
# crossing host and parasite doesn't improve? AIC is lower, but P is 1

# will compare fixed effects with BothSpecies
### Fixed effects ####

#### Geographic ####

#### Latitude
# ______________________________________________________________________________
# starting with Latitude because that's what I'm most interested in
# Latitude is accepted
compare_models(LM_IUCN_Species$BothSpecies_IO$Lat, 
               LM_IUCN_Species$BothSpecies$Null,
               "IUCN data, Latitude vs null, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$Lat, 
               LM_GBIF_Species$BothSpecies$Null,
               "GBIF data, Latitude vs null, random both species")


#### MedianProp
# ______________________________________________________________________________
# Next thing we're interested in is position within geographic range
# MedianProp is accepted in iucn but not in gbif
# Either way we want need to check it as a quadratic

compare_models(LM_IUCN_Species$BothSpecies_IO$LatMed,
               LM_IUCN_Species$BothSpecies_IO$Lat,
               "IUCN data, MedianProp added to Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatMed,
               LM_GBIF_Species$BothSpecies_IO$Lat,
               "GBIF data, MedianProp added to Latitude, random both species")

#### Quadratic MedianProp
# ______________________________________________________________________________
# Probably better to model it as a quadratic
# Definitely better as a quadratic, and included in both iucn and gbif
# Got the same thing of it switching signs, but same conclusion as previous

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQMed,
               LM_IUCN_Species$BothSpecies_IO$LatMed,
               "IUCN data, quadratic added to MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQMed,
               LM_GBIF_Species$BothSpecies_IO$LatMed,
               "GBIF data, quadratic added to MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQMed,
               LM_GBIF_Species$BothSpecies_IO$Lat,
               "GBIF data, quadratic MedianProp added to Latitude, random both species")

#### Asymmetric Quadratic MedianProp
# ______________________________________________________________________________
# Median prop might have an asymmetric pattern depending on which half of the range it's in
# fails to converge, but it would be included otherwise

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQMedAsym,
               LM_IUCN_Species$BothSpecies_IO$LatQMed,
               "IUCN data, asymmetry added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQMedAsym,
               LM_GBIF_Species$BothSpecies_IO$LatQMed,
               "GBIF data, asymmetry added to quadratic MedianProp and Latitude, random both species")



#### Latitude * Quadratic MedianProp
# ______________________________________________________________________________
# Testing out the interaction term 
# There is support for an interaction

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQMedInt,
               LM_IUCN_Species$BothSpecies_IO$LatQMed,
               "IUCN data, interaction added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQMedInt,
               LM_GBIF_Species$BothSpecies_IO$LatQMed,
               "GBIF data, interaction added to quadratic MedianProp and Latitude, random both species")

#### Climatic 10 ####
# ______________________________________________________________________________
#### Distance
# ______________________________________________________________________________
# Distance from niche centre
# yes
# But it should be modelled as a quadratic

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQIMedProp,
               LM_IUCN_Species$BothSpecies_IO$LatQMedInt,
               "IUCN data, NicheDist added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQIMedProp,
               LM_GBIF_Species$BothSpecies_IO$LatQMedInt,
               "GBIF data, NicheDist added to quadratic MedianProp and Latitude, random both species")


#### Quadratic Distance
# ______________________________________________________________________________
# Distance from niche centre as quadratic
# retain

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQIMedQProp,
               LM_IUCN_Species$BothSpecies_IO$LatQMedInt,
               "IUCN data, quadratic NicheDist added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQIMedQProp,
               LM_GBIF_Species$BothSpecies_IO$LatQMedInt,
               "GBIF data, quadratic NicheDist added to quadratic MedianProp and Latitude, random both species")

#### Latitude * (Quadratic MedianProp + Quadratic Distance)
# ______________________________________________________________________________
# Distance from niche centre interacting with latitude 

compare_models(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp,
               LM_IUCN_Species$BothSpecies_IO$LatQIMedQProp,
               "IUCN data, interaction added to quadratic NicheDist, quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$BothSpecies_IO$LatQIMedQIProp,
               LM_GBIF_Species$BothSpecies_IO$LatQIMedQProp,
               "GBIF data, interaction added to quadratic NicheDist, quadratic MedianProp and Latitude, random both species")


### Fixed effects Conclusions ####
# Final fixed effects model
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)

# in comparison to gbif there are some different effect size estimates...
# Nothing too drastic but they just look a little off
summary(LM_GBIF_Species$BothSpecies_IO$LatQIMedQIProp)
plot_latitude(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species)
plot_latitude(model_list = LM_GBIF_Species, model_data = GMPD_Both_Species)

plot_medianprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 25)
plot_medianprop(model_list = LM_GBIF_Species, model_data = GMPD_Both_Species, fixed_lat = 25)

plot_distprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 25)
plot_distprop(model_list = LM_GBIF_Species, model_data = GMPD_Both_Species, fixed_lat = 25)

# basically the same across raster resolutions 
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp_5)
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp_20)


plot_latitude(model_list = LM_IUCN_Species, 
                model_name = "LatQIMedQIProp_5", 
                model_data = GMPD_Both_Species,
                raster_res = "_5")
plot_latitude(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species)
plot_latitude(model_list = LM_IUCN_Species, 
                 model_name = "LatQIMedQIProp_20", 
                 model_data = GMPD_Both_Species,
                 raster_res = "_20")

plot_medianprop(model_list = LM_IUCN_Species, 
                model_name = "LatQIMedQIProp_5", 
                model_data = GMPD_Both_Species, 
                fixed_lat = 40,
                raster_res = "_5")
plot_medianprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40)
plot_medianprop(model_list = LM_IUCN_Species, 
                model_name = "LatQIMedQIProp_20", 
                model_data = GMPD_Both_Species, 
                fixed_lat = 40,
                raster_res = "_20")

plot_distprop(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIProp_5", 
              model_data = GMPD_Both_Species, 
              fixed_lat = 40,
              raster_res = "_5")
plot_distprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40)
plot_distprop(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIProp_20", 
              model_data = GMPD_Both_Species, 
              fixed_lat = 40,
              raster_res = "_20")

# diagnostics
sim_LatQIMedQIProp_BothSpecies    <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_BothSpecies)



### Random Effects double check and diagnostics ####
# When looking at the random effects the maximal fixed effects model was:
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)
# This repeatedly failed to converge when we tried to model random slopes so we stick to random intercepts


#### Host focus ####
##### Host
summary(LM_IUCN_Species$Host_IO$LatQIMedQIProp)

anova(LM_IUCN_Species$Host_IO$LatQIMedQIProp, LM_IUCN_Species$Fixed$LatQIMedQIProp)
# Host is definitely better than nothing 

sim_LatQIMedQIProp_Host_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Host_IO$LatQIMedQIProp)
sim_LatQIMedQIProp_Fixed <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Fixed)
testOverdispersion(sim_LatQIMedQIProp_Fixed)

plot(sim_LatQIMedQIProp_Host_IO)
testOverdispersion(sim_LatQIMedQIProp_Host_IO)
# Pretty strong residual trend with host

##### Group
{
summary(LM_IUCN_Species$Group_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$Group_IO$LatQIMedQIProp, LM_IUCN_Species$Fixed$LatQIMedQIProp)
# Also better than nothing
sim_LatQIMedQIProp_Group_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Group_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_Group_IO)
# Doesn't have the strong trend, but it's also very overdispersed

##### Group/Host
summary(LM_IUCN_Species$HostGroup_IO$LatQIMedQIProp)
# 
anova(LM_IUCN_Species$HostGroup_IO$LatQIMedQIProp, LM_IUCN_Species$Host_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$HostGroup_IO$LatQIMedQIProp, LM_IUCN_Species$Group_IO$LatQIMedQIProp)
# Not worth adding group alongside host
# But definitely worth adding host alongside group

sim_LatQIMedQIProp_HostGroup_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$HostGroup_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_HostGroup_IO)
}

#### Parasite focus ####
##### Parasite species
summary(LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$Parasite_IO$LatQIMedQIProp, LM_IUCN_Species$Fixed$LatQIMedQIProp)
# defo better than nothing

sim_LatQIMedQIProp_Parasite_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_Parasite_IO)
testOverdispersion(sim_LatQIMedQIProp_Parasite_IO)

# Much better with parasite
# The slope doesn't show up

##### Parasite Type
{
  summary(LM_IUCN_Species$Type_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$Type_IO$LatQIMedQIProp, LM_IUCN_Species$Fixed$LatQIMedQIProp)
# better than nothing

sim_LatQIMedQIProp_Type_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Type_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_Type_IO)
# still no slope introduced

##### ParType/Parasite
summary(LM_IUCN_Species$ParType_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Type_IO$LatQIMedQIProp)

anova(LM_IUCN_Species$ParType_IO$LatQIMedQIProp, LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$ParType_IO$LatQIMedQIProp, LM_IUCN_Species$Type_IO$LatQIMedQIProp)
# Same as before, type doesn't add to parasite, but parasite adds to type

sim_LatQIMedQIProp_ParType_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$ParType_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_ParType_IO)
plot(sim_LatQIMedQIProp_Parasite_IO)
plot(sim_LatQIMedQIProp_Type_IO)
# Looks nice, barely different from parasite
}

#### Host and Parasite ####
##### Host + Parasite
summary(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Host_IO$LatQIMedQIProp)

anova(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp, LM_IUCN_Species$Host_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp, LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
# Better than host, better than parasite

sim_LatQIMedQIProp_BothSpecies_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_Parasite_IO)
plot(sim_LatQIMedQIProp_BothSpecies_IO)
testDispersion(sim_LatQIMedQIProp_BothSpecies_IO)
# Doesn't fully get rid of the slope this time!

##### Host:Parasite Type
{
summary(LM_IUCN_Species$CrossSpecies_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
summary(LM_IUCN_Species$Host_IO$LatQIMedQIProp)

anova(LM_IUCN_Species$CrossSpecies_IO$LatQIMedQIProp, LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$CrossSpecies_IO$LatQIMedQIProp, LM_IUCN_Species$Host_IO$LatQIMedQIProp)
anova(LM_IUCN_Species$CrossSpecies_IO$LatQIMedQIProp, LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)

# Better than host, better than parasite species, not better than both

sim_LatQIMedQIProp_CrossSpecies_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$CrossSpecies_IO$LatQIMedQIProp)

plot(sim_LatQIMedQIProp_Host_IO)
plot(sim_LatQIMedQIProp_Parasite_IO)
plot(sim_LatQIMedQIProp_BothSpecies_IO)
plot(sim_LatQIMedQIProp_CrossSpecies_IO)
# Looks okay
# TODO: tidy up models that don't exist anymore
}

#### Random effect conclusions ####
if(plot_outputs){
  # Fixed
  sim_LatQIMedQIProp_Fixed <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed$LatQIMedQIProp)
  # Host
  sim_LatQIMedQIProp_Host_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Host_IO$LatQIMedQIProp)
  # Parasite 
  sim_LatQIMedQIProp_Parasite_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Parasite_IO$LatQIMedQIProp)
  # Both
  sim_LatQIMedQIProp_BothSpecies_IO <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$BothSpecies_IO$LatQIMedQIProp)
  # Doesn't fully get rid of the slope this time!
  

  pdf(here::here("Figures/random effect diagnostics.pdf"), width = 10, height = 7)
  plot(sim_LatQIMedQIProp_Fixed)
  mtext("A, no random terms", side=3, line = 0.5, adj = 0.43)
  plot(sim_LatQIMedQIProp_Host_IO)
  mtext("B, host random intercept", side=3, line = 0.5, adj = 0.43)
  plot(sim_LatQIMedQIProp_Parasite_IO)
  mtext("C, parasite random intercept", side=3, line = 0.5, adj = 0.43)
  plot(sim_LatQIMedQIProp_BothSpecies_IO)
  mtext("D, both random intercept", side=3, line = 0.5, adj = 0.43)
  dev.off()
}
### Host Phylogeny ####
# TODO: source Host Phylo script separately

library(U.PhyloMaker)
library(tidyverse)
library(DHARMa)
# library(phyr)

megatree <- read.tree('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_megatree.tre')
sp.list  <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_sample_species_list.csv', sep=",")
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_genus_list.csv', sep=",")

# using sp.list as example to get formatting right
names(sp.list)

phylo_model_data <- GMPD_Both_Species %>% 
  mutate(species = case_when(HostCorrectedName == "Cervus canadensis"      ~ "Cervus elaphus",
                             HostCorrectedName == "Lycalopex gymnocercus"  ~ "Pseudalopex gymnocercus",
                             TRUE ~ HostCorrectedName)) 

sp.list <- phylo_model_data %>% 
  mutate(genus = str_split_i(species, " ", 1)) %>% 
  select(species, genus) %>% 
  distinct()


result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result

# Create a distance matrix for phylo in residuals
phyloMat <- cophenetic.phylo(result) %>% as.data.frame.array() 

# rownames(phyloMat)[rownames(phyloMat) == "Cervus elaphus"] <- "Cervus canadensis"
# rownames(phyloMat)[rownames(phyloMat) == "Pseudalopex_gymnocercus"] <- "Lycalopex_gymnocercus"

# Matrix_Data <- phylo_model_data %>% 
#   dplyr::select(HostCorrectedName) %>% 
#   mutate(species = str_replace_all(HostCorrectedName, " ", "_")) %>% 
#   dplyr::select(species) %>% 
#   left_join(rownames_to_column(phyloMat), by = join_by(species == rowname)) 
# 
# distMat <- as.matrix(Matrix_Data[, Matrix_Data$species],
#                      dimnames = Matrix_Data$species)

phylo_models <- list()
phylo_models$Host_IO$LatQIMedQIProp        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistanceSquared, 2)) + (1|species), data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$Parasite_IO$LatQIMedQIProp    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistanceSquared, 2)) + (1|ParasiteCorrectedName), data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$BothSpecies_IO$LatQIMedQIProp <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistanceSquared, 2)) + (1|species) + (1|ParasiteCorrectedName), data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)

phylo_models$Host_IO$sim <- DHARMa::simulateResiduals(fittedModel = phylo_models$Host_IO$LatQIMedQIProp)
phylo_models$Host_IO$hostrecalc <- DHARMa::recalculateResiduals(phylo_models$Host_IO$sim, group = phylo_model_data$species)
if(plot_outputs){
  pdf(here::here("Figures/phylogeny host.pdf"), width = 10, height = 7)
  plot(phylo_models$Host_IO$hostrecalc)
  mtext("Host random intercepts", side=3, line = 0.5, adj = 0.43)
  dev.off()
}
DHARMa::testSpatialAutocorrelation(simulationOutput = phylo_models$Host_IO$hostrecalc, distMat = phyloMat)
# definitely no signal here, p-value = 0.9
# but there is heteroskedasticity

phylo_models$Parasite_IO$sim <- DHARMa::simulateResiduals(fittedModel = phylo_models$Parasite_IO$LatQIMedQIProp)
phylo_models$Parasite_IO$hostrecalc <- DHARMa::recalculateResiduals(phylo_models$Parasite_IO$sim, group = phylo_model_data$species)
if(plot_outputs){
  pdf(here::here("Figures/phylogeny parasite.pdf"), width = 10, height = 7)
  plot(phylo_models$Parasite_IO$hostrecalc)
  mtext("Parasite random intercepts", side=3, line = 0.5, adj = 0.43)
  dev.off()
}
DHARMa::testSpatialAutocorrelation(simulationOutput = phylo_models$Parasite_IO$hostrecalc, distMat = phyloMat)
# although there's no phylogenetic signal, parasites seem to be introducing something into the residuals vs predicted
# might be because different hosts have different sets of parasites
# either way, there is no reason to include phylogeny and the diagnostics without phylogeny look fine too

phylo_models$BothSpecies_IO$sim <- DHARMa::simulateResiduals(fittedModel = phylo_models$BothSpecies_IO$LatQIMedQIProp)
phylo_models$BothSpecies_IO$hostrecalc <- DHARMa::recalculateResiduals(phylo_models$BothSpecies_IO$sim, group = phylo_model_data$species)
if(plot_outputs){
  pdf(here::here("Figures/phylogeny both.pdf"), width = 10, height = 7)
  plot(phylo_models$BothSpecies_IO$hostrecalc)
  mtext("Both random intercepts", side=3, line = 0.5, adj = 0.43)
  dev.off()
}
DHARMa::testSpatialAutocorrelation(simulationOutput = phylo_models$BothSpecies_IO$hostrecalc, distMat = phyloMat)
# p-value is 0.7 so we accept the null hypothesis that there's no strong phylogenetic signal
# although there appears to be a slope in the plot


# Trying a phylogenetic model
# phylo_models$Phylo_IO   <- try(phyr::pglmm(formula = cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + 
#                                         poly(CorrectedDistPropSquared, 2)) + (1|species),
#                                       data = phylo_model_data, family = "binomial",
#                                       cov_ranef = list(species = result)), silent = FALSE)

# this model has too many parameters to run but we've confirmed above that there's no worry about host phylogeny

# Plots ########################################################################

### Fig 1, map ####
# map

if(plot_outputs){
  # source data
  IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))
  Projection_String <- sf::st_crs(IUCN_Native_Data)
  
  # Map
  sf::sf_use_s2(FALSE)
  tmap_mode("plot")
  data("World")
  
  IUCN_Map_Data <- IUCN_Native_Data %>% 
    mutate(Group = case_when(order_ == "CARNIVORA" ~ "Carnivores",
                             TRUE ~ "Ungulates")) 
  
  GMPD_Map_Data <- sf::st_as_sf(GMPD_Both_Species,
                                coords = c("Longitude", "Latitude"),
                                crs = Projection_String) 
  
  cairo_pdf(here::here("Figures/map.pdf"), width = 10, height = 5)
  tm_shape(World) + tm_fill() +
    tm_layout(bg.color = "white", fontfamily = "Outfit",
              legend.position = c("left", "center")) +
    tm_shape(IUCN_Map_Data) + tm_fill(col = "Group", palette = c("#8E8DBE", "#81F495"), alpha = 0.3) +
    tm_shape(GMPD_Map_Data) + tm_dots(col = "SampleSize", 
                                      palette = c("#222E50", 
                                                  "#D1495B",
                                                  "#EDAE49"
                                      ),
                                      style = "log10")
  
  # tm_add_legend(type = c("symbol"),
  #               col = c("#B99878", "#A1D6E2"),
  #               border.col = NA,
  #               labels = c("Carnivores", "Ungulates"),)
  dev.off()
}


## Fig 2, lat ####
# main plot
if(plot_outputs){

  Latitude_main <- plot_latitude(model_list = LM_IUCN_Species, 
                                 model_data = GMPD_Both_Species) +
    scale_x_continuous(breaks = c(0,25,40,45,55,90), labels = c("0","25","40","45","55","90")) +
    theme(axis.text.x = element_text(color = c("#656565", "black", "black", "#656565", "black", "#656565")),
          axis.ticks.x = element_line(color = c("#656565", "black", "black", "#656565", "black", "#656565"),
                                      linewidth = c(.5,1,1,.5,1,.5)))

  ggsave(here::here("Figures/latitude main.pdf"), Latitude_main, width = 10, height = 7, 
         device = cairo_pdf)
}


## Fig 3, main ####
if(plot_outputs){
  MedianProp_25  <- plot_medianprop(model_list = LM_IUCN_Species, 
                                    model_data = GMPD_Both_Species, 
                                    fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "A, latitude = 25",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_40 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_Both_Species, 
                                   fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "C, latitude = 40",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_55 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_Both_Species, 
                                   fixed_lat = 55) +
    # theme(axis.text.x = element_blank()) +
    labs(x = "Range position", y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E, latitude = 55",
             family = "Outfit", size = 5, hjust = 0)
  
  
  DistProp_25  <- plot_distprop(model_list = LM_IUCN_Species, 
                                model_data = GMPD_Both_Species, 
                                fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_40 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_data = GMPD_Both_Species, 
                               fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_55 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_data = GMPD_Both_Species, 
                               fixed_lat = 55) +
    # theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = "Niche position", y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "F",
             family = "Outfit", size = 5, hjust = 0)
  
  gg_y_axis <- cowplot::get_plot_component(ggplot() + 
                                             theme(text = element_text(family = "Outfit", size = 15)) + 
                                             labs(y = "Parasite prevalence"), 
                                           "ylab-l")
  gg_legend <- cowplot::get_plot_component(Latitude_main, 'guide-box-bottom', return_all = TRUE)                    


  range_niche_plots <- MedianProp_25 + DistProp_25 +
    MedianProp_40 + DistProp_40 +
    MedianProp_55 + DistProp_55 +
    gg_y_axis +
    gg_legend +
    plot_layout(design = "
              #AB
              GCD
              #EF
              #HH
              ",
              widths = c(1, 20, 20),
              heights = c(20,20,20,1)) 
  
  ggsave(here::here("Figures/range niche position.pdf"), range_niche_plots, width = 10, height = 10, 
         device = cairo_pdf)
}

## Supp fig 1, asym ####
# Range position asymmetry
if(plot_outputs) {
  MedianProp_asym_40  <- plot_medianprop_asym(model_list = LM_GBIF_Species, 
                                              model_data = GMPD_Both_Species, 
                                              fixed_lat = 40) 
  ggsave(here::here("Figures/range asymmetry.pdf"), MedianProp_asym_40, width = 10, height = 7, 
         device = cairo_pdf)
  
}

## Supp fig 2, gbif ####
# GBIF vs IUCN
if(plot_outputs){
  Latitude_IUCN <- plot_latitude(model_list = LM_IUCN_Species, 
                                 model_data = GMPD_Both_Species) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "A",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  Latitude_GBIF  <- plot_latitude(model_list = LM_GBIF_Species, 
                                  model_data = GMPD_Both_Species) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  
  MedianProp_IUCN <- plot_medianprop(model_list = LM_IUCN_Species, 
                                     model_data = GMPD_Both_Species, 
                                     fixed_lat = 40) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "C",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_GBIF  <- plot_medianprop(model_list = LM_GBIF_Species, 
                                      model_data = GMPD_Both_Species, 
                                      fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0) 
  
  DistProp_IUCN <- plot_distprop(model_list = LM_IUCN_Species, 
                                 model_data = GMPD_Both_Species, 
                                 fixed_lat = 40) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_GBIF  <- plot_distprop(model_list = LM_GBIF_Species, 
                                  model_data = GMPD_Both_Species, 
                                  fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "F",
             family = "Outfit", size = 5, hjust = 0)
  
  gg_lat_axis <- cowplot::get_plot_component(ggplot() + 
                                               theme(text = element_text(family = "Outfit", size = 15)) + 
                                               labs(x = "Latitude"), 
                                             "xlab-b")
  gg_range_axis <- cowplot::get_plot_component(ggplot() + 
                                                 theme(text = element_text(family = "Outfit", size = 15)) + 
                                                 labs(x = "Range position"), 
                                               "xlab-b")
  gg_niche_axis <- cowplot::get_plot_component(ggplot() + 
                                                 theme(text = element_text(family = "Outfit", size = 15)) + 
                                                 labs(x = "Niche position"), 
                                               "xlab-b")
  
  gg_y_axis <- cowplot::get_plot_component(ggplot() + 
                                             theme(text = element_text(family = "Outfit", size = 15)) + 
                                             labs(y = "Parasite prevalence"), 
                                           "ylab-l")
  
  gg_legend <- cowplot::get_plot_component(Latitude_main, 'guide-box-bottom', return_all = TRUE)                    
  
  
  gbif_iucn_plots <- Latitude_IUCN +  Latitude_GBIF +
    gg_lat_axis +
    MedianProp_IUCN + MedianProp_GBIF +
    gg_range_axis +
    DistProp_IUCN + DistProp_GBIF +
    gg_niche_axis +
    gg_y_axis +
    gg_legend +
    plot_layout(
      design = "
               #AB
               #CC
               JDE
               #FF
               #GH
               #II
               #KK
               ",
      heights = c(20,2,20,2,20,2,2),
    widths = c(2,40,40))


ggsave(here::here("Figures/iucn gbif.pdf"), gbif_iucn_plots, width = 10, height = 13, 
       device = cairo_pdf)
}

## Supp fig 3, raster ####
# Raster resolution
if(plot_outputs){
  Latitude_5  <- plot_latitude(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIProp_5", 
                               model_data = GMPD_Both_Species,
                               raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "A",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  Latitude_10 <- plot_latitude(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  Latitude_20 <- plot_latitude(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIProp_20", 
                               model_data = GMPD_Both_Species,
                               raster_res = "_20") +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "C",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  
  MedianProp_5  <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_name = "LatQIMedQIProp_5", 
                                   model_data = GMPD_Both_Species, 
                                   fixed_lat = 40,
                                   raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_10 <- plot_medianprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_20 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_name = "LatQIMedQIProp_20", 
                                   model_data = GMPD_Both_Species, 
                                   fixed_lat = 40,
                                   raster_res = "_20") +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "F",
             family = "Outfit", size = 5, hjust = 0)
  
  DistProp_5  <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIProp_5", 
                               model_data = GMPD_Both_Species, 
                               fixed_lat = 40,
                               raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "G",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_10 <- plot_distprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "H",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_20 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIProp_20", 
                               model_data = GMPD_Both_Species, 
                               fixed_lat = 40,
                               raster_res = "_20") +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "I",
             family = "Outfit", size = 5, hjust = 0)
  
  gg_lat_axis <- cowplot::get_plot_component(ggplot() + 
                                               theme(text = element_text(family = "Outfit", size = 15)) + 
                                               labs(x = "Latitude"), 
                                             "xlab-b")
  gg_range_axis <- cowplot::get_plot_component(ggplot() + 
                                                 theme(text = element_text(family = "Outfit", size = 15)) + 
                                                 labs(x = "Range position"), 
                                               "xlab-b")
  gg_niche_axis <- cowplot::get_plot_component(ggplot() + 
                                                 theme(text = element_text(family = "Outfit", size = 15)) + 
                                                 labs(x = "Niche position"), 
                                               "xlab-b")
  
  gg_y_axis <- cowplot::get_plot_component(ggplot() + 
                                             theme(text = element_text(family = "Outfit", size = 15)) + 
                                             labs(y = "Parasite prevalence"), 
                                           "ylab-l")
  
  gg_legend <- cowplot::get_plot_component(Latitude_main, 'guide-box-bottom', return_all = TRUE)
  
  
  raster_res_plots <- Latitude_5 +  Latitude_10 + Latitude_20 +
    gg_lat_axis +
    MedianProp_5 + MedianProp_10 + MedianProp_20 +
    gg_range_axis +
    DistProp_5 + DistProp_10 + DistProp_20 +
    gg_niche_axis +
    gg_y_axis +
    gg_legend +
    plot_layout(
      design = "
               #ABC
               ##D#
               MEFG
               ##H#
               #IJK
               ##L#
               #NNN
               ",
      heights = c(20,2,20,2,20,2),
      widths = c(2,40,40,40))
  
  
  ggsave(here::here("Figures/raster res.pdf"), raster_res_plots, width = 10, height = 9, 
         device = cairo_pdf)
}
