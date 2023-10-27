# Analysis and plotting of final dataset
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ###########################################################

library(here)               #
library(tidyverse)          #
library(corrplot)

source(here::here("Functions.R"))

GMPD_Climate_Data <- read.csv(here::here("Data/Data back ups/GMPD_Climate_Data_04.csv"), 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

# temp until I rerun earlier scripts
GMPD_Parasite_Data <- read.csv(here::here("Data/GMPD_datafiles/GMPD_main.csv"), 
                          header = TRUE, 
                          stringsAsFactors = FALSE) %>% 
  dplyr::select(ParasiteCorrectedName, ParType, ParPhylum, ParClass) %>% 
  dplyr::distinct() %>% 
  dplyr::right_join(GMPD_Climate_Data, by = join_by("ParasiteCorrectedName"))

# Prep and data exploration ####################################################
# Hosts ####
# Check whether there are any host groups which should be removed 

pdf(here::here("Figures/Host Group.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = Group, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Order.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostOrder, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Family.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = HostFamily, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Host Species.pdf"), width = 10, height = 25)

basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = Prevalence) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = Latitude) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
basic_barplot(dat = GMPD_Climate_Data, xvar = HostCorrectedName, yvar = AbsLatitude) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

dev.off()

## Group: uninformative considering we will include phylogeny
## Order: Doesn't give much more than Group. Arguably could remove Perissodactyla because of few rows, 
# but it's geographically restricted to regions where data is already lacking and will be accounted for in phylogeny
## Family: Too many factor levels to make much sense of it, and better to account for it in phylogeny

## Parasites ####
# Check whether there are any parasite groups which should be removed

pdf(here::here("Figures/Parasite Class.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParClass, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Parasite Type.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParType, yvar = AbsLatitude)

dev.off()

pdf(here::here("Figures/Parasite Phylum.pdf"), width = 8, height = 6)

basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = Prevalence)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = Latitude)
basic_barplot(dat = GMPD_Climate_Data, xvar = ParPhylum, yvar = AbsLatitude)

dev.off()

## Class: too many to make sense of
## Type: Fungi and prions are a bit sketchy. Both have relatively few rows and are geographically restricted
# Fungi are all in Finland, UK, and Czechia/Slovakia. Prions are all in USA in one area
# removing these two groups from all analyses
## Phylum: Ascomycota and prions are removed above. 
# Sarcomastigophora is a bit more acceptable because it's within Protozoa, can remove if I use Phylum

GMPD_Analysis_Data <- GMPD_Parasite_Data %>% 
  dplyr::filter(!ParType %in% c("Prion", "Fungus"))

## Correlations ####
GMPD_Analysis_Data %>% 
  dplyr::select(Prevalence, HostsSampled, Latitude, MedianProp, 
                CorrectedDensity, CorrectedDistance, CorrectedAngle,
                Axis1, Axis2, 
                wc2.1_10m_bio_5, wc2.1_10m_bio_6, wc2.1_10m_bio_12) %>% 
  cor() %>% 
  corrplot::corrplot()

# https://doi.org/10.1111/j.1600-0587.2012.07348.x
# None of the variables I plan on modelling 
# (Latitude, MedianProp, Axis1, Axis2, CorrectedDensity, CorrectedDistance, CorrectedAngle)
# have particularly high collinearity, though I will keep an eye on
# Latitude + Axis1, and CorrectedDistance + CorrectedDensity

# Analysis  ####################################################################
## Prepping data ####
GMPD_IUCN_Species <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "species") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                          TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                CorrectedDistanceLogged = base::scale(log(CorrectedDistance + 1)),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                HostSubgroup            = paste(HostCorrectedName, subgroup))

## Prepping data ####
GMPD_IUCN_Subgroup <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "iucn",
                RangeTaxonLvl == "subgroup") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                HostSubgroup            = paste(HostCorrectedName, subgroup))

## Prepping data ####
GMPD_GBIF_Species <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "gbif",
                RangeTaxonLvl == "species") %>%
  dplyr::mutate(LatitudeScaled          = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled        = base::scale(MedianProp),
                MedianPropSquared       = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled     = base::scale(MedianPropSquared),
                
                Axis1Scaled             = base::scale(Axis1),
                Axis2Scaled             = base::scale(Axis2),
                
                CorrectedDensityScaled  = base::scale(CorrectedDensity),
                CorrectedDistanceScaled = base::scale(CorrectedDistance),
                CorrectedAngleScaled    = base::scale(CorrectedAngle),
                
                ParasiteDetected        = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected      = HostsSampled - ParasiteDetected,
                HostSubgroup            = paste(HostCorrectedName, subgroup))

## Prepping data ####
GMPD_GBIF_Subgroup <- GMPD_Analysis_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod   == "gbif",
                RangeTaxonLvl == "subgroup") %>%
  dplyr::mutate(LatitudeScaled            = base::scale(abs(Latitude)),
                # EquatorwardsPropScaled  = base::scale(EquatorwardsProp),
                MedianPropScaled          = base::scale(MedianProp),
                MedianPropSquared         = dplyr::case_when(!AboveMedn ~ MedianProp * -1,
                                                           TRUE       ~ MedianProp),
                MedianPropSquScaled       = base::scale(MedianPropSquared),
                
                Axis1Scaled               = base::scale(Axis1),
                Axis2Scaled               = base::scale(Axis2),
                
                CorrectedDensityScaled    = base::scale(CorrectedDensity),
                CorrectedDistanceScaled   = base::scale(CorrectedDistance),
                CorrectedAngleScaled      = base::scale(CorrectedAngle),

                ParasiteDetected          = as.integer(round(HostsSampled * Prevalence, 0)),
                ParasiteUndetected        = HostsSampled - ParasiteDetected,
                HostSubgroup              = paste(HostCorrectedName, subgroup))

## Frequentist ####
source(here::here("05.1Frequentist analysis.R"))

# TODO: write a summary for frequentist analysis so far
## Explain why I chose to look at Type instead of Phylum
# Accept Latitude, quadratic MedianProp, and interaction, and AboveMedn
# Write notes on the ones I don't accept (climatic terms)
# Check if there are differences between corrected and uncorrected values
# Look at differences in these effects when I add Host
# Same again for Parasite Type
# Same again for both Type and Host, but say where it ran out of beans


### Conclusions ####

## Bayesian ####
source(here::here("05.3Bayesian analysis.R"))
### Conclusions ####

# Plots ########################################################################

