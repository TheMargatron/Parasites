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
library(MuMIn)

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

GMPD_Climate_Data %>% 
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
test_data <- GMPD_Climate_Data %>% 
  filter(RestrAll,
         CleanAll) %>% 
  mutate(LatitudeScaled    = as.vector(scale(abs(Latitude))),
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
# IUCN filter
GMPD_IUCN_Species <- GMPD_Climate_Data %>%
  dplyr::filter(RestrAll) %>%
  dplyr::mutate(LatitudeScaled            = as.vector(base::scale(abs(Latitude))),
                MedianPropScaled_iucn     = as.vector(base::scale(MedianProp_iucn)),
                MedianPropSquared_iucn    = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                                                TRUE       ~ MedianProp_iucn),
                MedianPropSquScaled_iucn  = as.vector(base::scale(MedianPropSquared_iucn)),
                
                Axis1Scaled               = as.vector(base::scale(Axis1)),
                Axis2Scaled               = as.vector(base::scale(Axis2)),
                
                CorrectedDensityScaled    = as.vector(base::scale(CorrectedDensity)),
                CorrectedDistanceScaled   = as.vector(base::scale(CorrectedDistance)),
                CorrectedAngleScaled      = as.vector(base::scale(CorrectedAngle)),
                
                CorrectedDensityScaled_5  = as.vector(base::scale(CorrectedDensity_5)),
                CorrectedDistanceScaled_5 = as.vector(base::scale(CorrectedDistance_5)),
                CorrectedAngleScaled_5    = as.vector(base::scale(CorrectedAngle_5)),
                
                CorrectedDensityScaled_20  = as.vector(base::scale(CorrectedDensity_20)),
                CorrectedDistanceScaled_20 = as.vector(base::scale(CorrectedDistance_20)),
                CorrectedAngleScaled_20    = as.vector(base::scale(CorrectedAngle_20)),
                
                ParasiteDetected          = as.integer(round(SampleSize * Prevalence, 0)),
                ParasiteUndetected        = SampleSize - ParasiteDetected,
                
                # HostSubgroup            = paste(HostCorrectedName, subgroup, sep = "_"),
                HostParasite              = paste(HostCorrectedName, ParasiteCorrectedName, sep = "_"),
                
                RowID                     = row.names(.)) %>% 
  
  # Proportional distance
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>%
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp),
                CorrectedDistPropSquScaled = as.vector(base::scale(CorrectedDistPropSquared))) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistProp_5, CorrectedAngle_5))) %>%
  dplyr::mutate(TopHalf_5 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_5 = case_when(!TopHalf ~ CorrectedDistProp_5 * -1,
                                                       TRUE ~ CorrectedDistProp_5),
                CorrectedDistPropSquScaled_5 = as.vector(base::scale(CorrectedDistPropSquared_5))) %>%
  dplyr::arrange(desc(pmax(CorrectedDistProp_20, CorrectedAngle_20))) %>%
  
  dplyr::mutate(TopHalf_20 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_20 = case_when(!TopHalf ~ CorrectedDistProp_20 * -1,
                                                        TRUE ~ CorrectedDistProp_20),
                CorrectedDistPropSquScaled_20 = as.vector(base::scale(CorrectedDistPropSquared_20)))


GMPD_IUCN_Species %>% 
  dplyr::select(Prevalence, SampleSize, Latitude, 
                MedianProp_iucn, MedianPropSquared_iucn,
                CorrectedDistProp, CorrectedDistPropSquared) %>% 
  mutate(abslatitude = abs(Latitude)) %>% 
  cor() %>% 
  corrplot::corrplot.mixed(diag = "n", tl.pos = "lt")

# GBIF filter
GMPD_GBIF_Species <- GMPD_Climate_Data %>%
  dplyr::filter(CleanAll) %>%
  dplyr::mutate(LatitudeScaled          = as.vector(base::scale(abs(Latitude))),
                MedianPropScaled_gbif   = as.vector(base::scale(MedianProp_gbif)),
                MedianPropSquared_gbif  = dplyr::case_when(!AboveMedn_gbif ~ MedianProp_gbif * -1,
                                                           TRUE       ~ MedianProp_gbif),
                MedianPropSquScaled_gbif = as.vector(base::scale(MedianPropSquared_gbif)),
                
                Axis1Scaled             = as.vector(base::scale(Axis1)),
                Axis2Scaled             = as.vector(base::scale(Axis2)),
                
                CorrectedDensityScaled  = as.vector(base::scale(CorrectedDensity)),
                CorrectedDistanceScaled = as.vector(base::scale(CorrectedDistance)),
                CorrectedAngleScaled    = as.vector(base::scale(CorrectedAngle)),
                
                CorrectedDensityScaled_5  = as.vector(base::scale(CorrectedDensity_5)),
                CorrectedDistanceScaled_5 = as.vector(base::scale(CorrectedDistance_5)),
                CorrectedAngleScaled_5    = as.vector(base::scale(CorrectedAngle_5)),
                
                CorrectedDensityScaled_20  = as.vector(base::scale(CorrectedDensity_20)),
                CorrectedDistanceScaled_20 = as.vector(base::scale(CorrectedDistance_20)),
                CorrectedAngleScaled_20    = as.vector(base::scale(CorrectedAngle_20)),
                
                ParasiteDetected        = as.integer(round(SampleSize * Prevalence, 0)),
                ParasiteUndetected      = SampleSize - ParasiteDetected,
                
                # HostSubgroup            = paste(HostCorrectedName, subgroup, sep = "_"),
                HostParasite            = paste(HostCorrectedName, ParasiteCorrectedName, sep = "_"),
                
                RowID                   = row.names(.)) %>% 
  
  # Proportional distance
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>%
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp),
                CorrectedDistPropSquScaled = as.vector(base::scale(CorrectedDistPropSquared))) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistProp_5, CorrectedAngle_5))) %>%
  dplyr::mutate(TopHalf_5 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_5 = case_when(!TopHalf ~ CorrectedDistProp_5 * -1,
                                                       TRUE ~ CorrectedDistProp_5),
                CorrectedDistPropSquScaled_5 = as.vector(base::scale(CorrectedDistPropSquared_5))) %>%
  dplyr::arrange(desc(pmax(CorrectedDistProp_20, CorrectedAngle_20))) %>%
  
  dplyr::mutate(TopHalf_20 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_20 = case_when(!TopHalf ~ CorrectedDistProp_20 * -1,
                                                        TRUE ~ CorrectedDistProp_20),
                CorrectedDistPropSquScaled_20 = as.vector(base::scale(CorrectedDistPropSquared_20)))


GMPD_GBIF_Species %>% 
  dplyr::select(Prevalence, SampleSize, Latitude, 
                MedianProp_gbif, MedianPropSquared_gbif,
                CorrectedDistProp, CorrectedDistPropSquared) %>% 
  mutate(abslatitude = abs(Latitude)) %>% 
  cor() %>% 
  corrplot::corrplot.mixed(diag = "n", tl.pos = "lt")

# both filters
GMPD_Both_Species <- GMPD_Climate_Data %>%
  dplyr::filter(RestrAll,
                CleanAll) %>%
  dplyr::mutate(LatitudeScaled                = as.vector(base::scale(abs(Latitude))),
                MedianPropScaled_iucn         = as.vector(base::scale(MedianProp_iucn)),
                MedianPropSquared_iucn        = dplyr::case_when(!AboveMedn_iucn ~ MedianProp_iucn * -1,
                                                                TRUE       ~ MedianProp_iucn),
                MedianPropSquScaled_iucn      = as.vector(base::scale(MedianPropSquared_iucn)),
                
                MedianPropScaled_gbif         = as.vector(base::scale(MedianProp_gbif)),
                MedianPropSquared_gbif        = dplyr::case_when(!AboveMedn_gbif ~ MedianProp_gbif * -1,
                                                           TRUE       ~ MedianProp_gbif),
                MedianPropSquScaled_gbif      = as.vector(base::scale(MedianPropSquared_gbif)),
                
                Axis1Scaled             = as.vector(base::scale(Axis1)),
                Axis2Scaled             = as.vector(base::scale(Axis2)),
                
                CorrectedDensityScaled  = as.vector(base::scale(CorrectedDensity)),
                CorrectedDistanceScaled = as.vector(base::scale(CorrectedDistance)),
                CorrectedAngleScaled    = as.vector(base::scale(CorrectedAngle)),
                
                CorrectedDensityScaled_5  = as.vector(base::scale(CorrectedDensity_5)),
                CorrectedDistanceScaled_5 = as.vector(base::scale(CorrectedDistance_5)),
                CorrectedAngleScaled_5    = as.vector(base::scale(CorrectedAngle_5)),
                
                CorrectedDensityScaled_20  = as.vector(base::scale(CorrectedDensity_20)),
                CorrectedDistanceScaled_20 = as.vector(base::scale(CorrectedDistance_20)),
                CorrectedAngleScaled_20    = as.vector(base::scale(CorrectedAngle_20)),
                
                ParasiteDetected        = as.integer(round(SampleSize * Prevalence, 0)),
                ParasiteUndetected      = SampleSize - ParasiteDetected,
                
                # HostSubgroup            = paste(HostCorrectedName, subgroup, sep = "_"),
                HostParasite            = paste(HostCorrectedName, ParasiteCorrectedName, sep = "_"),
                
                RowID                   = row.names(.)) %>% 
  
  # Proportional distance
  dplyr::arrange(desc(pmax(CorrectedDistProp, CorrectedAngle))) %>%
  dplyr::mutate(TopHalf = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared = case_when(!TopHalf ~ CorrectedDistProp * -1,
                                                     TRUE ~ CorrectedDistProp),
                CorrectedDistPropSquScaled = as.vector(base::scale(CorrectedDistPropSquared))) %>% 
  
  dplyr::arrange(desc(pmax(CorrectedDistProp_5, CorrectedAngle_5))) %>%
  dplyr::mutate(TopHalf_5 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_5 = case_when(!TopHalf ~ CorrectedDistProp_5 * -1,
                                                       TRUE ~ CorrectedDistProp_5),
                CorrectedDistPropSquScaled_5 = as.vector(base::scale(CorrectedDistPropSquared_5))) %>%
  dplyr::arrange(desc(pmax(CorrectedDistProp_20, CorrectedAngle_20))) %>%
  
  dplyr::mutate(TopHalf_20 = rep(c(TRUE, FALSE), length.out = nrow(.)),
                CorrectedDistPropSquared_20 = case_when(!TopHalf ~ CorrectedDistProp_20 * -1,
                                                        TRUE ~ CorrectedDistProp_20),
                CorrectedDistPropSquScaled_20 = as.vector(base::scale(CorrectedDistPropSquared_20)))


GMPD_Both_Species %>% 
  dplyr::select(Prevalence, SampleSize, Latitude, 
                MedianProp_iucn, MedianPropSquared_iucn,
                MedianProp_gbif, MedianPropSquared_gbif,
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
LM_IUCN_Species_B  <- list()
LM_GBIF_Species_B  <- list()

## Model sub lists
# Fixed         = Fixed effects only

# Host          
## Host         = Host species random effect
## Group        = Host group as a random effect
## HostGroup    = Host species nested within host group as random effects

# Parasite      
## Parasite     = Parasite species as random effect
## Type         = Parasite type as random effect
## ParType      = Parasite species nested within Parasite type as random effects

# Both
## Both  = Host species and parasite species as random effects
## Cross = Host species crossed with parasite species as random effects

if(rerun_models){
  LM_IUCN_Species   <- run_models(model_data = GMPD_IUCN_Species, method = "iucn")
  LM_GBIF_Species   <- run_models(model_data = GMPD_GBIF_Species, method = "gbif")
  LM_IUCN_Species_B <- run_models(model_data = GMPD_Both_Species, method = "iucn")
  LM_GBIF_Species_B <- run_models(model_data = GMPD_Both_Species, method = "gbif")
  
  saveRDS(LM_IUCN_Species,   here::here("Data/Model back ups/LM_IUCN_Species_05.rds"))
  saveRDS(LM_GBIF_Species,   here::here("Data/Model back ups/LM_GBIF_Species_05.rds"))
  saveRDS(LM_IUCN_Species_B, here::here("Data/Model back ups/LM_IUCN_Species_B_05.rds"))
  saveRDS(LM_GBIF_Species_B, here::here("Data/Model back ups/LM_GBIF_Species_B_05.rds"))
} else{
  
  LM_IUCN_Species   <- readRDS(here::here("Data/Model back ups/LM_IUCN_Species_05.rds"))
  LM_GBIF_Species   <- readRDS(here::here("Data/Model back ups/LM_GBIF_Species_05.rds"))
  LM_IUCN_Species_B <- readRDS(here::here("Data/Model back ups/LM_IUCN_Species_B_05.rds"))
  LM_GBIF_Species_B <- readRDS(here::here("Data/Model back ups/LM_GBIF_Species_B_05.rds"))
}

### Random effects quickly ####
null_random_anova <- list()

# Doing a quick comparison of hosts vs parasites
null_random_anova$Host_Fixed <- anova(LM_IUCN_Species$Host$Null, LM_IUCN_Species$Fixed$Null)
# host better than nothing

null_random_anova$Group_Fixed <- anova(LM_IUCN_Species$Group$Null, LM_IUCN_Species$Fixed$Null)
null_random_anova$HostGroup_Host <- anova(LM_IUCN_Species$HostGroup$Null, LM_IUCN_Species$Host$Null)
null_random_anova$HostGroup_Group <- anova(LM_IUCN_Species$HostGroup$Null, LM_IUCN_Species$Group$Null)
# group better than nothing, but doesn't add anything on top of host
# host adds things on top of group

null_random_anova$Parasite_Fixed <- anova(LM_IUCN_Species$Parasite$Null, LM_IUCN_Species$Fixed$Null)
# parasite better than nothing

null_random_anova$Type_Fixed <- anova(LM_IUCN_Species$Type$Null, LM_IUCN_Species$Fixed$Null)
null_random_anova$ParType_Parasite <- anova(LM_IUCN_Species$ParType$Null, LM_IUCN_Species$Parasite$Null)
null_random_anova$ParType_Type <- anova(LM_IUCN_Species$ParType$Null, LM_IUCN_Species$Type$Null)
# Type better than nothing, and adds a little on top of parasite
# parasite adds things on top of type

null_random_anova$Both_Parasite <- anova(LM_IUCN_Species$Both$Null, LM_IUCN_Species$Parasite$Null)
null_random_anova$Both_Host <- anova(LM_IUCN_Species$Both$Null, LM_IUCN_Species$Host$Null)
# Host adds more on top of parasite, parasite adds more on top of host

null_random_anova$Cross_Both <- anova(LM_IUCN_Species$Cross$Null, LM_IUCN_Species$Both$Null)
# crossing host and parasite doesn't improve? AIC is lower, but P is 1
# seems like they aren't nested from the way the model specifies them

null_random_effects <- c("Host",      "Fixed",
                         "Group",     "Fixed",
                         "HostGroup", "Host",
                         "HostGroup", "Group",
                         "Parasite",  "Fixed",
                         "Type",      "Fixed",
                         "ParType",   "Parasite",
                         "ParType",   "Type",
                         "Both",      "Parasite",
                         "Both",      "Host",
                         "Cross",     "Both") %>% 
  matrix(ncol = 2, byrow = TRUE) %>% 
  as.data.frame() %>% 
  rename("model1" = "V1", "model2" = "V2") %>% 
  mutate(paste_names = paste(model1, model2, sep = "_")) %>% 
  group_by(paste_names) %>% 
  mutate(deltaDF = null_random_anova[[paste_names]]$Df[[2]],
         m1 = str_split_i(rownames(null_random_anova[[paste_names]])[[1]], "\\$", 2),
         m2 = str_split_i(rownames(null_random_anova[[paste_names]])[[2]], "\\$", 2),
         
         npar1 = case_when(model1 == m1 ~ null_random_anova[[paste_names]]$npar[[1]],
                           model1 == m2 ~ null_random_anova[[paste_names]]$npar[[2]]),
         npar2 = case_when(model2 == m2 ~ null_random_anova[[paste_names]]$npar[[2]],
                           model2 == m1 ~ null_random_anova[[paste_names]]$npar[[1]]),
         
         loglike1 = case_when(model1 == m1 ~ null_random_anova[[paste_names]]$logLik[[1]],
                              model1 == m2 ~ null_random_anova[[paste_names]]$logLik[[2]]),
         loglike2 = case_when(model2 == m2 ~ null_random_anova[[paste_names]]$logLik[[2]],
                              model2 == m1 ~ null_random_anova[[paste_names]]$logLik[[1]]),
         
         Chisq = null_random_anova[[paste_names]]$Chisq[[2]],
         P = null_random_anova[[paste_names]]$`Pr(>Chisq)`[[2]],
         P = case_when(P < 0.001 ~ "< 0.001",
                       TRUE ~ as.character(round(P, 2)))) %>% 
  ungroup() %>% 
  select(-paste_names, -m1, -m2)

write.csv(null_random_effects, here::here("Tables", "null_random_effects.csv"), 
          row.names = FALSE, na = "")

null_random_stdev <- null_random_effects %>% 
  select(model1) %>%
  distinct() %>% 
  mutate(model1 = factor(model1, levels = unique(model1))) %>%
  group_by(model1) %>% 
  summarise(term = as.data.frame(VarCorr(LM_IUCN_Species[[as.character(model1)]]$Null))$grp,
            sdcor_null = as.data.frame(VarCorr(LM_IUCN_Species[[as.character(model1)]]$Null))$sdcor) %>% 
  arrange(nchar(term), .by_group = TRUE)
  

# will compare fixed effects with Both
### Fixed effects ####

#### Latitude
# ______________________________________________________________________________
# starting with Latitude because that's what I'm most interested in
# Latitude is accepted
compare_models(LM_IUCN_Species$Cross$Lat, 
               LM_IUCN_Species$Cross$Null,
               "IUCN data, Latitude vs null")

compare_models(LM_GBIF_Species$Cross$Lat, 
               LM_GBIF_Species$Cross$Null,
               "GBIF data, Latitude vs null")

#### MedianProp
# ______________________________________________________________________________
# Next thing we're interested in is position within geographic range
# MedianProp is accepted in iucn but not in gbif
# Either way we want need to check it as a quadratic

compare_models(LM_IUCN_Species$Cross$LatMed,
               LM_IUCN_Species$Cross$Lat,
               "IUCN data, MedianProp added to Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatMed,
               LM_GBIF_Species$Cross$Lat,
               "GBIF data, MedianProp added to Latitude, random both species")

#### Quadratic MedianProp
# ______________________________________________________________________________
# Probably better to model it as a quadratic
# Better as a quadratic, and included in both iucn and gbif
# Got the same thing of it switching signs, but same conclusion as previous

compare_models(LM_IUCN_Species$Cross$LatQMed,
               LM_IUCN_Species$Cross$LatMed,
               "IUCN data, quadratic added to MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQMed,
               LM_GBIF_Species$Cross$LatMed,
               "GBIF data, quadratic added to MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQMed,
               LM_GBIF_Species$Cross$Lat,
               "GBIF data, quadratic MedianProp added to Latitude, random both species")

#### Asymmetric Quadratic MedianProp
# ______________________________________________________________________________
# Median prop might have an asymmetric pattern depending on which half of the range it's in
# asymmetry included

compare_models(LM_IUCN_Species$Cross$LatQAMed,
               LM_IUCN_Species$Cross$LatQMed,
               "IUCN data, asymmetry added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQAMed,
               LM_GBIF_Species$Cross$LatQMed,
               "GBIF data, asymmetry added to quadratic MedianProp and Latitude, random both species")

#### Latitude * Quadratic MedianProp
# ______________________________________________________________________________
# Testing out the interaction term 
# There is support for an interaction

compare_models(LM_IUCN_Species$Cross$LatQIMed,
               LM_IUCN_Species$Cross$LatQMed,
               "IUCN data, interaction added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQIMed,
               LM_GBIF_Species$Cross$LatQMed,
               "GBIF data, interaction added to quadratic MedianProp and Latitude, random both species")

#### Latitude * Asymmetric Quadratic MedianProp
# ______________________________________________________________________________
# Interaction with the asymmetry
# failed to converge so have to compare which is more informative from previous models of them independently
# The interaction is more informative

summary(LM_IUCN_Species$Cross$LatQIAMed)
summary(LM_GBIF_Species$Cross$LatQIAMed)
summary(LM_IUCN_Species$Both$LatQIAMed)
summary(LM_GBIF_Species$Both$LatQIAMed)

compare_models(LM_IUCN_Species$Cross$LatQIMed,
               LM_IUCN_Species$Cross$LatQAMed,
               "IUCN data, interaction vs asymmetric, random both species")

compare_models(LM_GBIF_Species$Cross$LatQIMed,
               LM_GBIF_Species$Cross$LatQAMed,
               "GBIF data, interaction vs asymmetric, random both species")

#### Distance
# ______________________________________________________________________________
# Distance from niche centre
# yes, but it should be modelled as a quadratic

compare_models(LM_IUCN_Species$Cross$LatQIMedDist,
               LM_IUCN_Species$Cross$LatQIMed,
               "IUCN data, NicheDist added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQIMedDist,
               LM_GBIF_Species$Cross$LatQIMed,
               "GBIF data, NicheDist added to quadratic MedianProp and Latitude, random both species")


#### Quadratic Distance
# ______________________________________________________________________________
# Distance from niche centre as quadratic
# retain

compare_models(LM_IUCN_Species$Cross$LatQIMedQDist,
               LM_IUCN_Species$Cross$LatQIMed,
               "IUCN data, quadratic NicheDist added to quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQIMedQDist,
               LM_GBIF_Species$Cross$LatQIMed,
               "GBIF data, quadratic NicheDist added to quadratic MedianProp and Latitude, random both species")

#### Latitude * (Quadratic MedianProp + Quadratic Distance)
# ______________________________________________________________________________
# Distance from niche centre interacting with latitude 

compare_models(LM_IUCN_Species$Cross$LatQIMedQIDist,
               LM_IUCN_Species$Cross$LatQIMedQDist,
               "IUCN data, interaction added to quadratic NicheDist, quadratic MedianProp and Latitude, random both species")

compare_models(LM_GBIF_Species$Cross$LatQIMedQIDist,
               LM_GBIF_Species$Cross$LatQIMedQDist,
               "GBIF data, interaction added to quadratic NicheDist, quadratic MedianProp and Latitude, random both species")


#### Fixed effects Conclusions ####
# Final fixed effects model
summary(LM_IUCN_Species$Cross$LatQIMedQIDist)

fixed_effects <- data.frame(model1      = character(),
                            model2      = character(),
                            formula_m2  = character(),
                            AIC_m2      = numeric(),
                            delta_AIC   = numeric(),
                            rellik      = character(),
                            random      = character(),
                            method      = character(),
                            pseudor2m   = numeric(),
                            pseudor2c   = numeric(),
                            included    = logical())
raster_fixed_effects <- fixed_effects
asym_fixed_effects <- fixed_effects

fixed_effects <- fixed_effects %>% 
  record_fixef(m1 = "Null",          m2 = "Lat") %>% 
  record_fixef(m1 = "Lat",           m2 = "LatMed") %>% 
  record_fixef(m1 = "LatMed",        m2 = "LatQMed") %>% 
  record_fixef(m1 = "LatQMed",       m2 = "LatQIMed") %>% 
  record_fixef(m1 = "LatQIMed",      m2 = "LatQIMedDist") %>% 
  record_fixef(m1 = "LatQIMedDist",  m2 = "LatQIMedQDist") %>% 
  record_fixef(m1 = "LatQIMedQDist", m2 = "LatQIMedQIDist") 

# full gbif set, including base for easy comparison
gbif_fixed_effects <- fixed_effects %>% 
  record_fixef(m1 = "Null",          m2 = "Lat",            model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "Lat",           m2 = "LatMed",         model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "LatMed",        m2 = "LatQMed",        model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "LatQMed",       m2 = "LatQIMed",       model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "LatQIMed",      m2 = "LatQIMedDist",   model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "LatQIMedDist",  m2 = "LatQIMedQDist",  model_list = LM_GBIF_Species, method = "gbif") %>% 
  record_fixef(m1 = "LatQIMedQDist", m2 = "LatQIMedQIDist", model_list = LM_GBIF_Species, method = "gbif") 

gbif_fixed_effects <- gbif_fixed_effects %>% 
  record_fixef(m1 = "Null",          m2 = "Lat",            model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "Lat",           m2 = "LatMed",         model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "LatMed",        m2 = "LatQMed",        model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "LatQMed",       m2 = "LatQIMed",       model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "LatQIMed",      m2 = "LatQIMedDist",   model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "LatQIMedDist",  m2 = "LatQIMedQDist",  model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  record_fixef(m1 = "LatQIMedQDist", m2 = "LatQIMedQIDist", model_list = LM_IUCN_Species_B, method = "iucn_b") 

gbif_fixed_effects <- gbif_fixed_effects %>% 
  record_fixef(m1 = "Null",          m2 = "Lat",            model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "Lat",           m2 = "LatMed",         model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "LatMed",        m2 = "LatQMed",        model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "LatQMed",       m2 = "LatQIMed",       model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "LatQIMed",      m2 = "LatQIMedDist",   model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "LatQIMedDist",  m2 = "LatQIMedQDist",  model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  record_fixef(m1 = "LatQIMedQDist", m2 = "LatQIMedQIDist", model_list = LM_GBIF_Species_B, method = "gbif_b") 

# raster set
raster_fixed_effects <- raster_fixed_effects %>% 
  record_fixef(m1 = "LatQIMed",      m2 = "LatQIMedDist") %>% 
  record_fixef(m1 = "LatQIMedDist",  m2 = "LatQIMedQDist") %>% 
  record_fixef(m1 = "LatQIMedQDist", m2 = "LatQIMedQIDist") %>% 
  record_fixef(m1 = "LatQIMed",        m2 = "LatQIMedDist_5") %>% 
  record_fixef(m1 = "LatQIMedDist_5",  m2 = "LatQIMedQDist_5") %>% 
  record_fixef(m1 = "LatQIMedQDist_5", m2 = "LatQIMedQIDist_5") %>% 
  record_fixef(m1 = "LatQIMed",         m2 = "LatQIMedDist_20") %>% 
  record_fixef(m1 = "LatQIMedDist_20",  m2 = "LatQIMedQDist_20") %>% 
  record_fixef(m1 = "LatQIMedQDist_20", m2 = "LatQIMedQIDist_20") 


# asym set
asym_fixed_effects <- asym_fixed_effects %>% 
  record_fixef(m1 = "LatQMed",  m2 = "LatQAMed") %>% 
  record_fixef(m1 = "LatQMed",  m2 = "LatQIMed") %>% 
  record_fixef(m1 = "LatQAMed", m2 = "LatQIMed") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQAMed", model_list = LM_GBIF_Species, method = "gbif") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQIMed", model_list = LM_GBIF_Species, method = "gbif") %>% 
  # record_fixef(m1 = "LatQAMed", m2 = "LatQIMed", model_list = LM_GBIF_Species, method = "gbif") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQAMed", model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQIMed", model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  # record_fixef(m1 = "LatQAMed", m2 = "LatQIMed", model_list = LM_IUCN_Species_B, method = "iucn_b") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQAMed", model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  # record_fixef(m1 = "LatQMed",  m2 = "LatQIMed", model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  # record_fixef(m1 = "LatQAMed", m2 = "LatQIMed", model_list = LM_GBIF_Species_B, method = "gbif_b") %>% 
  mutate(m1m2 = paste(model1, model2, sep = " ")) 

# fixed_effects <- fixed_effects %>% 
#   mutate(modelname = case_when(model2 == "Lat"            ~ "Latitude",
#                                model2 == "LatMed"         ~ "Latitude $+$ Range position",
#                                model2 == "LatQMed"        ~ "Latitude $+$ Range position$^2$",
#                                model2 == "LatQIMed"       ~ "Latitude \\texttimes\\space Range position$^2$",
#                                model2 == "LatQIMedDist"   ~ "Latitude \\texttimes\\space Range position$^2$ $+$ Niche position",
#                                model2 == "LatQIMedQDist"  ~ "Latitude \\texttimes\\space Range position$^2$ $+$ Niche position$^2$",
#                                model2 == "LatQIMedQIDist" ~ "Latitude \\texttimes\\space (Range position$^2$ $+$ Niche position$^2$)")) %>% 
#   select(-formula_m2)
fixed_effects <- fixed_effects %>% 
  mutate(modelname = case_when(model2 == "Lat"            ~ "L",
                               model2 == "LatMed"         ~ "L $+$ R",
                               model2 == "LatQMed"        ~ "L $+$ R$^2$",
                               model2 == "LatQIMed"       ~ "L \\texttimes\\space R$^2$",
                               model2 == "LatQIMedDist"   ~ "L \\texttimes\\space R$^2$ $+$ N",
                               model2 == "LatQIMedQDist"  ~ "L \\texttimes\\space R$^2$ $+$ N$^2$",
                               model2 == "LatQIMedQIDist" ~ "L \\texttimes\\space (R$^2$ $+$ N$^2$)")) %>% 
  select(-formula_m2)

write.csv(fixed_effects,        here::here("Tables", "fixed_effects.csv"),        row.names = FALSE, na = "")
write.csv(gbif_fixed_effects,   here::here("Tables", "gbif_fixed_effects.csv"),   row.names = FALSE, na = "")
write.csv(raster_fixed_effects, here::here("Tables", "raster_fixed_effects.csv"), row.names = FALSE, na = "")
write.csv(asym_fixed_effects,   here::here("Tables", "asym_fixed_effects.csv"),   row.names = FALSE, na = "")

fixef_coef <- c("Null", "Lat", "LatMed", "LatQMed", "LatQIMed", 
                "LatQIMedDist", "LatQIMedQDist", "LatQIMedQIDist") %>% 
  sapply(function(name){
    fixef(LM_IUCN_Species$Cross[[name]])
  }, USE.NAMES=TRUE) %>% 
  bind_rows(.id = "model") %>% 
  rename(Intercept       = `(Intercept)`, Latitude = LatitudeScaled,
         MedianProp2     = `poly(MedianPropSquared, 2)2`,
         Lat.MedianProp1 = `LatitudeScaled:poly(MedianPropSquared, 2)1`,
         Lat.MedianProp2 = `LatitudeScaled:poly(MedianPropSquared, 2)2`,
         DistProp2       = `poly(CorrectedDistPropSquared, 2)2`,
         Lat.DistProp1   = `LatitudeScaled:poly(CorrectedDistPropSquared, 2)1`,
         Lat.DistProp2   = `LatitudeScaled:poly(CorrectedDistPropSquared, 2)2`) %>% 
  mutate(MedianProp1     = coalesce(MedianPropSquScaled, `poly(MedianPropSquared, 2)1`),
         DistProp1       = coalesce(CorrectedDistPropSquScaled, `poly(CorrectedDistPropSquared, 2)1`)) %>% 
  select(model, Intercept, Latitude, MedianProp1, MedianProp2, Lat.MedianProp1, Lat.MedianProp2,
         DistProp1, DistProp2, Lat.DistProp1, Lat.DistProp2) %>% 
  mutate(modelname = case_when(model == "Null"           ~ "Null",
                               model == "Lat"            ~ "L",
                               model == "LatMed"         ~ "L $+$ R",
                               model == "LatQMed"        ~ "L $+$ R$^2$",
                               model == "LatQIMed"       ~ "L \\texttimes\\space R$^2$",
                               model == "LatQIMedDist"   ~ "L \\texttimes\\space R$^2$ $+$ N",
                               model == "LatQIMedQDist"  ~ "L \\texttimes\\space R$^2$ $+$ N$^2$",
                               model == "LatQIMedQIDist" ~ "L \\texttimes\\space (R$^2$ $+$ N$^2$)"))

write.csv(fixef_coef, here::here("Tables", "fixef_coef1.csv"), row.names = FALSE, na = "")

fixef_coef <- fixef_coef %>% 
  select(-modelname) %>% 
  pivot_longer(!model, names_to = "col1", values_to = "col2") %>% 
  pivot_wider(names_from = "model", values_from = "col2") %>% 
  mutate(term = case_when(col1 == "MedianProp1"     ~ "R",
                          col1 == "MedianProp2"     ~ "R$^2$",
                          col1 == "Lat.MedianProp1" ~ "L:R",
                          col1 == "Lat.MedianProp2" ~ "L:R$^2$",
                          col1 == "DistProp1"       ~ "N",
                          col1 == "DistProp2"       ~ "N$^2$",
                          col1 == "Lat.DistProp1"   ~ "L:N",
                          col1 == "Lat.DistProp2"   ~ "L:N$^2$",
                          col1 == "Intercept"       ~ "I",
                          col1 == "Latitude"        ~ "L",
                          TRUE                      ~ col1))

write.csv(fixef_coef, here::here("Tables", "fixef_coef2.csv"), row.names = FALSE, na = "")



# in comparison to gbif there are some different effect size estimates...
# Nothing too drastic but they just look a little off
summary(LM_GBIF_Species$Cross$LatQIMedQIDist)
plot_latitude(model_list = LM_IUCN_Species, model_data = GMPD_IUCN_Species)
plot_latitude(model_list = LM_GBIF_Species, model_data = GMPD_GBIF_Species)

plot_medianprop(model_list = LM_IUCN_Species, model_data = GMPD_IUCN_Species, fixed_lat = 25)
plot_medianprop(model_list = LM_GBIF_Species, model_data = GMPD_GBIF_Species, fixed_lat = 25)

plot_distprop(model_list = LM_IUCN_Species, model_data = GMPD_IUCN_Species, fixed_lat = 25)
plot_distprop(model_list = LM_GBIF_Species, model_data = GMPD_GBIF_Species, fixed_lat = 25)

# basically the same across raster resolutions 
summary(LM_IUCN_Species$Both$LatQIMedQIDist)
summary(LM_IUCN_Species$Both$LatQIMedQIDist_5)
summary(LM_IUCN_Species$Both$LatQIMedQIDist_20)


plot_latitude(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIDist_5", 
              model_data = GMPD_Both_Species,
              raster_res = "_5")
plot_latitude(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species)
plot_latitude(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIDist_20", 
              model_data = GMPD_Both_Species,
              raster_res = "_20")

plot_medianprop(model_list = LM_IUCN_Species, 
                model_name = "LatQIMedQIDist_5", 
                model_data = GMPD_Both_Species, 
                fixed_lat = 40,
                raster_res = "_5")
plot_medianprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40)
plot_medianprop(model_list = LM_IUCN_Species, 
                model_name = "LatQIMedQIDist_20", 
                model_data = GMPD_Both_Species, 
                fixed_lat = 40,
                raster_res = "_20")

plot_distprop(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIDist_5", 
              model_data = GMPD_Both_Species, 
              fixed_lat = 40,
              raster_res = "_5")
plot_distprop(model_list = LM_IUCN_Species, model_data = GMPD_Both_Species, fixed_lat = 40)
plot_distprop(model_list = LM_IUCN_Species, 
              model_name = "LatQIMedQIDist_20", 
              model_data = GMPD_Both_Species, 
              fixed_lat = 40,
              raster_res = "_20")

# diagnostics
sim_LatQIMedQIDist_Both    <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Both$LatQIMedQIDist)

plot(sim_LatQIMedQIDist_Both)



### Random Effects double check and diagnostics ####
# When looking at the random effects the maximal fixed effects model was:
summary(LM_IUCN_Species$Cross$LatQIMedQIDist)
summary(LM_GBIF_Species$Cross$LatQIMedQIDist)
# This repeatedly failed to converge when we tried to model random slopes so we stick to random intercepts

# anova comparisons
full_random_anova <- list()

full_random_anova$Host_Fixed        <- anova(LM_IUCN_Species$Host$LatQIMedQIDist,       LM_IUCN_Species$Fixed$LatQIMedQIDist)
full_random_anova$Group_Fixed       <- anova(LM_IUCN_Species$Group$LatQIMedQIDist,      LM_IUCN_Species$Fixed$LatQIMedQIDist)
full_random_anova$HostGroup_Host    <- anova(LM_IUCN_Species$HostGroup$LatQIMedQIDist,  LM_IUCN_Species$Host$LatQIMedQIDist)
full_random_anova$HostGroup_Group   <- anova(LM_IUCN_Species$HostGroup$LatQIMedQIDist,  LM_IUCN_Species$Group$LatQIMedQIDist)

full_random_anova$Parasite_Fixed    <- anova(LM_IUCN_Species$Parasite$LatQIMedQIDist,   LM_IUCN_Species$Fixed$LatQIMedQIDist)
full_random_anova$Type_Fixed        <- anova(LM_IUCN_Species$Type$LatQIMedQIDist,       LM_IUCN_Species$Fixed$LatQIMedQIDist)
full_random_anova$ParType_Parasite  <- anova(LM_IUCN_Species$ParType$LatQIMedQIDist,    LM_IUCN_Species$Parasite$LatQIMedQIDist)
full_random_anova$ParType_Type      <- anova(LM_IUCN_Species$ParType$LatQIMedQIDist,    LM_IUCN_Species$Type$LatQIMedQIDist)

full_random_anova$Both_Host         <- anova(LM_IUCN_Species$Both$LatQIMedQIDist,       LM_IUCN_Species$Host$LatQIMedQIDist)
full_random_anova$Both_Parasite     <- anova(LM_IUCN_Species$Both$LatQIMedQIDist,       LM_IUCN_Species$Parasite$LatQIMedQIDist)
full_random_anova$Cross_Both        <- anova(LM_IUCN_Species$Cross$LatQIMedQIDist,      LM_IUCN_Species$Both$LatQIMedQIDist)

# simulate residuals
# TODO: put all models, tests, diagnostics in same list
sim_LatQIMedQIDist <- list()

sim_LatQIMedQIDist$Fixed     <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Fixed$LatQIMedQIDist)
sim_LatQIMedQIDist$Host      <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Host$LatQIMedQIDist)
# sim_LatQIMedQIDist$Group     <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Group$LatQIMedQIDist)
# sim_LatQIMedQIDist$HostGroup <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$HostGroup$LatQIMedQIDist)

sim_LatQIMedQIDist$Parasite  <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Parasite$LatQIMedQIDist)
# sim_LatQIMedQIDist$Type      <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Type$LatQIMedQIDist)
# sim_LatQIMedQIDist$ParType   <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$ParType$LatQIMedQIDist)

sim_LatQIMedQIDist$Both      <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Both$LatQIMedQIDist)
sim_LatQIMedQIDist$Cross     <- DHARMa::simulateResiduals(fittedModel = LM_IUCN_Species$Cross$LatQIMedQIDist)

# diagnostic tests

full_random_diagnostics <- list()

full_random_diagnostics$Dispersion$Fixed      <- capture.output(testDispersion(sim_LatQIMedQIDist$Fixed, alternative = "greater"))[[5]]
full_random_diagnostics$Dispersion$Host       <- capture.output(testDispersion(sim_LatQIMedQIDist$Host, alternative = "greater"))[[5]]
# full_random_diagnostics$Dispersion$Group      <- capture.output(testDispersion(sim_LatQIMedQIDist$Group))[[5]]
# full_random_diagnostics$Dispersion$HostGroup  <- capture.output(testDispersion(sim_LatQIMedQIDist$HostGroup))[[5]]

full_random_diagnostics$Dispersion$Parasite   <- capture.output(testDispersion(sim_LatQIMedQIDist$Parasite, alternative = "greater"))[[5]]
# full_random_diagnostics$Dispersion$Type       <- capture.output(testDispersion(sim_LatQIMedQIDist$Type))[[5]]
# full_random_diagnostics$Dispersion$ParType    <- capture.output(testDispersion(sim_LatQIMedQIDist$ParType))[[5]]

full_random_diagnostics$Dispersion$Both       <- capture.output(testDispersion(sim_LatQIMedQIDist$Both, alternative = "greater"))[[5]]
full_random_diagnostics$Dispersion$Cross      <- capture.output(testDispersion(sim_LatQIMedQIDist$Cross, alternative = "greater"))[[5]]

#-
full_random_diagnostics$ZeroInflation$Fixed      <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Fixed))[[5]]
full_random_diagnostics$ZeroInflation$Host       <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Host))[[5]]
# full_random_diagnostics$ZeroInflation$Group      <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Group))[[5]]
# full_random_diagnostics$ZeroInflation$HostGroup  <- capture.output(testZeroInflation(sim_LatQIMedQIDist$HostGroup))[[5]]

full_random_diagnostics$ZeroInflation$Parasite   <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Parasite))[[5]]
# full_random_diagnostics$ZeroInflation$Type       <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Type))[[5]]
# full_random_diagnostics$ZeroInflation$ParType    <- capture.output(testZeroInflation(sim_LatQIMedQIDist$ParType))[[5]]

full_random_diagnostics$ZeroInflation$Both       <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Both))[[5]]
full_random_diagnostics$ZeroInflation$Cross      <- capture.output(testZeroInflation(sim_LatQIMedQIDist$Cross))[[5]]

# general notes for planned reshuffle:
# overdispersion is bad at the start, then evens out and is okay once random effects are added
# zero inflation is always bad because we have basically no zeros
# I did not test for outliers because they would probably be a result of either dispersion or zero inflation

##### Fixed
summary(LM_IUCN_Species$Fixed$LatQIMedQIDist)

plotResiduals(sim_LatQIMedQIDist$Fixed)
testDispersion(sim_LatQIMedQIDist$Fixed, alternative = "greater")
testZeroInflation(sim_LatQIMedQIDist$Fixed)
# badly overdispersed
# zero deflated

##### Host
summary(LM_IUCN_Species$Host$LatQIMedQIDist)
full_random_anova$Host_Fixed
# Host is definitely better than nothing 

plotResiduals(sim_LatQIMedQIDist$Host)
testDispersion(sim_LatQIMedQIDist$Host, alternative = "greater")
testZeroInflation(sim_LatQIMedQIDist$Host)
# Pretty strong residual trend introduced by host, but dispersion is way better

##### Parasite species
summary(LM_IUCN_Species$Parasite$LatQIMedQIDist)
full_random_anova$Parasite_Fixed
# defo better than nothing

plotResiduals(sim_LatQIMedQIDist$Parasite)
testDispersion(sim_LatQIMedQIDist_Parasite, alternative = "greater")
testZeroInflation(sim_LatQIMedQIDist_Parasite)
# Much better with parasite
# The slope doesn't show up
# Dispersion is still greatly improved

##### Host + Parasite
summary(LM_IUCN_Species$Both$LatQIMedQIDist)
full_random_anova$Both_Host
full_random_anova$Both_Parasite

# Better than host, better than parasite

plotResiduals(sim_LatQIMedQIDist$Both)
testDispersion(sim_LatQIMedQIDist_Both, alternative = "greater")
testZeroInflation(sim_LatQIMedQIDist_Both)
# The slope in residuals returns, albeit less severely
# Dispersion is still fine though

##### Host:Parasite
summary(LM_IUCN_Species$Cross$LatQIMedQIDist)
full_random_anova$Cross_Both

# Better than host, better than parasite species, better than both

plotResiduals(sim_LatQIMedQIDist$Cross)
testDispersion(sim_LatQIMedQIDist_Cross, alternative = "greater")
testZeroInflation(sim_LatQIMedQIDist_Cross)

# Looks okay

#### Random effect conclusions ####
if(plot_outputs){
  pdf(here::here("Figures/random effect diagnostics.pdf"), width = 10, height = 7)
  plotResiduals(sim_LatQIMedQIDist$Fixed)
  mtext("No random terms", side=3, line = 0.5, adj = 0.5)
  plotResiduals(sim_LatQIMedQIDist$Host)
  mtext("Host random intercept", side=3, line = 0.5, adj = 0.5)
  plotResiduals(sim_LatQIMedQIDist$Parasite)
  mtext("Parasite random intercept", side=3, line = 0.5, adj = 0.5)
  plotResiduals(sim_LatQIMedQIDist$Both)
  mtext("Both host and parasite random intercept", side=3, line = 0.5, adj = 0.5)
  plotResiduals(sim_LatQIMedQIDist$Cross)
  mtext("Host, parasite, and interaction between host and parasite species random intercept", side=3, line = 0.5, adj = 0.5)
  dev.off()
}

# full_random_effects <- c("Host",      "Fixed",
#                          "Parasite",  "Fixed",
#                          "Both",      "Parasite",
#                          "Both",      "Host",
#                          "Cross",     "Both") %>% 
#   matrix(ncol = 2, byrow = TRUE) %>% 
#   as.data.frame() %>% 
#   rename("model1" = "V1", "model2" = "V2") %>% 
#   mutate(paste_names = paste(model1, model2, sep = "_")) %>% 
#   group_by(paste_names) %>% 
#   mutate(deltaDF = full_random_anova[[paste_names]]$Df[[2]],
#          m1 = str_split_i(rownames(full_random_anova[[paste_names]])[[1]], "\\$", 2),
#          m2 = str_split_i(rownames(full_random_anova[[paste_names]])[[2]], "\\$", 2),
#          
#          npar1 = case_when(model1 == m1 ~ full_random_anova[[paste_names]]$npar[[1]],
#                            model1 == m2 ~ full_random_anova[[paste_names]]$npar[[2]]),
#          npar2 = case_when(model2 == m2 ~ full_random_anova[[paste_names]]$npar[[2]],
#                            model2 == m1 ~ full_random_anova[[paste_names]]$npar[[1]]),
#          
#          loglike1 = case_when(model1 == m1 ~ full_random_anova[[paste_names]]$logLik[[1]],
#                               model1 == m2 ~ full_random_anova[[paste_names]]$logLik[[2]]),
#          loglike2 = case_when(model2 == m2 ~ full_random_anova[[paste_names]]$logLik[[2]],
#                               model2 == m1 ~ full_random_anova[[paste_names]]$logLik[[1]]),
#          
#          Chisq = full_random_anova[[paste_names]]$Chisq[[2]],
#          P = full_random_anova[[paste_names]]$`Pr(>Chisq)`[[2]]) %>% 
#   ungroup() %>% 
#   select(-paste_names, -m1, -m2)

full_random_effects <- null_random_effects %>% 
  select(model1, model2) %>% 
  mutate(paste_names = paste(model1, model2, sep = "_")) %>% 
  group_by(paste_names) %>% 
  mutate(deltaDF = full_random_anova[[paste_names]]$Df[[2]],
         m1 = str_split_i(rownames(full_random_anova[[paste_names]])[[1]], "\\$", 2),
         m2 = str_split_i(rownames(full_random_anova[[paste_names]])[[2]], "\\$", 2),
         
         npar1 = case_when(model1 == m1 ~ full_random_anova[[paste_names]]$npar[[1]],
                           model1 == m2 ~ full_random_anova[[paste_names]]$npar[[2]]),
         npar2 = case_when(model2 == m2 ~ full_random_anova[[paste_names]]$npar[[2]],
                           model2 == m1 ~ full_random_anova[[paste_names]]$npar[[1]]),
         
         loglike1 = case_when(model1 == m1 ~ full_random_anova[[paste_names]]$logLik[[1]],
                              model1 == m2 ~ full_random_anova[[paste_names]]$logLik[[2]]),
         loglike2 = case_when(model2 == m2 ~ full_random_anova[[paste_names]]$logLik[[2]],
                              model2 == m1 ~ full_random_anova[[paste_names]]$logLik[[1]]),
         
         Chisq = full_random_anova[[paste_names]]$Chisq[[2]],
         P = null_random_anova[[paste_names]]$`Pr(>Chisq)`[[2]],
         P = case_when(P < 0.001 ~ "< 0.001",
                       TRUE ~ as.character(round(P, 2)))) %>% 
  ungroup() %>% 
  select(-paste_names, -m1, -m2)

all_random_effects <- full_join(null_random_effects, full_random_effects,
                                by = join_by(model1, model2),
                                suffix = c("_Null", "_Full"))

all_random_effects <- all_random_effects %>% 
  mutate(name1 = case_when(model1 == "Parasite"  ~ "Ps",
                           model1 == "Host"      ~ "Hs",
                           model1 == "Group"     ~ "Hg",
                           model1 == "HostGroup" ~ "Hg/Hs",
                           model1 == "Type"      ~ "Pt",
                           model1 == "ParType"   ~ "Pt/Ps",
                           model1 == "Both"      ~ "Hs $+$ Ps",
                           model1 == "Cross"     ~ "Hs \\texttimes\\space Ps"),
         name2 = case_when(model2 == "Fixed"     ~ "None",
                           model2 == "Parasite"  ~ "Ps",
                           model2 == "Host"      ~ "Hs",
                           model2 == "Group"     ~ "Hg",
                           model2 == "HostGroup" ~ "Hg/Hs",
                           model2 == "Type"      ~ "Pt",
                           model2 == "ParType"   ~ "Pt/Ps",
                           model2 == "Both"      ~ "Hs $+$ Ps"))

write.csv(full_random_effects, here::here("Tables", "full_random_effects.csv"), row.names = FALSE, na = "")
write.csv(all_random_effects,  here::here("Tables", "all_random_effects.csv"),  row.names = FALSE, na = "")

full_random_stdev <- null_random_stdev %>% 
  select(model1) %>% 
  distinct() %>% 
  group_by(model1) %>% 
  summarise(term = as.data.frame(VarCorr(LM_IUCN_Species[[as.character(model1)]]$LatQIMedQIDist))$grp,
            sdcor_full = as.data.frame(VarCorr(LM_IUCN_Species[[as.character(model1)]]$LatQIMedQIDist))$sdcor) 

all_random_stdev <- full_join(null_random_stdev, full_random_stdev, by = c("model1", "term"))
all_random_stdev[match(unique(all_random_stdev$model1), all_random_stdev$model1), 
                 "keep_name"] <- TRUE
all_random_stdev <- all_random_stdev %>% 
  mutate(name1 = case_when(model1 == "Parasite"  ~ "Parasite species",
                           model1 == "Host"      ~ "Host species",
                           model1 == "Group"     ~ "Host group",
                           model1 == "HostGroup" ~ "Host group/species",
                           model1 == "Type"      ~ "Parasite type",
                           model1 == "ParType"   ~ "Parasite type/species",
                           model1 == "Both"      ~ "Host $+$ parasite species",
                           model1 == "Cross"     ~ "Host \\texttimes\\space parasite species"),
         name1 = case_when(keep_name ~ name1,
                            TRUE      ~ ""),
         term_name = str_replace_all(term, "CorrectedName", ""),
         term_name = str_replace_all(term_name, "ParType", "Type")) %>% 
  select(-keep_name)

write.csv(all_random_stdev,  here::here("Tables", "all_random_stdev.csv"),  row.names = FALSE, na = "")


full_diagnostics <- data.frame(modelRanef = c("Fixed",
                                              "Host",
                                              "Parasite",
                                              "Both",
                                              "Cross")) %>% 
  group_by(modelRanef) %>% 
  mutate(dispersion  = str_extract(full_random_diagnostics[["Dispersion"]][[modelRanef]], 
                                   pattern = "(?<== )[[:graph:]]*(?=,)"),
         dispersion_P = str_extract(full_random_diagnostics[["Dispersion"]][[modelRanef]], 
                                   pattern = "(?<= )[[:graph:]]*$"),
         
         zeros  = str_extract(full_random_diagnostics[["ZeroInflation"]][[modelRanef]], 
                              pattern = "(?<== )[[:graph:]]*(?=,)"),
         zeros_P = str_extract(full_random_diagnostics[["ZeroInflation"]][[modelRanef]], 
                              pattern = "(?<= )[[:graph:]]*$")) %>% 
  ungroup()

### Host Phylogeny ####
library(U.PhyloMaker)

megatree <- read.tree('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_megatree.tre')
sp.list  <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_sample_species_list.csv', sep=",")
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/mammal_20221117/main/mammal_genus_list.csv', sep=",")

phylo_model_data <- GMPD_IUCN_Species %>% 
  mutate(species = case_when(HostCorrectedName == "Cervus canadensis"      ~ "Cervus elaphus",
                             HostCorrectedName == "Lycalopex gymnocercus"  ~ "Pseudalopex gymnocercus",
                             TRUE ~ HostCorrectedName)) 

sp.list <- phylo_model_data %>% 
  mutate(genus = str_split_i(species, " ", 1)) %>% 
  select(species, genus) %>% 
  distinct()


result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)

# Create a distance matrix for phylo in residuals
phyloMat <- cophenetic.phylo(result) %>% as.data.frame.array() 

phylo_models <- list()
phylo_models$LatQIMedQIDist$Fixed       <- try(        glm(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistPropSquared, 2)),                                                                 data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$LatQIMedQIDist$Host        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistPropSquared, 2)) + (1|species),                                                                 data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$LatQIMedQIDist$Parasite    <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistPropSquared, 2)) + (1|ParasiteCorrectedName),                                                   data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$LatQIMedQIDist$Both        <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistPropSquared, 2)) + (1|species) + (1|ParasiteCorrectedName),                                     data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)
phylo_models$LatQIMedQIDist$Cross       <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * (poly(MedianPropSquared_iucn, 2) + poly(CorrectedDistPropSquared, 2)) + (1|species) + (1|ParasiteCorrectedName) + (1|species:ParasiteCorrectedName), data = phylo_model_data, family = binomial, weights = SampleSize), silent = TRUE)

phylo_models$sim$Fixed    <- DHARMa::simulateResiduals(fittedModel = phylo_models$LatQIMedQIDist$Fixed)
phylo_models$sim$Host     <- DHARMa::simulateResiduals(fittedModel = phylo_models$LatQIMedQIDist$Host)
phylo_models$sim$Parasite <- DHARMa::simulateResiduals(fittedModel = phylo_models$LatQIMedQIDist$Parasite)
phylo_models$sim$Both     <- DHARMa::simulateResiduals(fittedModel = phylo_models$LatQIMedQIDist$Both)
phylo_models$sim$Cross    <- DHARMa::simulateResiduals(fittedModel = phylo_models$LatQIMedQIDist$Cross)

phylo_models$hostrecalc$Fixed     <- DHARMa::recalculateResiduals(phylo_models$sim$Fixed,     group = phylo_model_data$species)
phylo_models$hostrecalc$Host      <- DHARMa::recalculateResiduals(phylo_models$sim$Host,      group = phylo_model_data$species)
phylo_models$hostrecalc$Parasite  <- DHARMa::recalculateResiduals(phylo_models$sim$Parasite,  group = phylo_model_data$species)
phylo_models$hostrecalc$Both      <- DHARMa::recalculateResiduals(phylo_models$sim$Both,      group = phylo_model_data$species)
phylo_models$hostrecalc$Cross     <- DHARMa::recalculateResiduals(phylo_models$sim$Cross,     group = phylo_model_data$species)

phylo_models$SpatialAuto$Fixed    <- capture.output(testSpatialAutocorrelation(simulationOutput = phylo_models$hostrecalc$Fixed,     distMat = phyloMat))[[5]]
phylo_models$SpatialAuto$Host     <- capture.output(testSpatialAutocorrelation(simulationOutput = phylo_models$hostrecalc$Host,      distMat = phyloMat))[[5]]
phylo_models$SpatialAuto$Parasite <- capture.output(testSpatialAutocorrelation(simulationOutput = phylo_models$hostrecalc$Parasite,  distMat = phyloMat))[[5]]
phylo_models$SpatialAuto$Both     <- capture.output(testSpatialAutocorrelation(simulationOutput = phylo_models$hostrecalc$Both,      distMat = phyloMat))[[5]]
phylo_models$SpatialAuto$Cross    <- capture.output(testSpatialAutocorrelation(simulationOutput = phylo_models$hostrecalc$Cross,     distMat = phyloMat))[[5]]

if(plot_outputs){
  pdf(here::here("Figures/phylogeny.pdf"), width = 10, height = 7)
  plotResiduals(phylo_models$hostrecalc$Fixed)
  mtext("No random intercepts",       side=3, line = 0.5, adj = 0.5)
  plotResiduals(phylo_models$hostrecalc$Host)
  mtext("Host random intercepts",     side=3, line = 0.5, adj = 0.5)
  plotResiduals(phylo_models$hostrecalc$Parasite)
  mtext("Parasite random intercepts", side=3, line = 0.5, adj = 0.5)
  plotResiduals(phylo_models$hostrecalc$Both)
  mtext("Both random intercepts",     side=3, line = 0.5, adj = 0.5)
  plotResiduals(phylo_models$hostrecalc$Cross)
  mtext("Cross random intercepts",    side=3, line = 0.5, adj = 0.5)
  dev.off()
}

##### Fixed
# No phylogenetic signal

##### Host
# definitely no signal here
# but there is heteroskedasticity

##### Parasite
# although there's no phylogenetic signal, parasites seem to be introducing something into the residuals vs predicted
# might be because different hosts have different sets of parasites
# either way, there is no reason to include phylogeny and the diagnostics without phylogeny look fine too

##### Both
# no strong phylogenetic signal
# although there appears to be a slope in the plot

##### Cross
# no strong phylogenetic signal
# although there appears to be a slope in the plot

full_diagnostics <- data.frame(model = c("Fixed", "Host", "Parasite", "Both", "Cross")) %>% 
  group_by(model) %>% 
  mutate(observed = str_extract(phylo_models$SpatialAuto[[model]],
                                pattern = "(?<=observed = )[[:graph:]]*(?=,)"),
         expected = str_extract(phylo_models$SpatialAuto[[model]],
                                pattern = "(?<=expected = )[[:graph:]]*(?=,)"),
         SD       = str_extract(phylo_models$SpatialAuto[[model]],
                                pattern = "(?<=sd = )[[:graph:]]*(?=,)"),
         phylo_P   = str_extract(phylo_models$SpatialAuto[[model]],
                                pattern = "(?<=value = )[[:graph:]]*$")) %>% 
  ungroup() %>% 
  full_join(full_diagnostics, by = join_by(model == modelRanef)) %>% 
  mutate(itemname  = case_when(model == "Fixed"    ~ "No",
                               model == "Host"     ~ "Host species",
                               model == "Parasite" ~ "Parasite species",
                               model == "Both"     ~ "Host and parasite species",
                               model == "Cross"    ~ "Host, parasite, and interaction between host and parasite species"),
         tablename = case_when(model == "Fixed"    ~ "None",
                               model == "Host"     ~ "Host species",
                               model == "Parasite" ~ "Parasite species",
                               model == "Both"     ~ "Host $+$ parasite species",
                               model == "Cross"    ~ "Host \\texttimes\\space parasite species")) %>% 
  mutate(phylo_P = case_when(as.numeric(phylo_P) < 0.01 ~ "< 0.01",
                             TRUE ~ as.character(round(as.numeric(phylo_P), 2))),
         dispersion_P = case_when(as.numeric(dispersion_P) < 0.01 ~ "< 0.01",
                             TRUE ~ as.character(round(as.numeric(dispersion_P), 2))),
         zeros_P = case_when(as.numeric(zeros_P) < 0.01 ~ "< 0.01",
                             TRUE ~ as.character(round(as.numeric(zeros_P), 2))),
         dispersion = as.numeric(dispersion))

write.table(full_diagnostics, here::here("Tables", "full_diagnostics.csv"), 
            row.names = FALSE, sep = ";", na = "")

# Plots ########################################################################

### Fig 1, map ####
# map

if(plot_outputs){
  IUCN_Native_Data <- readRDS(here::here("Data/Data back ups/IUCN_Native_Data_01"))
  Projection_String <- sf::st_crs(IUCN_Native_Data)
  Mollweide_String  <- sf::st_crs("+proj=moll")
  
  # Map
  raster.count <- raster::raster(resolution = 2/6)
  sf::sf_use_s2(FALSE)
  tmap::tmap_mode("plot")
  data("World", package = "tmap")
  bg_colour <- fasterize::fasterize(World,
                                    raster = raster.count, 
                                    background = 0) %>% 
    terra::values()
  World <- World %>%
    filter(continent != "Antarctica") %>%
    sf::st_transform(crs = Mollweide_String)
  
  GMPD_Map_Data <- sf::st_as_sf(GMPD_Climate_Data,
                                coords = c("Longitude", "Latitude"),
                                crs = Projection_String) %>% 
    sf::st_transform(crs = Mollweide_String)
  
  # populate base raster
  fasterized <- IUCN_Native_Data %>%
    mutate(Group = case_when(order_ == "CARNIVORA" ~ "Carnivores",
                             TRUE ~ "Ungulates")) %>%
    fasterize::fasterize(raster = raster.count,
                         fun = "count",
                         by = "Group",
                         background = 0)
  spatfaster <- as(fasterized, "SpatRaster")
  spatfaster[["total"]]    <- spatfaster[["Carnivores"]] + spatfaster[["Ungulates"]]
  spatfaster[["ung_prop"]] <- spatfaster[["Ungulates"]]/spatfaster[["total"]]
  spatfaster[["ung_prop"]][is.na(spatfaster[["ung_prop"]])] <- 0
  spatfaster[["opacity"]]  <- spatfaster[["total"]]/terra::global(spatfaster[["total"]], fun = "max", na.rm = TRUE)[[1]]

  # make my colours  
  My_colours <- colorRamp(c("#1B998B", "#BEEF9E"))
  col_tab <- terra::values(spatfaster[["ung_prop"]]) %>% 
    data.frame(ung_prop = ., 
               My_colours(.), 
               opacity = terra::values(spatfaster[["opacity"]]),
               bg = bg_colour) %>% 
    mutate(hex_col = rgb(X1, X2, X3, opacity*255, maxColorValue = 255),
           hex_col2 = rgb(X1, X2, X3, maxColorValue = 255)) %>% 
    mutate(bg2 = case_when(bg == 0 ~ "#FFFFFF",
                           bg == 1 ~ "#CCCCCC"),
           hex_final = mix_colours(hex_col2, bg2, 0.2+(opacity*(1-0.2)) )) 
  
  # make raster to plot
  plotting_raster <- as(raster.count, "SpatRaster")
  terra::values(plotting_raster) <- col_tab %>% 
    mutate(hex_final = case_when(opacity == 0 ~ NA_character_,
                                 TRUE         ~ hex_final),
           hex_final = as.factor(hex_final)) %>% 
    pull(hex_final) 
  
  # plot it
  tmap_options(max.categories = 250)
  
  cairo_pdf(here::here("Figures/map_moll2.pdf"), width = 10, height = 5)
  tm_shape(World) + tm_fill() +
    tm_layout(bg.color = "white", fontfamily = "Outfit",
              legend.position = c("left", "center"),
              frame = FALSE) +
    tm_shape(plotting_raster) + 
    tm_raster(palette = unique(sort(col_tab$hex_final)),
              legend.show = FALSE) +
    tm_shape(GMPD_Map_Data) + tm_dots(col = "SampleSize", 
                                      palette = c("#EDAE49",
                                                  "#D1495B",
                                                  "#222E50"
                                      ),
                                      title = "Sample size", 
                                      style = "log10") +
    tm_add_legend(type = "fill", 
                  col = rgb(My_colours(c(0, 0.5, 1)), maxColorValue = 255), 
                  title = "Taxonomic group",
                  labels = c("Ungulates", "Mixed", "Carnivores"),
                  border.alpha = 0) +
    tm_add_legend(type = "fill",
                  col = mix_colours(rgb(My_colours(0.5), maxColorValue = 255),
                                    "#FFFFFF", amount = c(1, 0.6, 0.22)),
                  title = "Number of\nhost ranges",
                  labels = c("27", "13", "1"),
                  border.alpha = 0)
  dev.off()
}


## Fig 2, lat ####
# main plot
if(plot_outputs){

  Latitude_main <- plot_latitude(model_list = LM_IUCN_Species, 
                                 model_data = GMPD_IUCN_Species) +
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
                                    model_data = GMPD_IUCN_Species, 
                                    fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "A, latitude = 25",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_40 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "C, latitude = 40",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_55 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 55) +
    # theme(axis.text.x = element_blank()) +
    labs(x = "Range position", y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E, latitude = 55",
             family = "Outfit", size = 5, hjust = 0)
  
  
  DistProp_25  <- plot_distprop(model_list = LM_IUCN_Species, 
                                model_data = GMPD_IUCN_Species, 
                                fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_40 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_data = GMPD_IUCN_Species, 
                               fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_55 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_data = GMPD_IUCN_Species, 
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

## Fig 4, range position ####
if(plot_outputs){
  MedianProp_25  <- plot_medianprop(model_list = LM_IUCN_Species, 
                                    model_data = GMPD_IUCN_Species, 
                                    fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "A, latitude = 25",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_40 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "C, latitude = 40",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_55 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 55) +
    # theme(axis.text.x = element_blank()) +
    labs(x = "Range position, full model", y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E, latitude = 55",
             family = "Outfit", size = 5, hjust = 0)
  
  
  MedianProp_25_2  <- plot_medianprop(model_list = LM_IUCN_Species,
                                      "LatQIMed",
                                      model_data = GMPD_IUCN_Species, 
                                      fixed_lat = 25) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_40_2 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                     "LatQIMed",
                                     model_data = GMPD_IUCN_Species, 
                                     fixed_lat = 40) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_55_2 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                     "LatQIMed",
                                     model_data = GMPD_IUCN_Species, 
                                     fixed_lat = 55) +
    # theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    labs(x = "Range position, no niche position", y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "F",
             family = "Outfit", size = 5, hjust = 0)
  
  gg_y_axis <- cowplot::get_plot_component(ggplot() + 
                                             theme(text = element_text(family = "Outfit", size = 15)) + 
                                             labs(y = "Parasite prevalence"), 
                                           "ylab-l")
  gg_legend <- cowplot::get_plot_component(Latitude_main, 'guide-box-bottom', return_all = TRUE)                    
  
  
  range_plots <- 
    MedianProp_25 + MedianProp_25_2 +
    MedianProp_40 + MedianProp_40_2 +
    MedianProp_55 + MedianProp_55_2 +
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
  
  ggsave(here::here("Figures/range position.pdf"), range_plots, width = 10, height = 10, 
         device = cairo_pdf)
}

## Supp fig 1, asym ####
# Range position asymmetry
if(plot_outputs) {
  MedianProp_asym_40  <- plot_medianprop_asym(model_list = LM_IUCN_Species, 
                                              model_data = GMPD_IUCN_Species, 
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
                               model_name = "LatQIMedQIDist_5", 
                               model_data = GMPD_IUCN_Species,
                               raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "A",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  Latitude_10 <- plot_latitude(model_list = LM_IUCN_Species, 
                               model_data = GMPD_IUCN_Species) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = 0.1*90, y = 0.85,
             label = "B",
             family = "Outfit", size = 5, hjust = 0) +
    theme(legend.position = "none")
  Latitude_20 <- plot_latitude(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIDist_20", 
                               model_data = GMPD_IUCN_Species,
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
                                   model_name = "LatQIMedQIDist_5", 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 40,
                                   raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "D",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_10 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_data = GMPD_IUCN_Species, 
                                   fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "E",
             family = "Outfit", size = 5, hjust = 0)
  MedianProp_20 <- plot_medianprop(model_list = LM_IUCN_Species, 
                                   model_name = "LatQIMedQIDist_20", 
                                   model_data = GMPD_IUCN_Species, 
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
                               model_name = "LatQIMedQIDist_5", 
                               model_data = GMPD_IUCN_Species, 
                               fixed_lat = 40,
                               raster_res = "_5") +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "G",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_10 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_data = GMPD_IUCN_Species, 
                               fixed_lat = 40) +
    theme(          
      axis.text.y = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    annotate(geom = "text", x = -0.8, y = 0.85,
             label = "H",
             family = "Outfit", size = 5, hjust = 0)
  DistProp_20 <- plot_distprop(model_list = LM_IUCN_Species, 
                               model_name = "LatQIMedQIDist_20", 
                               model_data = GMPD_IUCN_Species, 
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

## Supp fig 4, random intercept diagram ####
if(plot_outputs){
  # Set seed for reproducibility
  set.seed(123)
  
  # Generate x values
  x <- seq(-0.9, 0.9, length.out = 100)
  
  # Define coefficients
  a <- 0.4
  b <- 0
  c <- 0.3
  
  # Generate random noise
  epsilon <- rnorm(100, mean = 0, sd = 0.1)
  
  # Generate y values from a quadratic relationship
  y <- a * x^2 + b * x + c + epsilon
  
  # Plot the data
  # plot(x, y, main = "Quadratic Relationship", xlab = "x", ylab = "y", pch = 16, col = "blue")
  quad_data <- data.frame(x=x,y=y)
  quad_left <- ggplot(quad_data, aes(x = x, y = y)) +
    geom_point(color = "#EDAE49") +  # scatterplot
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), 
                se = FALSE, color = "#D1495B", linewidth = 2, lineend = "round") +
    labs(title = "a",
         x = "Range position",
         y = "Parasite prevalence") +
    lims(x = c(-1,1), y = c(0,1))+
    custom_theme
  
  tick_width <- 0.1
  quad_data$x_start <- quad_data$x - tick_width
  quad_data$x_end   <- quad_data$x + tick_width 
  
  # Plot
  quad_right <- ggplot(quad_data) +
    geom_segment(aes(x = x_start, xend = x_end, y = y), 
                 color = "#222E50", linewidth = 0.7, lineend = "round") +
    geom_segment(aes(x = min(x), xend = max(x), 
                     y = mean(y)),
                 color = "#D1495B", linewidth = 2, lineend = "round") +
    labs(title = "b",
         x = "Range position",
         y = "Parasite prevalence") +
    lims(x = c(-1,1), y = c(0,1))+
    custom_theme
  
  ranef_diagram <- quad_left + quad_right +
    plot_layout(
      design = "AB")
  
  
  ggsave(here::here("Figures/ranef diagram.pdf"), ranef_diagram, width = 10, height = 5, 
         device = cairo_pdf)
  
}

## test figures ####

plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))
  
  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Host[[model_name]])$HostCorrectedName %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("HostCorrectedName")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "HostCorrectedName", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           CorrectedDistPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Host[[model_name]], 
                               plot_data, 
                               re.form = ~(1|HostCorrectedName), 
                               type = "response"))
  
  
  ggplot(plot_data, aes(x = MedianPropSquared, y = model_out, group = HostCorrectedName)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Range position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

MedianProp_host <- plot_temp(model_list = LM_IUCN_Species, 
                            model_data = GMPD_IUCN_Species, 
                            fixed_lat = 25)

ggsave(here::here("Figures/range position host.pdf"), MedianProp_host, width = 10, height = 9, 
       device = cairo_pdf)

plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))

  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Both[[model_name]])$ParasiteCorrectedName %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("ParasiteCorrectedName")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName, ParasiteCorrectedName,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "ParasiteCorrectedName", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           CorrectedDistPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Both[[model_name]], 
                                plot_data, 
                                re.form = ~(1|ParasiteCorrectedName), 
                                type = "response"))
  
  
  ggplot(plot_data, aes(x = MedianPropSquared, y = model_out, group = ParasiteCorrectedName)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Range position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

MedianProp_par <- plot_temp(model_list = LM_IUCN_Species, 
                            model_data = GMPD_IUCN_Species, 
                            fixed_lat = 25)

ggsave(here::here("Figures/range position parasite.pdf"), MedianProp_par, width = 10, height = 9, 
       device = cairo_pdf)



plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))
  
  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Cross[[model_name]])[["HostCorrectedName:ParasiteCorrectedName"]] %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("HostParasite")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName, ParasiteCorrectedName, HostParasite,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "HostParasite", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           CorrectedDistPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Cross[[model_name]], 
                               plot_data, 
                               re.form = ~(1|HostCorrectedName:ParasiteCorrectedName), 
                               type = "response"))
  
  
  ggplot(plot_data, aes(x = MedianPropSquared, y = model_out, group = HostParasite)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Range position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

MedianProp_cross <- plot_temp(model_list = LM_IUCN_Species, 
                              model_data = GMPD_IUCN_Species, 
                              fixed_lat = 25)

ggsave(here::here("Figures/range position hostparasite.pdf"), MedianProp_cross, width = 10, height = 9, 
       device = cairo_pdf)

## test figures 2 ####

plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))
  
  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Host[[model_name]])$HostCorrectedName %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("HostCorrectedName")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "HostCorrectedName", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           MedianPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Host[[model_name]], 
                               plot_data, 
                               re.form = ~(1|HostCorrectedName), 
                               type = "response"))
  
  
  ggplot(plot_data, aes(x = CorrectedDistPropSquared, y = model_out, group = HostCorrectedName)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Niche position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

DistProp_host <- plot_temp(model_list = LM_IUCN_Species, 
                           model_data = GMPD_IUCN_Species, 
                           fixed_lat = 25)

ggsave(here::here("Figures/niche position host.pdf"), DistProp_host, width = 10, height = 9, 
       device = cairo_pdf)

plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))
  
  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Both[[model_name]])$ParasiteCorrectedName %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("ParasiteCorrectedName")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName, ParasiteCorrectedName,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "ParasiteCorrectedName", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           MedianPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Both[[model_name]], 
                               plot_data, 
                               re.form = ~(1|ParasiteCorrectedName), 
                               type = "response"))
  
  
  ggplot(plot_data, aes(x = CorrectedDistPropSquared, y = model_out, group = ParasiteCorrectedName)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Niche position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

DistProp_par <- plot_temp(model_list = LM_IUCN_Species, 
                          model_data = GMPD_IUCN_Species, 
                          fixed_lat = 25)

ggsave(here::here("Figures/niche position parasite.pdf"), DistProp_par, width = 10, height = 9, 
       device = cairo_pdf)



plot_temp <- function(model_list,
                      model_name = "LatQIMedQIDist", 
                      method = "iucn",
                      model_data,
                      fixed_lat, 
                      raster_res = ""){
  # Prep fixed values
  meanLat   <- mean(abs(model_data$Latitude))
  
  model_data <- rename_with(model_data, 
                            ~ gsub(paste0("_", method), "", .),
                            ends_with(method))
  
  fixed_lat.sc <- (fixed_lat - meanLat)/(sd(abs(model_data$Latitude) - meanLat))
  
  model_coefs <- coef(model_list$Cross[[model_name]])[["HostCorrectedName:ParasiteCorrectedName"]] %>% 
    rename(Intercept = `(Intercept)`) %>% 
    rownames_to_column("HostParasite")
  
  plot_data <- model_data %>% 
    select(HostCorrectedName, ParasiteCorrectedName, HostParasite,
           LatitudeScaled, MedianPropSquared, CorrectedDistPropSquared) %>% 
    left_join(model_coefs, by = "HostParasite", suffix = c(".data", "")) %>% 
    mutate(LatitudeScaled = fixed_lat.sc,
           MedianPropSquared = 0)
  
  plot_data <- plot_data %>% 
    mutate(model_out = predict(model_list$Cross[[model_name]], 
                               plot_data, 
                               re.form = ~(1|HostCorrectedName:ParasiteCorrectedName), 
                               type = "response"))
  
  
  ggplot(plot_data, aes(x = CorrectedDistPropSquared, y = model_out, group = HostParasite)) +
    geom_line(linewidth = 1, lineend = "round", color = "#222E50") +
    custom_theme +
    labs(x = "Niche position", y = "Parasite prevalence") +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1,1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
    theme(legend.position = "none")
}

DistProp_cross <- plot_temp(model_list = LM_IUCN_Species, 
                              model_data = GMPD_IUCN_Species, 
                              fixed_lat = 25)

ggsave(here::here("Figures/niche position hostparasite.pdf"), DistProp_cross, width = 10, height = 9, 
       device = cairo_pdf)
