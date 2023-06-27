# Analysing distance metrics (without climate)
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(MCMCglmm)
library(lme4)
library(here)
library(postMCMCglmm)
library(beepr)
source(here::here("Functions.R"))

GMPD_Distances_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), header = TRUE, stringsAsFactors = FALSE)

# resources for mcmcglmm:
# https://mran.microsoft.com/snapshot/2020-12-04/web/packages/MCMCglmm/vignettes/CourseNotes.pdf
# https://tomhouslay.files.wordpress.com/2017/02/indivvar_plasticity_tutorial_mcmcglmm1.pdf
# http://www2.uaem.mx/r-mirror/web/packages/MCMCglmm/vignettes/Overview.pdf
# https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf
# https://stackoverflow.com/questions/47598123/how-do-i-extract-random-effects-from-mcmcglmm

# other resources:
# https://www.researchgate.net/publication/314239202_Avoiding_the_misuse_of_BLUP_in_behavioural_ecology

# To check autocorrelation of means and variance across iterations : 
# autocorr(Model_Name$Sol) and autocorr(Model_Name$VCV)
# TODO: Write a function that takes the model and outputs lag 10 autocorr for each variable

# To look at random effects estimates, Model_Name$Sol 

# Single fixed effect:
# Fixed effect only:                            Response ~ FixedEffect
# Fixed effect and random intercept:            Response ~ FixedEffect, random = ~ RandomEffect
# Fixed effect and random slope:                Response ~ FixedEffect, random = ~ us(FixedEffect):RandomEffect
# Fixed effect and random slope and intercept:  Response ~ FixedEffect, random = ~ us(FixedEffect + 1):RandomEffect

# Multiple fixed effects:
# Fixed effect only:                                    Response ~ FixedEffect1 + FixedEffect2
# Fixed effect and random intercept:                    Response ~ FixedEffect1 + FixedEffect2, random = ~ RandomEffect
# Fixed effect and random slope:                        Response ~ FixedEffect1 + FixedEffect2, random = ~ us(FixedEffect1 + FixedEffect2):RandomEffect
# Same as above? (above may fit interaction?:           Response ~ FixedEffect1 + FixedEffect2, random = ~ us(FixedEffect1):RandomEffect + us(FixedEffect2):RandomEffect
# Fixed effect and random slope and intercept:          Response ~ FixedEffect1 + FixedEffect2, random = ~ us(FixedEffect1 + FixedEffect2 +1):RandomEffect
# Fixed effect and partial random slope:                Response ~ FixedEffect1 + FixedEffect2, random = ~ us(FixedEffect1):RandomEffect
# Fixed effect and partial random slope and intercept:  Response ~ FixedEffect1 + FixedEffect2, random = ~ us(FixedEffect1 + 1):RandomEffect

# Nested random effects:


# summary:
# linear models of parasitism rates using range position metrics and latitude
# first done using lme4 but later converted to mcmcglmm to cope with gappy data
# naming system of models:
## Md_[iucn/gbif]_[species/subgroup]_[...]

# IUCN species #################################################################

## Prepping data ####
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

## lme4 models ####
# TODO: summary of models in this section

### 1: Lat + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible 
# with host species as a random effect
Md_IUCN_Species_1 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName),
                           data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_1)
# Intercept:  -0.53746
# Latitude:   -0.42540
# MedianProp: 0.21591

# notes:
# Both latitude and proportional distance to range median are significant.
# Could model proportional distance to range median better as it's a bit simplistic here
# There are likely differences between range hemispheres and potentially a quadratic pattern

#### visualisation model 1
# TODO: model 1 visualisation

### 1a: Lat Med :sub ####
Md_IUCN_Species_1a <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                              (1|HostCorrectedName:subgroup) + 
                              (0+LatitudeScaled + MedianPropScaled|HostCorrectedName:subgroup),
                            data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_1a)

# notes:
# accounting for subgroups accounts for genetic variation and makes the slope stronger with a lower p value
# slope of model 1 is damped by local adaptation
# parasitism rate is Not driven by subgroups and genetic variation

### 2: + AboveMedn ####
# Latitude and proportional distance to range median modelled as in 1,
# This time including above or below median as a random effect for proportional distance to median
# (It's not relevant to include it for latitude)
Md_IUCN_Species_2 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName) + 
                             (0+MedianPropScaled|AboveMedn),
                   data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_2)

# notes:
# AIC increased from 8802.7 to 8804.7 so no improvement from adding abovemedn
# Differences from being above or below the median are likely best explained by latitude
# Still good to test out a quadratic model
# Above or below median is in relation to nearest pole rather than just north pole
# I.e. in southern hemisphere, "above median" refers to the (southern) polewards half of the range

#### visualisation model 2
# TODO: model 2 visualisation

### 3: quadratic ####
# Latitude with a quadratic term for proportional distance to range median
Md_IUCN_Species_3 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                      data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_3)

# notes:
# AIC is 8774.7 so definitely no improvement again, and the quadractic term is non-significant
# Should stick with the simplest, model 1
# I previously tried the same with only poly 1 in the random effects but AIC was comparable to m2

#### visualisation model 3
# TODO: model 3 visualisation

### 4: raw quadratic ####
# The same as model 3 but with raw rather than the orthogonal polynomial
Md_IUCN_Species_4 <- glmer(formula  = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE)|HostCorrectedName),
                           data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_4)

# notes:
# this one was just for fun so I could see the difference in estimates between poly and raw 
# Model fails to converge

### 5: scaled and squared ####
# Essentially the same as model 3 except proportional distance to median has been scaled
Md_IUCN_Species_5 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquScaled, 2) +
                             (1|HostCorrectedName) +
                             (0+LatitudeScaled + poly(MedianPropSquScaled, 2)|HostCorrectedName),
                           data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_5)

# notes:
# model fails to converge but gives basically the same output as model 3

### 6: interaction ####
Md_IUCN_Species_6a <- glmer(formula = Prevalence ~ LatitudeScaled*MedianPropScaled + 
                              (1|HostCorrectedName) + 
                              (0+LatitudeScaled + MedianPropScaled|HostCorrectedName),
                            data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_6a)

Md_IUCN_Species_6b <- glmer(formula = Prevalence ~ LatitudeScaled*MedianPropScaled + 
                             (1|HostCorrectedName:subgroup) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName:subgroup),
                           data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_6b)

Md_IUCN_Species_6c <- glm(formula = Prevalence ~ LatitudeScaled*MedianPropScaled,
                            data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_6c)


# Interaction term is not significant anywhere, but it is in the mcmc (4b) version
# But 4b has a very small posterior mean for the interaction and DIC is high in comparison to alternatives

## MCMC models ####

### Latitude practice ####
#### 1: Latitude fixed only ####
Mc_IUCN_Species_L <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                              data    = GMPD_IUCN_Species,
                              family  = "multinomial2",
                              nitt    = 130000,
                              thin    = 100,
                              burnin  = 30000)
summary(Mc_IUCN_Species_L)
# DIC:        561302.3
# Intercept:  -0.5790 
# Latitude:   -0.4827

autocorr(Mc_IUCN_Species_L$Sol)
autocorr(Mc_IUCN_Species_L$VCV)
# from the course notes page 22 we're looking at the second row and want it to be below 0.1
# for both the mean (sol) and variance (vcv)
# these ones are fine

#### 2a: Lat random intercept ####
Mc_IUCN_Species_LHa <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               random = ~ HostCorrectedName,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               nitt   = 130000,
                               thin   = 100,
                               burnin = 30000,
                               pr     = TRUE)
summary(Mc_IUCN_Species_LHa)
# DIC:        561138.2
# Intercept:  -0.1444
# Latitude:   -0.2114

autocorr(Mc_IUCN_Species_LHa$Sol)
autocorr(Mc_IUCN_Species_LHa$VCV)
# fine for the means, but not the variance in random effects
# TODO: run again with more iterations

#### 2b: Lat random slope ####
Mc_IUCN_Species_LHb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                               random = ~ us(LatitudeScaled):HostCorrectedName,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               pr = TRUE)
summary(Mc_IUCN_Species_LHb)
## DIC:        561178.2
# Intercept:  -0.6627
# Latitude:   -0.7245

autocorr(Mc_IUCN_Species_LH$Sol)
autocorr(Mc_IUCN_Species_LH$VCV)
# fine for the means, but not the variance in random effects
# Should probably ignore this one because it's random slope without random intercept

#### 3: Latitude host ####
Mc_IUCN_Species_L_H <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled,
                                random  = ~ us(1 + LatitudeScaled):HostCorrectedName,
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000,
                                pr      = TRUE)
summary(Mc_IUCN_Species_L_H)
## DIC:        561097.4
# Intercept:  -0.13652
# Latitude:   -0.29301

autocorr(Mc_IUCN_Species_L_H$Sol)
autocorr(Mc_IUCN_Species_L_H$VCV)
# again fine for means but not for variance 
# TODO: run with more iterations

# not an interesting model, just used to get to grips with random effects in mcmc

### G-structure ####
#### 4a: Lat + Med fixed ####
Mc_IUCN_Species_LMa <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000)
summary(Mc_IUCN_Species_LMa)
# DIC:        561268.3
# Intercept:  -0.5747
# Latitude:   -0.4492
# MedianProp: 0.3850

autocorr(Mc_IUCN_Species_LMa$Sol)
autocorr(Mc_IUCN_Species_LMa$VCV)
# all fine

#### 4b: Lat + Med : fixed ####
Mc_IUCN_Species_LMb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled*MedianPropScaled,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               nitt   = 130000,
                               thin   = 100,
                               burnin = 30000)
summary(Mc_IUCN_Species_LMb)
# DIC:        561268.3
# Intercept:  -0.5747
# Latitude:   -0.4492
# MedianProp: 0.3850

autocorr(Mc_IUCN_Species_LMb$Sol)
autocorr(Mc_IUCN_Species_LMb$VCV)
# all fine

# mentioned above (Md 6) that the interaction may be significant here but it has a very small estimate and the DIC is high

#### 5: Lat + Med rand intercept####
Mc_IUCN_Species_LMH <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               random = ~ HostCorrectedName,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               nitt   = 130000,
                               thin   = 100,
                               burnin = 30000,
                               pr     = TRUE)
summary(Mc_IUCN_Species_LMH)
# DIC:        561146.4
# Intercept:  -0.22402
# Latitude:   -0.28840
# MedianProp: 0.12705

autocorr(Mc_IUCN_Species_LMH$Sol)
autocorr(Mc_IUCN_Species_LMH$VCV)
# A bit high for host variance
# TODO: iterations

#### 6a: Lat + Med rand : slope ####
Mc_IUCN_Species_LM_Ha <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                random  = ~ us(LatitudeScaled + MedianPropScaled):HostCorrectedName,
                                data    = GMPD_IUCN_Species,
                                family  = "multinomial2",
                                nitt    = 130000,
                                thin    = 100,
                                burnin  = 30000,
                                pr      = TRUE)
summary(Mc_IUCN_Species_LM_Ha)
# from "Avoiding the misuse of BLUP..."
# "The labels correspond to how we would set up a covariance matrix for the intercept and slope variance" p15
# Latitide:Latitude is the variance of latitude
# MedianProp:Latitude is the covariance between medianprop and latitude
# So we don't want this random effect specification (or G-structure) because 
# we're not (yet?) interested in the covariance between latitude and medianprop

# DIC:        561144.6
# Intercept:  -0.7016
# Latitude:   -0.6076 
# MedianProp: 0.3050 
autocorr(Mc_IUCN_Species_LM_Ha$Sol)
autocorr(Mc_IUCN_Species_LM_Ha$VCV)
# same
# TODO: iterations

#### 6b: Lat + Med rand slope ####
Mc_IUCN_Species_LM_Hb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                 random = ~ us(LatitudeScaled):HostCorrectedName + us(MedianPropScaled):HostCorrectedName,
                                 data   = GMPD_IUCN_Species,
                                 family = "multinomial2",
                                 nitt   = 130000,
                                 thin   = 100,
                                 burnin = 30000,
                                 pr     = TRUE)
summary(Mc_IUCN_Species_LM_Hb)
# This G-structure is more appropriate because there is no covariance between Latitude and MedianProp

# DIC:        561127 <3
# Intercept:  -0.7042
# Latitude:   -0.6146 
# MedianProp: 0.3034 
autocorr(Mc_IUCN_Species_LM_Hb$Sol)
autocorr(Mc_IUCN_Species_LM_Hb$VCV)
# TODO: iterations

#### 7a: Lat + Med random :)* ####
Mc_IUCN_Species_L_M_Ha <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                  random  = ~ HostCorrectedName + us(LatitudeScaled):HostCorrectedName + us(MedianPropScaled):HostCorrectedName,
                                  data    = GMPD_IUCN_Species,
                                  family  = "multinomial2",
                                  nitt    = 130000,
                                  thin    = 100,
                                  burnin  = 30000,
                                  pr      = TRUE)
summary(Mc_IUCN_Species_L_M_Ha)
# DIC:        561091 
# Intercept:  -0.39864 **
# Latitude:   -0.43493 ** 
# MedianProp: 0.20899 * 
autocorr(Mc_IUCN_Species_L_M_Ha$Sol)
autocorr(Mc_IUCN_Species_L_M_Ha$VCV)
# TODO: iterations

# https://stats.stackexchange.com/questions/86958/variance-covariance-structure-for-random-effects-in-lme4
# "a model with random intercept and slope where the intercept and slope are uncorrelated"


#### 7b: Lat + Med random :)* ####
Mc_IUCN_Species_L_M_Hb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                  random  = ~ us(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName,
                                  data    = GMPD_IUCN_Species,
                                  family  = "multinomial2",
                                  nitt    = 130000,
                                  thin    = 100,
                                  burnin  = 30000,
                                  pr      = TRUE)
summary(Mc_IUCN_Species_L_M_Hb)
# DIC:        561086.3 
# Intercept:  -0.445467 **
# Latitude:   -0.521344 ***
# MedianProp: 0.205312 .
autocorr(Mc_IUCN_Species_L_M_Hb$Sol)
autocorr(Mc_IUCN_Species_L_M_Hb$VCV)
# TODO: iterations

# Back to being uncertain on whether I'm specifying the random effects correctly
# Should I allow them to covary?
# Should I try a different covariance matrix (idh)?
# us() means the covariance is unstructured (i.e. I haven't imposed any constraints on the variance/covariance so it has to be estimated)
# I could impose fixed variance 
# https://www.theanalysisfactor.com/unstructured-covariance-matrix-when-it-does-and-doesn%E2%80%99t-work/

# I think what is happening is that the covariance is estimated in the G-matrix and it's using up df
# "The only other common structure for a G matrix is a variance components structure, 
# "which fits different variance estimates, but 0 covariances"
# https://stackoverflow.com/questions/62650072/what-do-the-us-idh-idv-stand-for-in-the-mcmcglmm-in-r
# next thing to do is fit a model with idh which fits variance but not covariance and see if it's the same as 7a

#### 7c: Lat + Med random :)* ####
Mc_IUCN_Species_L_M_Hc <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random = ~ idh(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName,
                                   data   = GMPD_IUCN_Species,
                                   family = "multinomial2",
                                   nitt   = 130000,
                                   thin   = 100,
                                   burnin = 30000,
                                   pr     = TRUE)
summary(Mc_IUCN_Species_L_M_Hc)
# DIC:        561090.2 <3
# Intercept:  -0.3949
# Latitude:   -0.4310 
# MedianProp: 0.2055 
autocorr(Mc_IUCN_Species_L_M_Hc$Sol)
autocorr(Mc_IUCN_Species_L_M_Hc$VCV)
# TODO: iterations

# it is the same as 7a and the DIC of us is lowest
# but what's my justification of using idh? Why wouldn't I expect covariance between the intercept and slopes?
# 

### Separate median? ####

#### 8a: Lat + Med rand Med ####
Mc_IUCN_Species_L_MHa <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                 random = ~ us(MedianPropScaled):HostCorrectedName,
                                 data   = GMPD_IUCN_Species,
                                 family = "multinomial2",
                                 nitt   = 130000,
                                 thin   = 100,
                                 burnin = 30000,
                                 pr     = TRUE)
summary(Mc_IUCN_Species_L_MHa)
# DIC:        561209.9
# Intercept:  -0.6290
# Latitude:   -0.4114 
# MedianProp: 0.3780 
autocorr(Mc_IUCN_Species_L_MHa$Sol)
autocorr(Mc_IUCN_Species_L_MHa$VCV)
# TODO: iterations

#### 8b: Lat + Med rand Med int ####
Mc_IUCN_Species_L_MHb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                 random = ~ HostCorrectedName + us(MedianPropScaled):HostCorrectedName,
                                 data   = GMPD_IUCN_Species,
                                 family = "multinomial2",
                                 nitt   = 130000,
                                 thin   = 100,
                                 burnin = 30000,
                                 pr     = TRUE)
summary(Mc_IUCN_Species_L_MHb)
# DIC:        561111 
# Intercept:  -0.36787
# Latitude:   -0.30979 
# MedianProp: 0.21061
autocorr(Mc_IUCN_Species_L_MHb$Sol)
autocorr(Mc_IUCN_Species_L_MHb$VCV)
# TODO: iterations

#### 8c: Lat + Med rand Med full ####
Mc_IUCN_Species_L_MHc <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                  random  = ~ us(1+MedianPropScaled):HostCorrectedName,
                                  data    = GMPD_IUCN_Species,
                                  family  = "multinomial2",
                                  nitt    = 130000,
                                  thin    = 100,
                                  burnin  = 30000,
                                  pr      = TRUE)
summary(Mc_IUCN_Species_L_MHc)
# DIC:        561111 
# Intercept:  -0.36787
# Latitude:   -0.30979 
# MedianProp: 0.21061
autocorr(Mc_IUCN_Species_L_MHc$Sol)
autocorr(Mc_IUCN_Species_L_MHc$VCV)
# TODO: iterations

# same as for 7a:c allowing covariance between random effects takes away the effect of median prop quite drastically
# I think that in the frequentist version variance or covariance or something is fixed
# I can't remember exactly but might be worth looking up if it explains the differences between mcmclmm and glmer

### Subgroups ####

#### 9a: subgroup intercept ####
Mc_IUCN_Species_LM_HSa <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random = ~ HostSubgroup,
                                   data   = GMPD_IUCN_Species,
                                   family = "multinomial2",
                                   nitt   = 130000,
                                   thin   = 100,
                                   burnin = 30000,
                                   pr     = TRUE)
summary(Mc_IUCN_Species_LM_HSa)

#### 9b: subgroup full ####
Mc_IUCN_Species_LM_HSb <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random = ~ us(LatitudeScaled + MedianPropScaled + 1):HostSubgroup,
                                   data   = GMPD_IUCN_Species,
                                   family = "multinomial2",
                                   nitt   = 130000,
                                   thin   = 100,
                                   burnin = 30000,
                                   pr     = TRUE)
summary(Mc_IUCN_Species_LM_HSb)

# One run gave the error "ill-conditioned G/R structure"
# Have a look at page 23/24 of the course notes
# This indicates a reducible chain (in that run) which can arise from variance being zero
# I should probably define uninformative priors

#### 9b: subgroup idh ####
Mc_IUCN_Species_LM_HSc <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                                   random = ~ idh(LatitudeScaled + MedianPropScaled + 1):(HostSubgroup),
                                   data   = GMPD_IUCN_Species,
                                   family = "multinomial2",
                                   nitt   = 130000,
                                   thin   = 100,
                                   burnin = 30000,
                                   pr     = TRUE)
summary(Mc_IUCN_Species_LM_HSc)


# notes:
# DIC of model we selected previously using glmer is the best again :)
# and estimates are very similar

# The first three (latitude only) models tell me that host isn't super important for latitude effects
# because the dic was lowest for the one with only host intercept as random, not slope
# Intercept being different between hosts makes sense because they will each have a base level difference in parasitism rates
# Adding hosts in for slopes takes up a lot of variance so they have to be informative to reduce the DIC
# 

# Reading "Avoiding the misuse of BLUPs in behavioural ecology"
# Would I assume that my host/subgroup slopes are independent?
# How would I model phylogeny of subgroups with subgroup as a random effect?
# Presumably similarly to how I modelled subgroup in 9, though I would have to modify the phylogeny to include subgroups
# They would all have to have same length branches because I have no knowledge of their divergence 

#### 9c: subgroup half full ####
host_and_subgroup3 <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                               random = ~ idh(LatitudeScaled + MedianPropScaled + 1):HostCorrectedName + idh(MedianPropScaled + 1):HostSubgroup,
                               data   = GMPD_IUCN_Species,
                               family = "multinomial2",
                               nitt   = 130000,
                               thin   = 100,
                               burnin = 30000,
                               pr     = TRUE)
summary(host_and_subgroup3)


#########

Mc_IUCN_Species_0 <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                              data = GMPD_IUCN_Species,
                              family = "multinomial2",
                              verbose = FALSE)

Mc_IUCN_Species_1 <- MCMCglmm(cbind(ParasiteDetected, ParasiteUndetected) ~ LatitudeScaled + MedianPropScaled,
                              random = ~ HostCorrectedName,
                              data = GMPD_IUCN_Species,
                              family = "multinomial2",
                              verbose = FALSE)

# IUCN subgroups ###############################################################

## Prepping data ####
GMPD_IUCN_Subgroup <- GMPD_Distances_Data %>%
  dplyr::filter(RestrSub, 
                RangeMethod == "iucn",
                RangeTaxonLvl == "subgroup") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected)

## lme4 models ####
# TODO: summary of models in this section

### 1: Latitude + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible 
# with host species as a random effect and subgroup nested within host species
Md_IUCN_Subgroup_1 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName:subgroup) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName:subgroup),
                           data = GMPD_IUCN_Subgroup, family = binomial)

summary(Md_IUCN_Subgroup_1)

# notes:
# In comparison to the species level analysis, quite different results
# Basic conclusion is that subgroup level analysis is irrelevant overall
# Though there may be some species for which subspecies is relevant
# This is backed up by the gbif verison of this model 
# May still be worth investigating a little further:
# TODO: Look at sample size and range size effects for subspecies

# subspecies groups absorb variation from median position (genetic? coincidental?)
# variance used up by random or fixed (posthoc)
# what would it mean to add

#### visualisation model 1
# TODO: model 1 visualisation

# GBIF species #################################################################

## Prepping data ####
GMPD_GBIF_Species <- GMPD_Distances_Data %>%
  dplyr::filter(CleanAll, 
                RangeMethod == "gbif",
                RangeTaxonLvl == "species",
                MedianProp <= 1) %>% # temporary fix # TODO: tidy
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected)

## lme4 models ####
# TODO: summary of models in this section

### 1: Latitude + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible 
# with host species as a random effect
Md_GBIF_Species_1 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_1)

# notes:
# Latitude is significant with a similar estimate
# but proportional distance to median is not significant. 
# Maybe because of sample sizes of gbif species
# TODO: redo with limited sample sizes

#### visualisation model 1
# TODO: model 1 visualisation

### 1: Latitude + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible 
# with host species as a random effect
Md_GBIF_Species_1a <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName:subgroup) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName:subgroup),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_1a)

# notes:
# GBIF subgroups ###############################################################

## Prepping data ####
GMPD_GBIF_Subgroup <- GMPD_Distances_Data %>%
  dplyr::filter(CleanSub, 
                RangeMethod == "gbif",
                RangeTaxonLvl == "subgroup",
                MedianProp <= 1) %>% # temporary fix # TODO: tidy
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared),
         ParasiteDetected       = as.integer(round(HostsSampled * Prevalence, 0)),
         ParasiteUndetected     = HostsSampled - ParasiteDetected)

## lme4 models ####
# TODO: summary of models in this section

### 1: Latitude + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible with host species as a random effect
Md_GBIF_Subgroup_1 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName:subgroup) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName:subgroup),
                           data = GMPD_GBIF_Subgroup, family = binomial)

summary(Md_GBIF_Subgroup_1)

# notes:
# singular fit
# Median prop is closer to being significant than in the iucn equivalent
# Talked about relevance of subgrouping in this context
## Could be because of data coverage among subgroups being quite variable
## But that correlates witht the availability of data in gmpd so likely not this because it wouldn't hold much sway in the model

## Could be because of positioning of subgroups: those that are latitudinally aligned will not differ from species much in their estimates of medianprop
## e.g. Alces alces or Canis aureus
## in comparison to ones that are latitudinally scattered
## e.g. Cervus nippon

## Or could be something else that we haven't thought of

# notes:
# Both latitude and proportional distance to range median are significant.
# Could model proportional distance to range median better as it's a bit simplistic here
# There are likely differences between range hemispheres and potentially a quadratic pattern

#### visualisation model 1
# TODO: model 1 visualisation

### 2: + AboveMedn ####
# Latitude and proportional distance to range median modelled as in 1,
# This time including above or below median as a random effect for proportional distance to median
# (It's not relevant to include it for latitude)
Md_GBIF_Subgroup_2 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName) + 
                             (0+MedianPropScaled|AboveMedn),
                           data = GMPD_GBIF_Subgroup, family = binomial)

summary(Md_GBIF_Subgroup_2)

# notes:
# AIC increased from 8802.7 to 8804.7 so no improvement from adding abovemedn
# Differences from being above or below the median are likely best explained by latitude
# Still good to test out a quadratic model
# Above or below median is in relation to nearest pole rather than just north pole
# I.e. in southern hemisphere, "above median" refers to the (southern) polewards half of the range

#### visualisation model 2
# TODO: model 2 visualisation

### 3: quadratic ####
# Latitude with a quadratic term for proportional distance to range median
Md_GBIF_Subgroup_3 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                           data = GMPD_GBIF_Subgroup, family = binomial)

summary(Md_GBIF_Subgroup_3)

# notes:
# AIC is 8781.6 so definitely no improvement again, and the quadractic term is non-significant
# Should stick with the simplest, model 1
# I previously tried the same with only poly 1 in the random effects but AIC was comparable to m2

#### visualisation model 3
# TODO: model 3 visualisation

### 4: raw quadratic ####
# The same as model 3 but with raw rather than the orthogonal polynomial
Md_GBIF_Subgroup_4 <- glmer(formula  = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE)|HostCorrectedName),
                           data = GMPD_GBIF_Subgroup, family = binomial)

summary(Md_GBIF_Subgroup_4)

# notes:
# this one was just for fun so I could see the difference in estimates between poly and raw 
# Model fails to converge

### 5: scaled and squared ####
# Essentially the same as model 3 except proportional distance to median has been scaled
Md_GBIF_Subgroup_5 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquScaled, 2) +
                             (1|HostCorrectedName) +
                             (0+LatitudeScaled + poly(MedianPropSquScaled, 2)|HostCorrectedName),
                           data = GMPD_GBIF_Subgroup, family = binomial)

summary(Md_GBIF_Subgroup_5)

# notes:
# model fails to converge but gives basically the same output as model 3

#old############################################################################

# sorting out data scaling etc. ####

# using data restricted by iucn species polygons
# scale variables 
medDatTemp <- Distances_Data_res_all[["DistanceMetrics"]]
medDatTemp$LatitudeScaled <- base::scale(abs(medDatTemp$Latitude))
medDatTemp$EquatorwardsProp.sc <- base::scale(medDatTemp$EquatorwardsProp)
medDatTemp$MedianPropScaled <- base::scale(medDatTemp$MedianProp)
medDatTemp$MedianPropSquared <- medDatTemp$MedianProp 
medDatTemp[!medDatTemp$AboveMedn, "MedianPropSquared"] <- medDatTemp[!medDatTemp$AboveMedn, "MedianPropSquared"] * -1
medDatTemp$MedianProp.sqsc <- base::scale(medDatTemp$MedianPropSquared)
#medDatTemp$MedianPropSquared <- medDatTemp$MedianPropScaled * medDatTemp$MedianPropScaled




# running the same again but only including species with above 10 samples
medDatTemp10 <- medDatTemp %>% 
  group_by(HostCorrectedName) %>%
  filter(n() > 10)
medDatTemp10$LatitudeScaled <- base::scale(abs(medDatTemp10$Latitude))
medDatTemp10$EquatorwardsProp.sc <- base::scale(medDatTemp10$EquatorwardsProp)
medDatTemp10$MedianPropScaled <- base::scale(medDatTemp10$MedianProp)
medDatTemp10$MedianPropSquared <- medDatTemp10$MedianProp 
medDatTemp10[!medDatTemp10$AboveMedn, "MedianPropSquared"] <- medDatTemp10[!medDatTemp10$AboveMedn, "MedianPropSquared"] * -1
medDatTemp10$MedianProp.sqsc <- base::scale(medDatTemp10$MedianPropSquared)

medModelquad10p <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + (1|HostCorrectedName) + (0+LatitudeScaled + poly(MedianPropSquared, 1)|HostCorrectedName),
                        data = medDatTemp10, family = binomial)

medModelquad10r <- glmer(formula  = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE) + (1|HostCorrectedName) + (0+LatitudeScaled + MedianPropSquared|HostCorrectedName),
                         data = medDatTemp10, family = binomial)

medModelquad10ra <- glmer(formula  = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE) + (1|HostCorrectedName) + (0+LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE)|HostCorrectedName),
                          data = medDatTemp10, family = binomial)

## model plots without aboveMedn ####
summary(medModel)
medPlotTemp <- medDatTemp
medPlotTemp$fit <- fitted(medModel)

medIntercept <- fixef(medModel)[[1]]
medLatitude <- fixef(medModel)[[2]]
medMedian <- fixef(medModel)[[3]]

# bigline
latrange <- seq(min(abs(medDatTemp$Latitude)), max(abs(medDatTemp$Latitude)), length.out = 2000)
medrange <- seq(min(medDatTemp$MedianProp), max(medDatTemp$MedianProp), length.out = 2000)

latrange.sc <- (latrange - 
                  mean(abs(medDatTemp$Latitude)))/(sd(abs(medDatTemp$Latitude) - mean(abs(medDatTemp$Latitude))))
medrange.sc <- (medrange - 
                  mean(medDatTemp$MedianProp))/(sd(medDatTemp$MedianProp - mean(medDatTemp$MedianProp)))

bigline <- data.frame("latrange" = latrange, "medrange" = medrange, "latrange.sc" = latrange.sc, "medrange.sc" = medrange.sc)

logitlat <- medIntercept + (mean(medDatTemp$MedianPropScaled)*medMedian) + (bigline$latrange.sc*medLatitude)
logitpol <- medIntercept + (mean(medDatTemp$LatitudeScaled)*medLatitude) + (bigline$medrange.sc*medMedian)
bigline$linelat <- 1/(1+exp(-logitlat))
bigline$linepol <- 1/(1+exp(-logitpol))

# small lines
plotranef <- ranef(medModel)[["HostCorrectedName"]]
medPlotTemp$smallatlines <- NA
medPlotTemp$smalmedlines <- NA

for(i in 1:length(plotranef[[1]])){
  p.rows <- which(medPlotTemp$HostCorrectedName == rownames(plotranef)[[i]])
  c.meanmed <- mean(medPlotTemp[p.rows, "MedianPropScaled"])
  c.meanlat <- mean(abs(medPlotTemp[p.rows,"LatitudeScaled"]))
  c.logitlat <- medIntercept + plotranef[[1]][[i]]+ (c.meanmed*(medMedian + plotranef[[3]][[i]])) + (medPlotTemp[p.rows,"LatitudeScaled"]*(medLatitude + plotranef[[2]][[i]]))
  c.logitmed <- medIntercept + plotranef[[1]][[i]]+ (c.meanlat*(medLatitude + plotranef[[2]][[i]])) + (medPlotTemp[p.rows,"MedianPropScaled"]*(medMedian + plotranef[[3]][[i]]))
  
  medPlotTemp[p.rows,"smallatlines"] <- 1/(1 + exp(-c.logitlat))
  medPlotTemp[p.rows,"smalmedlines"] <- 1/(1 + exp(-c.logitmed))
}

ggplot(medPlotTemp, aes(MedianProp, Prevalence), col = HostCorrectedName) +
  geom_line(aes(y = smalmedlines, col = HostCorrectedName), size = 0.7) +
  geom_line(data = bigline, aes(y = linepol, x = medrange), size = 1) +
  #geom_point(aes(col = HostCorrectedName), alpha = 0.3) +
  labs(title = "Effect of position within range at mean latitude", subtitle = "Random host intercept and slope", x = "Distance to median", y = "Parasite prevalence")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")

## both directions, two plot method ####
bigline$rvrsmedrange <- bigline$medrange*-1

smallinespos <- medPlotTemp[medPlotTemp$AboveMedn, ]
smallinesneg <- medPlotTemp[!medPlotTemp$AboveMedn, ]
smallinesneg$MedianProp <- smallinesneg$MedianProp*-1

neg <- ggplot(smallinesneg, aes(MedianProp, Prevalence), col = HostCorrectedName) +
  geom_line(aes(y = smalmedlines, col = HostCorrectedName), size = 0.7) +
  geom_line(data = bigline, aes(y = linepol, x = rvrsmedrange), size = 1.2) +
  labs(title = "Effect of position within range at mean latitude", subtitle = "Random host intercept and slope", x = "Distance to median", y = "Parasite prevalence")+
  theme_bw() + ylim(0,1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")

pos <- ggplot(smallinespos, aes(MedianProp, Prevalence), col = HostCorrectedName) +
  geom_line(aes(y = smalmedlines, col = HostCorrectedName), size = 0.7) +
  geom_line(data = bigline, aes(y = linepol, x = medrange), size = 1.2) +
  labs(title = "", subtitle = "", x = "", y = "")+
  theme_bw() + ylim(0,1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")

## both directions, single plot method ####
biglineboth <- data.frame(medrange = c(bigline$medrange, bigline$medrange * -1),
                          linepol = rep(bigline$linepol, times = 2))

smallinesboth <- medPlotTemp %>%
  mutate(HostCorrectedName = case_when(AboveMedn ~ paste0(HostCorrectedName, "1"),
                                       !AboveMedn ~ HostCorrectedName)) %>%
  mutate(MedianProp = case_when(AboveMedn ~ MedianProp,
                                !AboveMedn ~ MedianProp * -1))

ggplot(smallinesboth, aes(MedianProp, Prevalence), col = HostCorrectedName) +
  geom_line(aes(y = smalmedlines, col = HostCorrectedName), size = 0.71) +
  geom_line(data = biglineboth, aes(y = linepol, x = medrange), size = 1.2) +
  labs(title = "Effect of position within range at mean latitude", subtitle = "Random host intercept and slope", x = "Median", y = "Parasite prevalence")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")



## quadratic plot ####
# model values
summary(medModelquad2)
medPlotTemp <- medDatTemp
medPlotTemp$fit <- fitted(medModelquad2)

medIntercept <- fixef(medModelquad2)[[1]]
medLatitude <- fixef(medModelquad2)[[2]]
medMedian <- fixef(medModelquad2)[[3]]
medMedian2 <- fixef(medModelquad2)[[4]]

polyMedian <- poly(medDatTemp$MedianPropScaled, 2, raw = TRUE)

# bigline
latrange <- seq(min(abs(medDatTemp$Latitude)), max(abs(medDatTemp$Latitude)), length.out = 2000)
medrange <- seq(min(medDatTemp$MedianProp), max(medDatTemp$MedianProp), length.out = 2000)

latrange.sc <- (latrange - mean(abs(medDatTemp$Latitude)))/(sd(abs(medDatTemp$Latitude) - mean(abs(medDatTemp$Latitude))))
medrange.sc <- (medrange - mean(medDatTemp$MedianProp))/(sd(medDatTemp$MedianProp - mean(medDatTemp$MedianProp)))

bigline <- data.frame("latrange" = latrange, "medrange" = medrange, "latrange.sc" = latrange.sc, "medrange.sc" = medrange.sc)

logitlat <- medIntercept + 
  (mean(polyMedian[,1])*medMedian) + (mean(polyMedian[,2])*medMedian2) + 
  (bigline$latrange.sc*medLatitude)

logitpol <- medIntercept + 
  (mean(medDatTemp$LatitudeScaled)*medLatitude) + 
  (poly(bigline$medrange.sc, 2, coefs = attr(polyMedian, "coefs"))[,1]*medMedian) + (poly(bigline$medrange.sc, 2, coefs = attr(polyMedian, "coefs"))[,2]*medMedian2)

bigline$linelat <- 1/(1+exp(-logitlat))
bigline$linepol <- 1/(1+exp(-logitpol))

# small lines
plotranef <- ranef(medModelquad2)[["HostCorrectedName"]]
medPlotTemp$smallatlines <- NA
medPlotTemp$smalmedlines <- NA
medPlotTemp$medproppoly <- poly(medPlotTemp$MedianPropScaled, 1)

for(i in 1:length(plotranef[[1]])){
  p.rows <- which(medPlotTemp$HostCorrectedName == rownames(plotranef)[[i]])
  c.meanmed <- mean(medPlotTemp[p.rows, "medproppoly"])
  c.meanlat <- mean(abs(medPlotTemp[p.rows,"LatitudeScaled"]))
  c.logitlat <- medIntercept + plotranef[[1]][[i]]+ (c.meanmed*(medMedian + plotranef[[3]][[i]])) + (medPlotTemp[p.rows,"LatitudeScaled"]*(medLatitude + plotranef[[2]][[i]]))
  c.logitmed <- medIntercept + plotranef[[1]][[i]]+ (c.meanlat*(medLatitude + plotranef[[2]][[i]])) + (medPlotTemp[p.rows,"medproppoly"]*(medMedian + plotranef[[3]][[i]]))
  
  medPlotTemp[p.rows,"smallatlines"] <- 1/(1 + exp(-c.logitlat))
  medPlotTemp[p.rows,"smalmedlines"] <- 1/(1 + exp(-c.logitmed))
}

ggplot(medPlotTemp, aes(MedianProp, Prevalence), col = HostCorrectedName) +
  geom_line(aes(y = smalmedlines, col = HostCorrectedName), size = 0.7) +
  geom_line(data = bigline, aes(y = linepol, x = medrange), size = 1) +
  #geom_point(aes(col = HostCorrectedName), alpha = 0.3) +
  labs(title = "Effect of position within range at mean latitude", subtitle = "Random host intercept and slope", x = "Distance to median", y = "Parasite prevalence")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")
