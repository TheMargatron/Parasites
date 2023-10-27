# Frequentist stats
# Written by Margaret Bolton mb804(at)exeter.ac.uk 

# libraries and data ###########################################################
library(lme4)
library(here)

# Model naming convention ######################################################
## Model main lists:
# LM = Frequentist rather than Bayesian
# IUCN/GBIF = method used for restricting data and calculating position within range
# Species/Subgroup = position within range calculated from position within species (SP) or subgroup (SG) range
LM_IUCN_Species  <- list()
LM_IUCN_Subgroup <- list()
LM_GBIF_Species  <- list()
LM_GBIF_Subgroup <- list()

## Model sub lists
# Fixed     = Fixed effects only
# Fixed_T   = Fixed effects using ParType data
# Host      = Host random effects
# Host_T    = Host random effects using ParType data
# ParType   = ParType random effects
# Both      = Host and ParType random effects
# Phylogeny = Host with phylogeny
LM_IUCN_Species$Fixed     <- list()
# LM_IUCN_Species$Fixed_T   <- list()
LM_IUCN_Species$Host      <- list()
# LM_IUCN_Species$Host_T    <- list()
LM_IUCN_Species$ParType   <- list()
# LM_IUCN_Species$ParPhylum <- list()
LM_IUCN_Species$Both      <- list()
LM_IUCN_Species$Phylogeny <- list()

LM_IUCN_Subgroup$Fixed     <- list()
LM_IUCN_Subgroup$Host      <- list()
LM_IUCN_Subgroup$ParType   <- list()
LM_IUCN_Subgroup$ParPhylum <- list()
LM_IUCN_Subgroup$Both      <- list()
LM_IUCN_Subgroup$Phylogeny <- list()
  
LM_GBIF_Species$Fixed     <- list()
LM_GBIF_Species$Host      <- list()
LM_GBIF_Species$ParType   <- list()
LM_GBIF_Species$ParPhylum <- list()
LM_GBIF_Species$Both      <- list()
LM_GBIF_Species$Phylogeny <- list()
  
LM_GBIF_Subgroup$Fixed     <- list()
LM_GBIF_Subgroup$Host      <- list()
LM_GBIF_Subgroup$ParType   <- list()
LM_GBIF_Subgroup$ParPhylum <- list()
LM_GBIF_Subgroup$Both      <- list()
LM_GBIF_Subgroup$Phylogeny <- list()

## Model names:
# 0 = Latitude 
# 1 = Latitude + MedianProp
# 2 = Latitude * MedianProp
# 3 = Latitude + poly(MedianPropSquared)
# 4 = Latitude + poly(MedianPropSquared):AboveMedn

# Letter order for random effects:
# no letter = estimated full vcv matrix for host
# a = No covariance estimated, only variance in host
# b = Misc alterations to host
# c = estimated full vcv matrix for host:subgroup
# d = No covariance estimated, only variance in host:subgroup
# e = Misc alterations to host:subgroup

# IUCN with Species ############################################################
## Fixed effects exploration ###################################################
### Geographic niche ###########################################################
# Looking at fixed effects first because there are lots of random effects that absorb a lot of variance and limit flexibility
#_______________________________________________________________________________
LM_IUCN_Species$Fixed$Lat <- glm(formula = Prevalence ~ LatitudeScaled,
                                       data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Lat)
# For comparisons later on
LM_IUCN_Species$Fixed$Lat_Weighted <- glm(formula = Prevalence ~ LatitudeScaled,
                                          data = GMPD_IUCN_Species, family = binomial,
                                          weights = HostsSampled)

summary(LM_IUCN_Species$Fixed$Lat_Weighted)
# For comparisons later on

LM_IUCN_Species$Fixed$Med <- glm(formula = Prevalence ~ MedianPropScaled,
                                       data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Med)
# For comparisons later on

LM_IUCN_Species$Fixed$QMed <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2),
                                 data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QMed)
# For comparisons later on

LM_IUCN_Species$Fixed$QMedAsym <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2):AboveMedn,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QMedAsym)
# For comparisons later on

LM_IUCN_Species$Fixed$LatMed <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled,
                                          data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$Lat)
# Accept Latitude and MedianProp

LM_IUCN_Species$Fixed$LatMedAsym <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled:AboveMedn,
                                        data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatMedAsym)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatMedAsym)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatMedAsym)
# Accept Latitude and MedianProp

#_______________________________________________________________________________
LM_IUCN_Species$Fixed$Quad <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                        data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatMed)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatMed)

# When modelling the quadratic I haven't considered including AboveMedn with it
# This posed a bit of a puzzle initially
# But let's say I have some data that show a different slope on either side of the median
# the quadratic term may put the vertex of the parabola away from the centre because of the differences in slope
# (It may also do this because the niche centre does not match with the centre of the range)

# If I were to use just the scaled data without converting those below the median to negative
# I might end up with a W or V shape to predicted values when I include abovemedn as a factor

# Also, if I do not include abovemedn I will not be explaining the variance as well as I could
# particularly if the vertex is not central
# One side would appear more linear while the other would be clearly quadratic, 
# but these would not be aligned and would show as unexplained variance

# this is easier to illustrate on paper, but it's also visible in model outputs
# Compare the following models 

# This is the original quadratic model with transformed data
summary(LM_IUCN_Species$Fixed$Quad)

# This is the model with untransformed, scaled data
LM_IUCN_Species$Fixed$Quad_Drop <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropScaled, 2),
                                             data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Quad_Drop)
AIC(LM_IUCN_Species$Fixed$Quad_Drop, LM_IUCN_Species$Fixed$Quad)

# The AIC jumps up between Quad and Quad_Drop because the two sides of the quadratic are not aligned

# This is the model with untransformed, scaled data with abovemedn to account for asymmetry
LM_IUCN_Species$Fixed$QuadAsym_Drop <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropScaled, 2):AboveMedn,
                                                 data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QuadAsym_Drop)
AIC(LM_IUCN_Species$Fixed$QuadAsym_Drop, LM_IUCN_Species$Fixed$Quad_Drop, LM_IUCN_Species$Fixed$Quad)

# The AIC is reduced again by adding abovemedn as a factor to account for misalignment, 
# but this approach means that most of the work that AboveMedn is doing is accounting for this mismatch
# and the AIC isn't even as low as the original Quad
# If we add AboveMedn to the correctly formatted data with the appropriate quadratic (and therefore vertex)
# it can now do what we really want it to which is test is there are differences in slope between the two sides

LM_IUCN_Species$Fixed$QuadAsym <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2):AboveMedn,
                                            data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QuadAsym)
AIC(LM_IUCN_Species$Fixed$QuadAsym, LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$QuadAsym_Drop)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$QuadAsym)

# When we model the data this way, the AIC is lowered further but only marginally
# Also the relative likelihood of Quad is quite high (0.73)
# Don't include AboveMedn

#_______________________________________________________________________________
LM_IUCN_Species$Fixed$LatMedInt <- glm(formula = Prevalence ~ LatitudeScaled * MedianPropScaled,
                                             data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatMedInt)
AIC(LM_IUCN_Species$Fixed$LatMedInt, LM_IUCN_Species$Fixed$LatMed)
relative_likelihood(LM_IUCN_Species$Fixed$LatMedInt, LM_IUCN_Species$Fixed$LatMed)
# AIC is marginally lowered and relative likelihood is not great (0.46)
# Only accept interaction between Latitude and MedianProp if quadratic interaction is alright

LM_IUCN_Species$Fixed$QuadInt <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2),
                                           data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QuadInt)
AIC(LM_IUCN_Species$Fixed$QuadInt, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$QuadInt, LM_IUCN_Species$Fixed$Quad)
# The difference in AIC is a bit better here, but the relative likelihood is not great

LM_IUCN_Species$Fixed$QuadAsymInt <- glm(formula = Prevalence ~ LatitudeScaled * poly(MedianPropSquared, 2):AboveMedn,
                                               data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QuadAsymInt)
AIC(LM_IUCN_Species$Fixed$QuadInt, LM_IUCN_Species$Fixed$QuadAsymInt)
relative_likelihood(LM_IUCN_Species$Fixed$QuadInt, LM_IUCN_Species$Fixed$QuadAsymInt)
# Marginal difference in AIC again, and relative likelihood is high

### Climatic niche #############################################################
# I expect that Latitude and MedianProp effects are driven by climatic niche
# So there's potential that the pca axes will perform better than geographic variables

#_______________________________________________________________________________
# First lets compare Latitude and PCA axes as conceptual twins
LM_IUCN_Species$Fixed$Axis1 <- glm(formula = Prevalence ~ Axis1Scaled,
                                  data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Axis1)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axis1, LM_IUCN_Species$Fixed$Lat)
# Latitude outperforms PCA axes 

LM_IUCN_Species$Fixed$Axis2 <- glm(formula = Prevalence ~ Axis2Scaled,
                                   data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Axis2)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axis2, LM_IUCN_Species$Fixed$Lat)
# Latitude outperforms PCA axes 

LM_IUCN_Species$Fixed$Axes <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled,
                                        data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Axes)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$Axes, LM_IUCN_Species$Fixed$Lat)
# Latitude outperforms PCA axes 

#_______________________________________________________________________________
# What about alongside each other to see how much overlap there is?
LM_IUCN_Species$Fixed$LatAxis1 <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled,
                                            data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatAxis1)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$LatAxis1, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$LatAxis1, LM_IUCN_Species$Fixed$Lat)
# reduces AIC marginally, but relative likelihood is not good enough

LM_IUCN_Species$Fixed$LatAxis2 <- glm(formula = Prevalence ~ LatitudeScaled + Axis2Scaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatAxis2)
summary(LM_IUCN_Species$Fixed$Lat)
AIC(LM_IUCN_Species$Fixed$LatAxis2, LM_IUCN_Species$Fixed$Lat)
# Axis 2 does not

LM_IUCN_Species$Fixed$LatAxes <- glm(formula = Prevalence ~ LatitudeScaled + Axis1Scaled + Axis2Scaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatAxes)
summary(LM_IUCN_Species$Fixed$Lat)
summary(LM_IUCN_Species$Fixed$LatAxis1)
AIC(LM_IUCN_Species$Fixed$LatAxes, LM_IUCN_Species$Fixed$Lat)
relative_likelihood(LM_IUCN_Species$Fixed$LatAxes, LM_IUCN_Species$Fixed$Lat)
# reduces AIC marginally, but relative likelihood is not good enough

#_______________________________________________________________________________
# Just to confirm that they don't do better alongside medianprop
LM_IUCN_Species$Fixed$MedAxis1 <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled,
                                            data = GMPD_IUCN_Species, family = binomial)

LM_IUCN_Species$Fixed$MedAxis2 <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis2Scaled,
                                            data = GMPD_IUCN_Species, family = binomial)

LM_IUCN_Species$Fixed$MedAxes <- glm(formula = Prevalence ~ poly(MedianPropSquared, 2) + Axis1Scaled + Axis2Scaled,
                                           data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$MedAxis1)
summary(LM_IUCN_Species$Fixed$MedAxis2)
summary(LM_IUCN_Species$Fixed$MedAxes)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$MedAxis1,
    LM_IUCN_Species$Fixed$MedAxis2,
    LM_IUCN_Species$Fixed$MedAxes,
    LM_IUCN_Species$Fixed$Quad)
# Latitude is definitely better than pca axes

#_______________________________________________________________________________
# Is the same true if we include Distance as a conceptual parallel to MedianProp?
LM_IUCN_Species$Fixed$Dist <- glm(formula = Prevalence ~ CorrectedDistanceScaled,
                                  data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Dist)
summary(LM_IUCN_Species$Fixed$Med)
AIC(LM_IUCN_Species$Fixed$Med, LM_IUCN_Species$Fixed$Dist)
relative_likelihood(LM_IUCN_Species$Fixed$Med, LM_IUCN_Species$Fixed$Dist)
# Med does better than Dist

LM_IUCN_Species$Fixed$LatDist <- glm(formula = Prevalence ~ LatitudeScaled + CorrectedDistanceScaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatDist)
summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatDist)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatDist)
# Med still does better

LM_IUCN_Species$Fixed$AxesDist <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + CorrectedDistanceScaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesDist)
summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$AxesDist)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$AxesDist)
# Med still does better

#_______________________________________________________________________________
# Maybe distance should be a quadratic, like medianprop
LM_IUCN_Species$Fixed$QDist <- glm(formula = Prevalence ~ poly(CorrectedDistance, 2),
                                  data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QDist)
summary(LM_IUCN_Species$Fixed$QMed)
AIC(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$QDist)
relative_likelihood(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$QDist)
# Med still does better than Dist

LM_IUCN_Species$Fixed$LatQDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(CorrectedDistance, 2),
                                     data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatQDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatQDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatQDist)
# Med still does better

LM_IUCN_Species$Fixed$AxesQDist <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + poly(CorrectedDistance, 2),
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesQDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesQDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesQDist)
# Med still does better

#_______________________________________________________________________________
# Trying Distance logged CorrectedDistanceLogged
LM_IUCN_Species$Fixed$LDist <- glm(formula = Prevalence ~ CorrectedDistanceLogged,
                                   data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LDist)
summary(LM_IUCN_Species$Fixed$QMed)
AIC(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$LDist)
relative_likelihood(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$LDist)
# Med still does better than Dist

LM_IUCN_Species$Fixed$LatLDist <- glm(formula = Prevalence ~ LatitudeScaled + CorrectedDistanceLogged,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatLDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatLDist)
# Med still does better

LM_IUCN_Species$Fixed$AxesLDist <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + CorrectedDistanceLogged,
                                       data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesLDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesLDist)
# Med still does better

#_______________________________________________________________________________
# Trying Distance logged and quadratic
LM_IUCN_Species$Fixed$QLDist <- glm(formula = Prevalence ~ poly(CorrectedDistanceLogged, 2),
                                   data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QLDist)
summary(LM_IUCN_Species$Fixed$QMed)
AIC(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$QLDist)
relative_likelihood(LM_IUCN_Species$Fixed$QMed, LM_IUCN_Species$Fixed$LDist)
# Med still does better than Dist

LM_IUCN_Species$Fixed$LatQLDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(CorrectedDistanceLogged, 2),
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatQLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatQLDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$LatQLDist)
# Med still does better

LM_IUCN_Species$Fixed$AxesQLDist <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + poly(CorrectedDistanceLogged, 2),
                                       data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesQLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesQLDist)
relative_likelihood(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$AxesQLDist)
# Med still does better

#_______________________________________________________________________________
# What about density?
LM_IUCN_Species$Fixed$Dens <- glm(formula = Prevalence ~ CorrectedDensityScaled,
                                  data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$Dens)
summary(LM_IUCN_Species$Fixed$Med)
AIC(LM_IUCN_Species$Fixed$Med, LM_IUCN_Species$Fixed$Dens)
relative_likelihood(LM_IUCN_Species$Fixed$Med, LM_IUCN_Species$Fixed$Dens)
# Med does better than Dens

LM_IUCN_Species$Fixed$LatDens <- glm(formula = Prevalence ~ LatitudeScaled + CorrectedDensityScaled,
                                     data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatDens)
summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatDens)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$LatDens)
# Med still does better

LM_IUCN_Species$Fixed$AxesDens <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + CorrectedDensityScaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesDens)
summary(LM_IUCN_Species$Fixed$LatMed)
AIC(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$AxesDens)
relative_likelihood(LM_IUCN_Species$Fixed$LatMed, LM_IUCN_Species$Fixed$AxesDens)
# Geographic variables do better all round

#_______________________________________________________________________________
# Maybe density should be a quadratic, like medianprop
LM_IUCN_Species$Fixed$QDens <- glm(formula = Prevalence ~ poly(CorrectedDensity, 2),
                                   data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$QDens)
summary(LM_IUCN_Species$Fixed$QMed)
AIC(LM_IUCN_Species$Fixed$QDens, LM_IUCN_Species$Fixed$QMed)
relative_likelihood(LM_IUCN_Species$Fixed$QDens, LM_IUCN_Species$Fixed$QMed)
# Med still does better than Dist

LM_IUCN_Species$Fixed$LatQDens <- glm(formula = Prevalence ~ LatitudeScaled + poly(CorrectedDensity, 2),
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$LatQDens)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$LatQDens, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$LatQDens, LM_IUCN_Species$Fixed$Quad)
# Med still does better

LM_IUCN_Species$Fixed$AxesQDens <- glm(formula = Prevalence ~ Axis1Scaled + Axis2Scaled + poly(CorrectedDensity, 2),
                                       data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$AxesQDens)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$AxesQDens, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$AxesQDens, LM_IUCN_Species$Fixed$Quad)
# Med still does better

#_______________________________________________________________________________
# What about axes in the full model

LM_IUCN_Species$Fixed$GeoAxis1 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis1Scaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoAxis1)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoAxis1, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoAxis1, LM_IUCN_Species$Fixed$Quad)
# Axis1 increases AIC

LM_IUCN_Species$Fixed$GeoAxis2 <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + Axis2Scaled,
                                            data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoAxis2)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoAxis2, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoAxis2, LM_IUCN_Species$Fixed$Quad)
# Same goes for Axis2 (just)

#_______________________________________________________________________________
# Trying Distance in the full model
LM_IUCN_Species$Fixed$GeoDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled,
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoDist, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoDist, LM_IUCN_Species$Fixed$Quad)
# rejected

LM_IUCN_Species$Fixed$GeoQDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistance, 2),
                                     data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoQDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoQDist, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoQDist, LM_IUCN_Species$Fixed$Quad)
# Distance does stuff!

LM_IUCN_Species$Fixed$GeoQDistInt <- glm(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2),
                                         data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoQDist)
summary(LM_IUCN_Species$Fixed$GeoQDistInt)
AIC(LM_IUCN_Species$Fixed$Quad, LM_IUCN_Species$Fixed$GeoQDist, LM_IUCN_Species$Fixed$GeoQDistInt)
relative_likelihood(LM_IUCN_Species$Fixed$GeoQDist, LM_IUCN_Species$Fixed$GeoQDistInt)
# And it interacts with Latitude!

LM_IUCN_Species$Fixed$VeryFull <- glm(formula = Prevalence ~ LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)),
                                         data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$VeryFull)
summary(LM_IUCN_Species$Fixed$GeoQDistInt)
AIC(LM_IUCN_Species$Fixed$Quad, 
    LM_IUCN_Species$Fixed$GeoQDist, 
    LM_IUCN_Species$Fixed$GeoQDistInt,
    LM_IUCN_Species$Fixed$VeryFull)
relative_likelihood(LM_IUCN_Species$Fixed$VeryFull, LM_IUCN_Species$Fixed$GeoQDistInt)
# And the interaction between Lat and MedianProp can now be accepted!

LM_IUCN_Species$Fixed$VeryFull_Weighted <- glm(formula = Prevalence ~ LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)),
                                      data = GMPD_IUCN_Species, family = binomial,
                                      weights = HostsSampled)

summary(LM_IUCN_Species$Fixed$VeryFull_Weighted)
summary(LM_IUCN_Species$Fixed$GeoQDistInt)
AIC(LM_IUCN_Species$Fixed$Quad, 
    LM_IUCN_Species$Fixed$GeoQDist, 
    LM_IUCN_Species$Fixed$GeoQDistInt,
    LM_IUCN_Species$Fixed$VeryFull)
relative_likelihood(LM_IUCN_Species$Fixed$VeryFull, LM_IUCN_Species$Fixed$GeoQDistInt)
# And the interaction between Lat and MedianProp can now be accepted!

#_______________________________________________________________________________
# Trying Logged Distance in the full model
LM_IUCN_Species$Fixed$GeoLDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceLogged,
                                     data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoLDist, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoLDist, LM_IUCN_Species$Fixed$Quad)
# rejected

LM_IUCN_Species$Fixed$GeoQLDist <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistanceLogged, 2),
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoQLDist)
summary(LM_IUCN_Species$Fixed$Quad)
AIC(LM_IUCN_Species$Fixed$GeoQLDist, LM_IUCN_Species$Fixed$Quad)
relative_likelihood(LM_IUCN_Species$Fixed$GeoQLDist, LM_IUCN_Species$Fixed$Quad)
# Distance doesn't do stuff when logged :'(

#_______________________________________________________________________________
LM_IUCN_Species$Fixed$GeoDens <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDensityScaled,
                                     data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoDens)
summary(LM_IUCN_Species$Fixed$Quad)
# AIC isn't lowered

LM_IUCN_Species$Fixed$GeoQDens <- glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDensity, 2),
                                      data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Fixed$GeoQDens)
summary(LM_IUCN_Species$Fixed$Quad)
# Density doesn't do stuff



## Random effects exploration ##################################################
# Working up the chain of fixed effect structures (in order of AIC) 
# with increasingly complex random effect structures
# Doing it this way to be more thorough and check that fixed effects hold with random effects added
# Lat < \LatMed\ < Quad < \GeoDist\ < GeoQDist < GeoQDistInt < VeryFull

# Little reminder of random effect specification:
# 1 is slope, 0 is intercept, || means model does not estimate covariance, / or : are nested effects
# https://stats.stackexchange.com/questions/608525/choosing-random-effects-to-include-in-a-linear-mixed-model
# Sometimes better to only estimate intercept if the random effects don't fit the data well 
# TODO: Check I'm specifying nested effects right

### Hosts ######################################################################
#_______________________________________________________________________________
LM_IUCN_Species$Host$Lat <- lme4::glmer(formula = Prevalence ~ LatitudeScaled +
                                          (1 + LatitudeScaled|HostCorrectedName),
                                        data = GMPD_IUCN_Species, family = binomial)

LM_IUCN_Species$Host$LatIO <- lme4::glmer(formula = Prevalence ~ LatitudeScaled +
                                          (1|HostCorrectedName),
                                        data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Host$Lat)
summary(LM_IUCN_Species$Host$LatIO)
summary(LM_IUCN_Species$Fixed$Lat)

anova(LM_IUCN_Species$Host$Lat, LM_IUCN_Species$Fixed$Lat)

#_______________________________________________________________________________
LM_IUCN_Species$Host$LatMed <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                                             (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName),
                                           data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Host$LatMed)
# Singular fit already
# removing slopes to cope with hosts

LM_IUCN_Species$Host$LatMedIO <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                                             (1|HostCorrectedName),
                                        data = GMPD_IUCN_Species, family = binomial)

summary(LM_IUCN_Species$Host$LatMedIO)
summary(LM_IUCN_Species$Fixed$LatMed)

anova(LM_IUCN_Species$Host$LatMedIO, LM_IUCN_Species$Fixed$LatMed)
# Works okay

#_______________________________________________________________________________
LM_IUCN_Species$Host$Quad <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                               (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                                             data = GMPD_IUCN_Species, family = binomial),
                                 silent = TRUE)

summary(LM_IUCN_Species$Host$Quad)
# Failed to converge

LM_IUCN_Species$Host$QuadIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                               (1|HostCorrectedName),
                                             data = GMPD_IUCN_Species, family = binomial),
                                 silent = TRUE)

summary(LM_IUCN_Species$Host$QuadIO)
summary(LM_IUCN_Species$Fixed$Quad)

anova(LM_IUCN_Species$Host$QuadIO, LM_IUCN_Species$Fixed$Quad)

#_______________________________________________________________________________
# LM_IUCN_Species$Host$GeoDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
#                                                   poly(MedianPropSquared, 2) +
#                                                   CorrectedDistanceScaled +
#                                                   (1 + LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled|HostCorrectedName),
#                                                 data = GMPD_IUCN_Species, family = binomial),
#                                     silent = TRUE)
# 
# summary(LM_IUCN_Species$Host$GeoDist)
# # Failed to converge
# 
# LM_IUCN_Species$Host$GeoDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
#                                                     poly(MedianPropSquared, 2) +
#                                                     CorrectedDistanceScaled +
#                                                    (1|HostCorrectedName),
#                                                  data = GMPD_IUCN_Species, family = binomial),
#                                      silent = TRUE)
# 
# summary(LM_IUCN_Species$Host$GeoDistIO)
# summary(LM_IUCN_Species$Fixed$GeoDist)
# # okay

#_______________________________________________________________________________
LM_IUCN_Species$Host$GeoQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                   poly(MedianPropSquared, 2) +
                                                   poly(CorrectedDistance, 2) +
                                                   (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistance, 2)|HostCorrectedName),
                                                 data = GMPD_IUCN_Species, family = binomial),
                                     silent = TRUE)

summary(LM_IUCN_Species$Host$GeoQDist)
# Fail

LM_IUCN_Species$Host$GeoQDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                     poly(MedianPropSquared, 2) +
                                                     poly(CorrectedDistance, 2) +
                                                     (1|HostCorrectedName),
                                                   data = GMPD_IUCN_Species, family = binomial),
                                       silent = TRUE)

summary(LM_IUCN_Species$Host$GeoQDistIO)
summary(LM_IUCN_Species$Fixed$GeoQDist)

anova(LM_IUCN_Species$Host$GeoQDistIO, LM_IUCN_Species$Fixed$GeoQDist)
# Okay

#_______________________________________________________________________________
LM_IUCN_Species$Host$GeoQDistInt <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                      poly(CorrectedDistance, 2) + 
                                                      poly(MedianPropSquared, 2) +
                                                      (1 + LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|HostCorrectedName),
                                                    data = GMPD_IUCN_Species, family = binomial),
                                        silent = TRUE)

summary(LM_IUCN_Species$Host$GeoQDistInt)
# fail

LM_IUCN_Species$Host$GeoQDistIntIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2) +
                                                        (1|HostCorrectedName),
                                                      data = GMPD_IUCN_Species, family = binomial),
                                          silent = TRUE)

summary(LM_IUCN_Species$Host$GeoQDistIntIO)
summary(LM_IUCN_Species$Fixed$GeoQDistInt)
# singular fit 

#_______________________________________________________________________________
#  VeryFull
LM_IUCN_Species$Host$VeryFull <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                   (poly(CorrectedDistance, 2) + 
                                                      poly(MedianPropSquared, 2)) +
                                                   (1 + LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2))|HostCorrectedName),
                                                 data = GMPD_IUCN_Species, family = binomial),
                                     silent = TRUE)

summary(LM_IUCN_Species$Host$VeryFull)
# fail

LM_IUCN_Species$Host$VeryFullIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                     (poly(CorrectedDistance, 2) + 
                                                        poly(MedianPropSquared, 2)) +
                                                     (1|HostCorrectedName),
                                                   data = GMPD_IUCN_Species, family = binomial),
                                          silent = TRUE)

summary(LM_IUCN_Species$Host$VeryFullIO)
summary(LM_IUCN_Species$Fixed$VeryFull)
# singular fit

### Parasites ###################################################################
# Different approach to hosts because I'm trying to choose which grouping factor is most informative
# Still building up along the same route, but comparing only between ParType and ParPhylum

# https://stats.stackexchange.com/questions/117497/comparing-between-random-effects-structures-in-a-linear-mixed-effects-model
# Can't compare AIC as I usually would for fixed effects so could use anova instead
# https://stats.stackexchange.com/questions/418593/understanding-anova-to-compare-mixed-model-with-a-gzlm

# Can't use anova because there is no difference in the number of parameters between models
# Could then return to comparing with AIC? The reason it's not used for random effects is 
# because it penalises according to the number of parameters and that's not the point for random
# But then Phylum will always come out on top just because it has more levels

# I'm going to go back to taking a more thoughtful approach

GMPD_IUCN_Species %>% 
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

# First need to do some data prep by removing rejected groups
# Although Sarcomastigophora could be included in Type, I wouldn't be able to compare AIC with Phylum
GMPD_IUCN_Species_ParType <- GMPD_IUCN_Species %>% 
  dplyr::filter(!ParType %in% c("Prion", "Fungus") & !is.na(ParPhylum) )


#_______________________________________________________________________________
LM_IUCN_Species$ParType$Lat <- lme4::glmer(formula = Prevalence ~ LatitudeScaled +
                                             (1 + LatitudeScaled|ParType),
                                           data = GMPD_IUCN_Species_ParType, family = binomial)

summary(LM_IUCN_Species$ParType$Lat)
# Boundary fit is singular

LM_IUCN_Species$ParType$LatIO <- lme4::glmer(formula = Prevalence ~ LatitudeScaled +
                                               (1|ParType),
                                             data = GMPD_IUCN_Species_ParType, family = binomial)

summary(LM_IUCN_Species$ParType$LatIO)

#_______________________________________________________________________________
# LM_IUCN_Species$ParType$LatMed <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
#                                                     (1 + LatitudeScaled + MedianPropScaled|ParType),
#                                                   data = GMPD_IUCN_Species_ParType, family = binomial),
#                                       silent = TRUE)
# 
# summary(LM_IUCN_Species$ParType$LatMed)
# # singular fit
# 
# LM_IUCN_Species$ParType$LatMedIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
#                                                     (1|ParType),
#                                                   data = GMPD_IUCN_Species_ParType, family = binomial),
#                                       silent = TRUE)
# 
# 
# summary(LM_IUCN_Species$ParType$LatMedIO)

#_______________________________________________________________________________
LM_IUCN_Species$ParType$Quad <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                                  (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType),
                                                data = GMPD_IUCN_Species_ParType, family = binomial),
                                    silent = TRUE)

summary(LM_IUCN_Species$ParType$Quad)
# failure to converge

LM_IUCN_Species$ParType$QuadIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                                    (1|ParType),
                                                  data = GMPD_IUCN_Species_ParType, family = binomial),
                                      silent = TRUE)

summary(LM_IUCN_Species$ParType$QuadIO)
# singular fit for parasite type

#_______________________________________________________________________________
LM_IUCN_Species$ParType$GeoDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                     poly(MedianPropSquared, 2) +
                                                     CorrectedDistanceScaled +
                                                     (1 + LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled|ParType),
                                                   data = GMPD_IUCN_Species_ParType, family = binomial),
                                       silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoDist)
# fail

LM_IUCN_Species$ParType$GeoDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                       poly(MedianPropSquared, 2) +
                                                       CorrectedDistanceScaled +
                                                       (1|ParType),
                                                     data = GMPD_IUCN_Species_ParType, family = binomial),
                                         silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoDistIO)
# singular fit

#_______________________________________________________________________________
LM_IUCN_Species$ParType$GeoQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                      poly(MedianPropSquared, 2) +
                                                      poly(CorrectedDistance, 2) +
                                                      (1 + LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistance, 2)|ParType),
                                                    data = GMPD_IUCN_Species_ParType, family = binomial),
                                        silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoQDist)
# fail

LM_IUCN_Species$ParType$GeoQDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                        poly(MedianPropSquared, 2) +
                                                        poly(CorrectedDistance, 2) +
                                                        (1|ParType),
                                                      data = GMPD_IUCN_Species_ParType, family = binomial),
                                          silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoQDistIO)
# singular fit

#_______________________________________________________________________________
LM_IUCN_Species$ParType$GeoQDistInt <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                         poly(CorrectedDistance, 2) + 
                                                         poly(MedianPropSquared, 2) +
                                                         (1 + LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|ParType),
                                                       data = GMPD_IUCN_Species_ParType, family = binomial),
                                           silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoQDistInt)
# fail

LM_IUCN_Species$ParType$GeoQDistIntIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                           poly(CorrectedDistance, 2) + 
                                                           poly(MedianPropSquared, 2) +
                                                           (1|ParType),
                                                         data = GMPD_IUCN_Species_ParType, family = binomial),
                                             silent = TRUE)

summary(LM_IUCN_Species$ParType$GeoQDistIntIO)
# singular fit 

#_______________________________________________________________________________
LM_IUCN_Species$ParType$VeryFull <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                      (poly(CorrectedDistance, 2) + 
                                                         poly(MedianPropSquared, 2)) +
                                                      (1 + LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2))|ParType),
                                                    data = GMPD_IUCN_Species_ParType, family = binomial),
                                        silent = TRUE)

summary(LM_IUCN_Species$ParType$VeryFull)
# fail

LM_IUCN_Species$ParType$VeryFullIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                        (poly(CorrectedDistance, 2) + 
                                                           poly(MedianPropSquared, 2)) +
                                                        (1|ParType),
                                                      data = GMPD_IUCN_Species_ParType, family = binomial),
                                          silent = TRUE)

summary(LM_IUCN_Species$ParType$VeryFullIO)
# fail 

### Both hosts and parasites ####
# Not likely to get far with this but worth a try anyway

#_______________________________________________________________________________
LM_IUCN_Species$Both$Lat <- lme4::glmer(formula = Prevalence ~ LatitudeScaled +
                                          (1 + LatitudeScaled|ParType) +
                                          (1 + LatitudeScaled|HostCorrectedName),
                                        data = GMPD_IUCN_Species_ParType, family = binomial)
# Not a singular fit!

LM_IUCN_Species$Host$LatType <- lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                              (1 + LatitudeScaled|HostCorrectedName),
                                            data = GMPD_IUCN_Species_ParType, family = binomial)

LM_IUCN_Species$Fixed$LatType <- glm(formula = Prevalence ~ LatitudeScaled,
                                     data = GMPD_IUCN_Species_ParType, family = binomial)

summary(LM_IUCN_Species$Both$Lat)
anova(LM_IUCN_Species$Both$Lat, 
      LM_IUCN_Species$Host$LatType)
relative_likelihood(LM_IUCN_Species$Both$Lat, 
                    LM_IUCN_Species$Host$LatType)

#_______________________________________________________________________________
LM_IUCN_Species$Both$LatMed <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                                                 (1 + LatitudeScaled + MedianPropScaled|ParType) +
                                                 (1 + LatitudeScaled + MedianPropScaled|HostCorrectedName),
                                               data = GMPD_IUCN_Species_ParType, family = binomial),
                                   silent = TRUE)

summary(LM_IUCN_Species$Both$LatMed)
# failed to converge

LM_IUCN_Species$Both$LatMedIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                                                   (1|ParType) +
                                                   (1|HostCorrectedName),
                                                 data = GMPD_IUCN_Species_ParType, family = binomial),
                                     silent = TRUE)

LM_IUCN_Species$Host$LatMedType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled +
                                                     (1|HostCorrectedName),
                                                   data = GMPD_IUCN_Species_ParType, family = binomial),
                                       silent = TRUE)

LM_IUCN_Species$Fixed$LatMedType <- glm(formula = Prevalence ~ LatitudeScaled + MedianPropScaled,
                                        data = GMPD_IUCN_Species_ParType, family = binomial)

summary(LM_IUCN_Species$Both$LatMedIO)
anova(LM_IUCN_Species$Both$LatMedIO, 
      LM_IUCN_Species$Host$LatMedType)
relative_likelihood(LM_IUCN_Species$Both$LatMedIO, 
                    LM_IUCN_Species$Host$LatMedType)

#_______________________________________________________________________________
LM_IUCN_Species$Both$Quad <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                               (1 + LatitudeScaled + poly(MedianPropSquared, 2)|ParType) +
                                               (1 + LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                                             data = GMPD_IUCN_Species_ParType, family = binomial),
                                 silent = TRUE)

summary(LM_IUCN_Species$Both$Quad)
# failure to converge

# TODO: tidy up order of terms
LM_IUCN_Species$Both$QuadIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                                 (1|ParType) +
                                                 (1|HostCorrectedName),
                                               data = GMPD_IUCN_Species_ParType, family = binomial),
                                   silent = TRUE)

LM_IUCN_Species$Host$QuadType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) +
                                                   (1|HostCorrectedName),
                                                 data = GMPD_IUCN_Species_ParType, family = binomial),
                                     silent = TRUE)

LM_IUCN_Species$Fixed$QuadType <- try(glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2),
                                          data = GMPD_IUCN_Species_ParType, family = binomial),
                                      silent = TRUE)

summary(LM_IUCN_Species$Both$QuadIO)
anova(LM_IUCN_Species$Both$QuadIO, 
      # LM_IUCN_Species$Fixed$QuadType, 
      LM_IUCN_Species$Host$QuadType)
relative_likelihood(LM_IUCN_Species$Both$QuadIO, 
                    LM_IUCN_Species$Host$QuadType)

#_______________________________________________________________________________
LM_IUCN_Species$Both$GeoDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                  poly(MedianPropSquared, 2) + CorrectedDistanceScaled +
                                                  (1 + LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled|ParType) +
                                                  (1 + LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled|HostCorrectedName),
                                                data = GMPD_IUCN_Species_ParType, family = binomial),
                                    silent = TRUE)

summary(LM_IUCN_Species$Both$GeoDist)
# fail

LM_IUCN_Species$Both$GeoDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled +
                                                    (1|ParType) +
                                                    (1|HostCorrectedName),
                                                  data = GMPD_IUCN_Species_ParType, family = binomial),
                                      silent = TRUE)

LM_IUCN_Species$Host$GeoDistType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + CorrectedDistanceScaled +
                                                      (1|HostCorrectedName),
                                                    data = GMPD_IUCN_Species_ParType, family = binomial),
                                        silent = TRUE)

LM_IUCN_Species$Fixed$GeoDistType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                       poly(MedianPropSquared, 2) + 
                                                       CorrectedDistanceScaled,
                                                     data = GMPD_IUCN_Species_ParType, family = binomial),
                                         silent = TRUE)

summary(LM_IUCN_Species$Both$GeoDistIO)

#_______________________________________________________________________________
LM_IUCN_Species$Both$GeoQDist <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                   poly(MedianPropSquared, 2) +
                                                   poly(CorrectedDistance, 2) +
                                                   (1 + LatitudeScaled + poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|ParType) +
                                                   (1 + LatitudeScaled + poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|HostCorrectedName),
                                                 data = GMPD_IUCN_Species_ParType, family = binomial),
                                     silent = TRUE)

summary(LM_IUCN_Species$Both$GeoQDist)
# fail

LM_IUCN_Species$Both$GeoQDistIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                     poly(MedianPropSquared, 2) + 
                                                     poly(CorrectedDistance, 2) +
                                                     (1|ParType) +
                                                     (1|HostCorrectedName),
                                                   data = GMPD_IUCN_Species_ParType, family = binomial),
                                       silent = TRUE)

LM_IUCN_Species$Host$GeoQDistType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled + 
                                                       poly(MedianPropSquared, 2) + 
                                                       poly(CorrectedDistance, 2) +
                                                       (1|HostCorrectedName),
                                                     data = GMPD_IUCN_Species_ParType, family = binomial),
                                         silent = TRUE)

LM_IUCN_Species$Fixed$GeoQDistType <- try(glm(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + poly(CorrectedDistance, 2),
                                              data = GMPD_IUCN_Species_ParType, family = binomial),
                                          silent = TRUE)

summary(LM_IUCN_Species$Both$GeoQDistIO)
summary(LM_IUCN_Species$Host$GeoQDistType)
anova(LM_IUCN_Species$Both$QuadIntIO, 
      LM_IUCN_Species$Host$QuadIntType,
      LM_IUCN_Species$Fixed$QuadIntType)

#_______________________________________________________________________________
LM_IUCN_Species$Both$GeoQDistInt <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                      poly(CorrectedDistance, 2) + 
                                                      poly(MedianPropSquared, 2) +
                                                      (1 + LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|ParType) +
                                                      (1 + LatitudeScaled * poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2)|HostCorrectedName),
                                                    data = GMPD_IUCN_Species_ParType, family = binomial),
                                        silent = TRUE)

summary(LM_IUCN_Species$Both$GeoQDistInt)
# fail

LM_IUCN_Species$Both$GeoQDistIntIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                        poly(CorrectedDistance, 2) + 
                                                        poly(MedianPropSquared, 2) +
                                                        (1|ParType) +
                                                        (1|HostCorrectedName),
                                                      data = GMPD_IUCN_Species_ParType, family = binomial),
                                          silent = TRUE)
# Also fail

LM_IUCN_Species$Host$GeoQDistIntType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                          poly(CorrectedDistance, 2) + 
                                                          poly(MedianPropSquared, 2) +
                                                          (1|HostCorrectedName),
                                                        data = GMPD_IUCN_Species_ParType, family = binomial),
                                            silent = TRUE)
# Also fail

LM_IUCN_Species$Fixed$GeoQDistIntType <- try(glm(formula = Prevalence ~ LatitudeScaled * 
                                                   poly(CorrectedDistance, 2) + 
                                                   poly(MedianPropSquared, 2),
                                                 data = GMPD_IUCN_Species_ParType, family = binomial),
                                             silent = TRUE)
# Also fail


summary(LM_IUCN_Species$Both$GeoQDistIntIO)

#_______________________________________________________________________________
LM_IUCN_Species$Both$VeryFull <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                      (poly(CorrectedDistance, 2) + 
                                                         poly(MedianPropSquared, 2)) +
                                                      (1 + LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2))|ParType) +
                                                      (1 + LatitudeScaled * (poly(CorrectedDistance, 2) + poly(MedianPropSquared, 2))|HostCorrectedName),
                                                    data = GMPD_IUCN_Species_ParType, family = binomial),
                                        silent = TRUE)

summary(LM_IUCN_Species$Both$VeryFull)
# fail

LM_IUCN_Species$Both$VeryFullIO <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                        (poly(CorrectedDistance, 2) + 
                                                           poly(MedianPropSquared, 2)) +
                                                        (1|ParType) +
                                                        (1|HostCorrectedName),
                                                      data = GMPD_IUCN_Species_ParType, family = binomial),
                                          silent = TRUE)
# Also fail

LM_IUCN_Species$Host$VeryFullType <- try(lme4::glmer(formula = Prevalence ~ LatitudeScaled * 
                                                          (poly(CorrectedDistance, 2) + 
                                                             poly(MedianPropSquared, 2)) +
                                                          (1|HostCorrectedName),
                                                        data = GMPD_IUCN_Species_ParType, family = binomial),
                                            silent = TRUE)
# Also fail

LM_IUCN_Species$Fixed$VeryFullType <- try(glm(formula = Prevalence ~ LatitudeScaled * 
                                                   (poly(CorrectedDistance, 2) + 
                                                      poly(MedianPropSquared, 2)),
                                                 data = GMPD_IUCN_Species_ParType, family = binomial),
                                             silent = TRUE)
# Also fail


summary(LM_IUCN_Species$Both$VeryFullIO)

# TODO: check all models are run
# TODO: sort out all summaries etc.



