# Analysing distance metrics (without climate)
# Written by Margaret Bolton mb804(at)exeter.ac.uk

# Libraries and data ###########################################################
library(MCMCglmm)
library(lme4)
library(here)
source(here::here("Functions.R"))

GMPD_Distances_Data <- read.csv(here::here("Data/Data back ups/GMPD_Distances_Data_03.csv"), header = TRUE, stringsAsFactors = FALSE)

# summary:
# linear models of parasitism rates using range position metrics and latitude
# first done using lme4 but later converted to mcmcglmm to cope with gappy data
# naming system of models:
## Md_[iucn/gbif]_[species/subgroup]_[...]

# IUCN species #################################################################

## Prepping data ####
GMPD_IUCN_Species <- GMPD_Distances_Data %>%
  dplyr::filter(RestrAll, 
                RangeMethod == "iucn",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared))

## lme4 models ####
# TODO: summary of models in this section

### 1: Latitude + MedianProp #### 
# Latitude and proportional distance to range median modelled as simply as possible 
# with host species as a random effect
Md_IUCN_Species_1 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName),
                  data = GMPD_IUCN_Species, family = binomial)

summary(Md_IUCN_Species_1)

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
# AIC is 8781.6 so definitely no improvement again, and the quadractic term is non-significant
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
         MedianPropSquScaled    = base::scale(MedianPropSquared))

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
# May still be worth investigating a little further:
# TODO: Look at sample size and range size effects for subspecies

#### visualisation model 1
# TODO: model 1 visualisation

# GBIF species #################################################################

## Prepping data ####
GMPD_GBIF_Species <- GMPD_Distances_Data %>%
  dplyr::filter(CleanAll, 
                RangeMethod == "gbif",
                RangeTaxonLvl == "species") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared))

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
# Both latitude and proportional distance to range median are significant.
# Could model proportional distance to range median better as it's a bit simplistic here
# There are likely differences between range hemispheres and potentially a quadratic pattern

#### visualisation model 1
# TODO: model 1 visualisation

### 2: + AboveMedn ####
# Latitude and proportional distance to range median modelled as in 1,
# This time including above or below median as a random effect for proportional distance to median
# (It's not relevant to include it for latitude)
Md_GBIF_Species_2 <- glmer(formula = Prevalence ~ LatitudeScaled + MedianPropScaled + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + MedianPropScaled|HostCorrectedName) + 
                             (0+MedianPropScaled|AboveMedn),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_2)

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
Md_GBIF_Species_3 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2)|HostCorrectedName),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_3)

# notes:
# AIC is 8781.6 so definitely no improvement again, and the quadractic term is non-significant
# Should stick with the simplest, model 1
# I previously tried the same with only poly 1 in the random effects but AIC was comparable to m2

#### visualisation model 3
# TODO: model 3 visualisation

### 4: raw quadratic ####
# The same as model 3 but with raw rather than the orthogonal polynomial
Md_GBIF_Species_4 <- glmer(formula  = Prevalence ~ LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE) + 
                             (1|HostCorrectedName) + 
                             (0+LatitudeScaled + poly(MedianPropSquared, 2, raw = TRUE)|HostCorrectedName),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_4)

# notes:
# this one was just for fun so I could see the difference in estimates between poly and raw 
# Model fails to converge

### 5: scaled and squared ####
# Essentially the same as model 3 except proportional distance to median has been scaled
Md_GBIF_Species_5 <- glmer(formula = Prevalence ~ LatitudeScaled + poly(MedianPropSquScaled, 2) +
                             (1|HostCorrectedName) +
                             (0+LatitudeScaled + poly(MedianPropSquScaled, 2)|HostCorrectedName),
                           data = GMPD_GBIF_Species, family = binomial)

summary(Md_GBIF_Species_5)

# notes:
# model fails to converge but gives basically the same output as model 3

# GBIF subgroups ###############################################################

## Prepping data ####
GMPD_GBIF_Subgroup <- GMPD_Distances_Data %>%
  dplyr::filter(CleanSub, 
                RangeMethod == "gbif",
                RangeTaxonLvl == "subgroup") %>%
  mutate(LatitudeScaled         = base::scale(abs(Latitude)),
         EquatorwardsPropScaled = base::scale(EquatorwardsProp),
         MedianPropScaled       = base::scale(MedianProp),
         MedianPropSquared      = case_when(!AboveMedn ~ MedianProp * -1,
                                            TRUE       ~ MedianProp),
         MedianPropSquScaled    = base::scale(MedianPropSquared))

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
