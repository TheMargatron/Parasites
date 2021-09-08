# Climate niche functions and analysis
# Written by Olivier Broennimann and Blaise Petitpierre. Departement of Ecology and Evolution (DEE).
# University of Lausanne. Switzerland. April 2012.
# Modified by Margaret Bolton mb804(at)exeter.ac.uk 

# Libraries and data ##########################################################################################
library(ade4)
library(adehabitatMA)
library(here)
library(raster)
library(tidyverse)

GBIF_Data <- read.csv(here::here("Data/Data back ups/GBIF_Data.csv"), header = TRUE, stringsAsFactors = FALSE)
GMPD_Data <- read.csv(here::here("Data/Data back ups/GMPD_Data.csv"), header = TRUE, stringsAsFactors = FALSE)

# IUCN and GBIF names are not always matching
GBIF_Data <- GBIF_Data %>%
  mutate(species = case_when(GBIFName == "Neovison vison"  ~ "Mustela vison",
                            GBIFName == "Tragelaphus oryx" ~ "Taurotragus oryx",
                            GBIFName == "Martes pennanti"  ~ "Pekania pennanti",
                            TRUE                           ~ species))

bio.dat <- getData('worldclim', var = 'bio', res = 10)

BIO_050612 <- subset(bio.dat, c(5, 6, 12))

# PCA and climatic niche ##########################################################################################

gmpd.cols <- c("ParasiteCorrectedName", 
               "HostsSampled", 
               "Longitude", 
               "Latitude", 
               "Prevalence", 
               "EquatorwardsProp",
               "EquatorwardsDist",
               "MedianDist",
               "MedianProp")

# extract all bio data
clim.xy <- xyFromCell(BIO_050612, 1:ncell(BIO_050612))
clim <- as.data.frame(na.omit(cbind(clim.xy, raster::extract(BIO_050612, clim.xy)))) 
xVar <- c(3:5)
rm(clim.xy)

pca.cal <- dudi.pca(clim[xVar], center = T, scale = T, scannf = F, nf = 2)

# species.name <- "Genetta genetta" #example for test runs
# spat.dat <- GBIF_Data
# clim.dat <- BIO_050612
# par.dat <- GMPD_Data

ks.test <- function(species.name, spat.dat, clim.dat, par.dat){
  
  spat.dat <- spat.dat %>%
    filter(species == species.name) %>%
    dplyr::select(c("decimalLongitude", "decimalLatitude"))
  
  # (partially) account for sampling bias by gridding species occurence data and extracting filled cells
  rastr10 <- raster(resolution = 10/60)
  rast.sp <- rasterize(spat.dat, rastr10, fun = 'count')
  rast.filled <- Which(rast.sp, cells = TRUE)
  
  # extract bioclim variables for points where species occurs
  occ.xy <- xyFromCell(rast.sp, rast.filled)
  occ.xy <- as.data.frame(na.omit(cbind(occ.xy, raster::extract(x = clim.dat, y = occ.xy))))
  
  nvar <- length(xVar)     
  R = 100     
  
  #rbind occ data to glob data and add binary weighting column/vector
  #row.w.env <- c(rep(1, 1:nrow(clim)), rep(0, (nrow(clim) + 1):(nrow(clim) + nrow(occ.sp))))   #1 for background clim, zero for occ clim
  #Is this necessary? 
  
  scores.clim <- suprow(pca.cal, clim[, xVar])$lisup     #The pca scores for all climate
  scores.occ <- suprow(pca.cal, occ.xy[, xVar])$lisup     #The pca scores for current species
  
  gmpd.sp <- par.dat %>%
    filter(HostCorrectedName == species.name) %>%
    dplyr::select(gmpd.cols)
  
  gmpd.xy <- as.data.frame(na.omit(cbind(gmpd.sp, raster::extract(clim.dat, gmpd.sp[, c("Longitude", "Latitude")]))))
  scores.gmpd <- cbind(suprow(pca.cal, gmpd.xy[, c("bio5", "bio6", "bio12")])$lisup, gmpd.xy)      # pca scores for current species from gmpd
  
  z2 <- grid.clim(scores.clim, scores.occ, scores.gmpd, R)    #Niche_dyn_funcs_myversion
  #z2$sp.scores <- scores.gmpd
  return(z2)
}

Hostlist <- unique(GMPD_Data$HostCorrectedName)

ks.out <- lapply(Hostlist, ks.test, spat.dat = GBIF_Data, clim.dat = BIO_050612, par.dat = GMPD_Data)
names(ks.out) <- Hostlist

saveRDS(ks.out, file = "KS_data")
saveRDS(ks.out2, file = "KS_data2")

#ks.out.save <- ks.out
ks.out <- ks.out2

######################################################analysis of ks output
library(lme4)


ks.host.list <- lapply(hostlist,function(host,dat){
  ks <- data.frame("HostCorrectedName" = host,
                   "ClimSpan" = nrow(dat[[host]]$sp),
                   "KSSpan" = length(which(dat[[host]]$w == 1)),
                   "KSPeak.uncor.x" = dat[[host]]$zp.uncor[1],
                   "KSPeak.uncor.y" = dat[[host]]$zp.uncor[2],
                   "KSPeak.cor.x" = dat[[host]]$zp.cor[1],
                   "KSPeak.cor.y" = dat[[host]]$zp.cor[2],
                   "dist.cor.sd" = sd(dat[[host]]$gmpd.cor$distance),
                   "dist.uncor.sd" = sd(dat[[host]]$gmpd.uncor$distance))
}, ks.out)
names(ks.host.list) <- hostlist
ks.host <- do.call(rbind,ks.host.list)


ks.cor.list <- lapply(hostlist, function(host,dat,h.dat){
  HostCorrectedName <- rep(host,nrow(dat[[host]]$gmpd.cor))
  distance.prop <- dat[[host]]$gmpd.cor$distance/h.dat[[host]]$KSSpan
  cor <- cbind(HostCorrectedName,dat[[host]]$gmpd.cor,distance.prop)
},dat=ks.out,h.dat=ks.host.list)
names(ks.cor.list) <- hostlist
ks.cor <- do.call(rbind,ks.cor.list)

ks.uncor.list <- lapply(hostlist, function(host,dat,h.dat){
  HostCorrectedName <- rep(host,nrow(dat[[host]]$gmpd.uncor))
  distance.prop <- dat[[host]]$gmpd.uncor$distance/h.dat[[host]]$KSSpan
  uncor <- cbind(HostCorrectedName,dat[[host]]$gmpd.uncor,distance.prop)
},dat=ks.out,h.dat=ks.host.list)
names(ks.uncor.list) <- hostlist
ks.uncor <- do.call(rbind,ks.uncor.list)

ks.cor.2 <- merge(ks.cor,unique(Pardata60[,c("ParasiteCorrectedName","MicMac")]),by="ParasiteCorrectedName")
ks.cor.2 <- ks.cor.2[which(ks.cor.2$MicMac=="Macro"),]

summary(lm(distance ~ KSSpan,data = merge(ks.cor,ks.host,by="HostCorrectedName")))
summary(glm(Prevalence~HostsSampled,data=ks.cor,family=binomial))
summary(lm(dist.cor.sd~KSSpan,data=ks.host))
summary(lm(HostsSampled~Latitude,data=ks.cor))

temp2 <- lapply(unique(temp$Citation),function(cit,dat){
  if(length(which(dat[,"Citation"]==cit))>1){
    h <- summary(dat[which(dat[,"Citation"]==cit),c("HostAge","HostSex","HostsSampled")])
  }
},dat=temp)

temp <- Pardata[,c("Citation","HostAge","HostSex","HostsSampled")]
temp$Citation <- as.factor(temp$Citation)
temp$HostAge <- as.factor(temp$HostAge)
temp$HostSex <- as.factor(temp$HostSex)

ks.cor$Latitude.sc <- scale(ks.cor$Latitude)
ks.cor$distance.sc <- scale(ks.cor$distance)
ks.cor$prop2equa.sc <- scale(ks.cor$prop2equa)
ks.cor$zp.cor.sc <- scale(ks.cor$zp.cor)
ks.cor$distance.prop.sc <- scale(ks.cor$distance.prop)

ks.model1 <- glm(formula = Prevalence ~ (zp.cor.sc),data=ks.cor,family = binomial)
ks.model1w <- glm(formula = Prevalence ~ (zp.cor.sc),data=ks.cor,family = binomial,weights=HostsSampled)

ks.model2 <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc),data=ks.cor,family = binomial)
ks.model2w <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc),data=ks.cor,family = binomial,weights=HostsSampled)

ks.model3 <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc),data=ks.cor,family = binomial)
ks.model3w <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc),data=ks.cor,family = binomial,weights=HostsSampled)


#ks.model1 <- glm(formula = Prevalence ~ (zp.cor.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial)
ks.model1w <- glm(formula = Prevalence ~ (zp.cor.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)
#ks.model1h <- glmer(formula = Prevalence ~ (zp.cor.sc)+(1+zp.cor.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial)
ks.model1hw <- glmer(formula = Prevalence ~ (zp.cor.sc)+(1+zp.cor.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)

#ks.model2 <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc),data=ks.cor.10,family = binomial)
ks.model2w <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)
#ks.model2h <- glmer(formula = Prevalence ~ (Latitude.sc+zp.cor.sc)+(1+Latitude.sc+zp.cor.sc|HostCorrectedName),data=ks.cor.10,family = binomial)
ks.model2hw <- glmer(formula = Prevalence ~ (Latitude.sc+zp.cor.sc)+(1+Latitude.sc+zp.cor.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)

#ks.model3 <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc),data=ks.cor.10,family = binomial)
ks.model3w <- glm(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)
#ks.model3h <- glmer(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc)+(1+Latitude.sc+zp.cor.sc+prop2equa.sc|HostCorrectedName),data=ks.cor.10,family = binomial)
ks.model3hw <- glmer(formula = Prevalence ~ (Latitude.sc+zp.cor.sc+prop2equa.sc)+(1+Latitude.sc+zp.cor.sc+prop2equa.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)

ks.model4 <- glm(formula = Prevalence ~ (Latitude.sc+prop2equa.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial)
ks.model4w <- glm(formula = Prevalence ~ (Latitude.sc+prop2equa.sc),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)
ks.model4h <- glmer(formula = Prevalence ~ (Latitude.sc+prop2equa.sc)+(1+Latitude.sc+prop2equa.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial)
ks.model4hw <- glmer(formula = Prevalence ~ (Latitude.sc+prop2equa.sc)+(1+Latitude.sc+prop2equa.sc|HostCorrectedName),data=ks.cor.10[which(ks.cor.10$HostsSampled<10000),],family = binomial,weights=HostsSampled)

summary(ks.model1w)
summary(ks.model1hw)
summary(ks.model2w)
summary(ks.model2hw)
summary(ks.model3w)
summary(ks.model3hw)
summary(ks.model4w)
summary(ks.model4hw)
summary(ks.model4)
summary(ks.model4h)



# ks.model1 <- glm(formula = Prevalence ~ (distance.sc),data=ks.cor,family = binomial)
# ks.model1p <- glm(formula = Prevalence ~ (distance.prop.sc),data=ks.cor,family = binomial)
# ks.model1pw <- glm(formula = Prevalence ~ (distance.prop.sc),data=ks.cor,family = binomial,weights = HostsSampled)
# 
# ks.model2 <- glm(formula = Prevalence ~ (Latitude.sc + distance.sc),data=ks.cor,family = binomial)
# ks.model2p <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc),data=ks.cor,family = binomial)
# ks.model2pw <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model3 <- glm(formula = Prevalence ~ (Latitude.sc + distance.sc + prop2equa.sc),data=ks.cor,family = binomial)
# ks.model3p <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc),data=ks.cor,family = binomial)
# ks.model3pw <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model4 <- glm(formula = Prevalence ~ (Latitude.sc + prop2equa.sc),data=ks.cor,family = binomial)
# ks.model4w <- glm(formula = Prevalence ~ (Latitude.sc + prop2equa.sc),data=ks.cor,family = binomial,weights=HostsSampled)
# ks.model5 <- glm(formula = Prevalence ~ (Latitude.sc),data=ks.cor,family = binomial)
# ks.model5w <- glm(formula = Prevalence ~ (Latitude.sc),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model1h <- glmer(formula = Prevalence ~ (distance.sc)+(1|HostCorrectedName)+(0+distance.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model1hp <- glmer(formula = Prevalence ~ (distance.prop.sc)+(1|HostCorrectedName)+(0+distance.prop.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model1hpw <- glmer(formula = Prevalence ~ (distance.prop.sc)+(1|HostCorrectedName)+(0+distance.prop.sc|HostCorrectedName),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model2h <- glmer(formula = Prevalence ~ (Latitude.sc + distance.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model2hp <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model2hpw <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc|HostCorrectedName),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model3h <- glmer(formula = Prevalence ~ (Latitude.sc + distance.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.sc+prop2equa.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model3hp <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc+prop2equa.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model3hpw <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc+prop2equa.sc|HostCorrectedName),data=ks.cor,family = binomial,weights=HostsSampled)
# 
# ks.model4h <- glmer(formula = Prevalence ~ (Latitude.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+prop2equa.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model4hw <- glmer(formula = Prevalence ~ (Latitude.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+prop2equa.sc|HostCorrectedName),data=ks.cor,family = binomial,weights=HostSampled)
# 
# ks.model5h <- glmer(formula = Prevalence ~ (Latitude.sc)+(1|HostCorrectedName)+(0+Latitude.sc|HostCorrectedName),data=ks.cor,family = binomial)
# ks.model5hw <- glmer(formula = Prevalence ~ (Latitude.sc)+(1|HostCorrectedName)+(0+Latitude.sc|HostCorrectedName),data=ks.cor,family = binomial,weights=HostsSampled)

############test plot

plotdataPB<-ks.cor
plotdataPB$fit<-fitted(ks.model1hp0)

b0PB<-fixef(ks.model1hp0)[[1]] #overall intercept
b1PB<-fixef(ks.model1hp0)[[2]] #overall latitude effect
#b2PB<-fixef(ks.model1hp0)[[3]] #overall prop2pole effect

range <- seq(min(abs(ks.cor$distance.prop)),max(abs(ks.cor$distance.prop)),length.out=200)
#p2prange <- seq(min(ks.cor$prop2pole),max(ks.cor$prop2pole),length.out=200)

range.sc <- (seq(min(abs(ks.cor$distance.prop.sc)),max(abs(ks.cor$distance.prop.sc)),length.out=200)-mean(abs(ks.cor$distance.prop.sc)))/(sd(abs(ks.cor$distance.prop.sc)-mean(abs(ks.cor$distance.prop.sc))))
#p2prange.sc <- (seq(min(ks.cor$prop2pole),max(ks.cor$prop2pole),length.out=200)-mean(ks.cor$prop2pole))/(sd(ks.cor$prop2pole-mean(ks.cor$prop2pole)))

biglinedatPB<-data.frame("range"=range,"range.sc"=range.sc)


logitPBlat <- b0PB + (biglinedatPB$range.sc*b1PB)
#logitPBpol <- b0PB + (mean(ks.cor$Latitude.sc)*b1PB) + (biglinedatPB$p2prange.sc*b2PB)
biglinedatPB$linelat <- 1/(1+exp(-logitPBlat))
#biglinedatPB$linepol <- 1/(1+exp(-logitPBpol))

#Now try this
plotPBranef <- ranef(ks.model1hp0)[["HostCorrectedName"]]
plotdataPB$smallatlines <- NA
#plotdataPB$smalp2plines <- NA

mlatPB <- mean(abs(plotdataPB$Latitude))
#mp2pPB <- mean(abs(plotdataPB$prop2pole))

for(i in 1:length(plotPBranef[[1]])){
  p.rows <- which(plotdataPB$ParasiteCorrectedName==rownames(plotPBranef)[[i]])
  #c.meanp2p <- mean(plotdataPB[p.rows,"prop2pole.sc"])
  c.meanlat <- mean(abs(plotdataPB[p.rows,"Latitude.sc"]))
  c.logitPBlat <- b0PB+plotPBranef[[1]][[i]]+  (ks.cor[p.rows,"Latitude.sc"]*(b1PB+plotPBranef[[2]][[i]]))
  #c.logitPBp2p <- b0PB+plotPBranef[[1]][[i]]+ (c.meanlat*(b1PB+plotPBranef[[2]][[i]])) + (ks.cor[p.rows,"prop2pole.sc"]*(b2PB+plotPBranef[[3]][[i]]))
  
  plotdataPB[p.rows,"smallatlines"] <- 1/(1+exp(-c.logitPBlat))
  #plotdataPB[p.rows,"smalp2plines"] <- 1/(1+exp(-c.logitPBp2p))
}


ggplot(data=plotdataPB,aes(abs(Latitude), Prevalence),col=ParasiteCorrectedName)+
  geom_line(aes(y=smallatlines,col=ParasiteCorrectedName),size=0.8)+
  geom_point(aes(col=ParasiteCorrectedName),alpha=0.3)+
  geom_line(data=biglinedatPB,aes(y=linelat,x=latrange),size=1) +
  labs(title="Effect of latitude at mean relative location",subtitle="Random parasite intercept and slope",x="Latitude",y="Parasite prevalence")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none")







ks.uncor$Latitude.sc <- scale(ks.uncor$Latitude)
ks.uncor$distance.sc <- scale(ks.uncor$distance)
ks.uncor$prop2equa.sc <- scale(ks.uncor$prop2equa)
ks.uncor$zp.cor.sc <- scale(ks.uncor$zp.uncor)
ks.uncor$distance.prop.sc <- scale(ks.uncor$distance.prop)


ks.model1u <- glm(formula = Prevalence ~ (distance.sc),data=ks.uncor,family = binomial)
ks.model1pu <- glm(formula = Prevalence ~ (distance.prop.sc),data=ks.uncor,family = binomial)

ks.model2u <- glm(formula = Prevalence ~ (Latitude.sc + distance.sc),data=ks.uncor,family = binomial)
ks.model2pu <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc),data=ks.uncor,family = binomial)

ks.model3u <- glm(formula = Prevalence ~ (Latitude.sc + distance.sc + prop2equa.sc),data=ks.uncor,family = binomial)
ks.model3pu <- glm(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc),data=ks.uncor,family = binomial)


ks.model1hu <- glmer(formula = Prevalence ~ (distance.sc)+(1|HostCorrectedName)+(0+distance.sc|HostCorrectedName),data=ks.uncor,family = binomial)
ks.model1hpu <- glmer(formula = Prevalence ~ (distance.prop.sc)+(1|HostCorrectedName)+(0+distance.prop.sc|HostCorrectedName),data=ks.uncor,family = binomial)

ks.model2hu <- glmer(formula = Prevalence ~ (Latitude.sc + distance.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.sc|HostCorrectedName),data=ks.uncor,family = binomial)
ks.model2hpu <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc|HostCorrectedName),data=ks.uncor,family = binomial)

ks.model3hu <- glmer(formula = Prevalence ~ (Latitude.sc + distance.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.sc+prop2equa.sc|HostCorrectedName),data=ks.uncor,family = binomial)
ks.model3hpu <- glmer(formula = Prevalence ~ (Latitude.sc + distance.prop.sc + prop2equa.sc)+(1|HostCorrectedName)+(0+Latitude.sc+distance.prop.sc+prop2equa.sc|HostCorrectedName),data=ks.uncor,family = binomial)




###########################################
#niche breadth using ade4 and subniche

host.count.stack <- rastr10 <- raster(resolution=10/60)

for (j in hostlist) {
  lyr <- rasterize(as.data.frame(GMPD_dists2[[j]][, c("Longitude","Latitude")]), rastr10, fun = 'count')
  names(lyr) <- j
  host.count.stack <- stack(host.count.stack, lyr)
}

host.count.df <- cbind(clim[, c(1,2)], raster::extract(host.count.stack, clim[, c(1,2)]))
names(host.count.df) <- c("x", "y", hostlist)
host.count.df <- replace(host.count.df, is.na(host.count.df), 0)

global.niche <- niche(pca.cal, host.count.df[, 3:52], scannf = F, nf = 2)
View(global.niche)

#subniche.test <- subniche(global.niche, factor= rep(1,nrow(host.count.df)))
#Can't do this now as I don't have factors to split it into
#But would be nice to do it after imputing parasite richness



