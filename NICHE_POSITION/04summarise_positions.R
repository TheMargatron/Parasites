##### Summary stats and graphs niche and range position of parasite sampling locations #####
library(tidyverse)
dat.wd <- "E:/NON_PROJECT/BI_CLIMATE/MAMMAL_PARASITE/DATA/27thSept2021"

dat <- readRDS(paste0(dat.wd, "/niche_pos_bio1"))

##### Number of parasite sample points per host-parasite combination #####
numpar <- vector()
for(i in 1:length(dat)) {
  if(length(dat[[i]])>1) {
  numpar <- c(numpar, as.vector(table(dat[[i]]$pardist.uncor$ParasiteCorrectedName)))
  }
}
hist(numpar, breaks=seq(min(numpar)-1, max(numpar)+5, 5),
     main="Number of sampling locations\nfor each host-parasite combo") ## Most host-parasite combinations have two or three sampling points. 
          
##### Upper and lower latitudinal boundaries of parasite sampling locations for eac host-parasite combo #####
minlong.par <- maxlong.par <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("par", "lat", "hostpar"))

for(i in names(dat)) {
  if(length(dat[[i]])>1) {
    lat <- dat[[i]]$pardist.uncor[,c("ParasiteCorrectedName", "Latitude")]
    
    lwr <- aggregate(Latitude ~ ParasiteCorrectedName, lat, function(x) min(x))
    lwr[,"host-par"] <- paste0(lwr$ParasiteCorrectedName, "-", i)
    colnames(lwr) <- colnames(minlong.par)
    minlong.par <- rbind(minlong.par, lwr)
    
    upr <- aggregate(Latitude ~ ParasiteCorrectedName, lat, function(x) max(x))
    upr[,"host-par"] <- paste0(upr$ParasiteCorrectedName, "-", i)
    colnames(upr) <- colnames(maxlong.par)
    maxlong.par <- rbind(maxlong.par, upr)
  }
}  
minmaxlong.par <- merge(minlong.par, maxlong.par, by.x=c("par", "hostpar"), by.y=c("par", "hostpar"))
minmaxlong.par <- minmaxlong.par[order(minmaxlong.par$lat.x),]
minmaxlong.par$y <- c(1:nrow(minmaxlong.par))

par(mfrow=c(1,2))
ggplot() +
  geom_segment(data = minmaxlong.par, aes(x=lat.x, y=y, xend=lat.y, yend=y), size = 1) +
  labs(x="Parasite sampling latitude", y="", title="Ordered by lowest latitude")

minmaxlong.par$lat.rng <- minmaxlong.par$lat.y - minmaxlong.par$lat.x
minmaxlong.par <- minmaxlong.par[order(minmaxlong.par$lat.rng),]
minmaxlong.par$y <- c(1:nrow(minmaxlong.par))

ggplot() +
  geom_segment(data = minmaxlong.par, aes(x=lat.x, y=y, xend=lat.y, yend=y), size = 1) +
  labs(x="Parasite sampling latitude", y="", title="Ordered by latitudinal range")

##### Anthro position of the parasite sampling locations#####
anthro.par <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("par", "long", "lat", "ghsl", "hii", "hostpar"))

for(i in names(dat)) {
  if(length(dat[[i]])>1) {
    if(nrow(dat[[i]]$pste.anthro)>0) {
      
      anthro <- as.data.frame(lapply(dat[[i]]$pste.anthro[c("Latitude", "Longitude", "ghsl", "hii")],function(x) {as.numeric(as.character(x))})) ## For some reason numeric columns were saved as factors by code 03, and I can't change it - may be the saveRDS function.
      anthro$ParasiteCorrectedName <- dat[[i]]$pste.anthro$ParasiteCorrectedName
      
      x <- aggregate(ghsl ~ ParasiteCorrectedName + Longitude + Latitude, anthro, function(x) min(x))
      y <- aggregate(hii ~ ParasiteCorrectedName + Longitude + Latitude, anthro, function(x) min(x))
      
      anthro <- merge(x,y)
      anthro[,"hostpar"] <- paste0(anthro$ParasiteCorrectedName, "_", i)
      anthro.par <- rbind(anthro.par, anthro)
      
    }
  }
}  

ggplot(data = anthro.par, aes(x=ParasiteCorrectedName, y=ghsl)) + 
  stat_summary(fun.min=min, fun=median, fun.max=max) # + coord_flip()

ggplot(data = anthro.par, aes(x=ParasiteCorrectedName, y=hii)) + 
  stat_summary(fun.min=min, fun=median, fun.max=max) # + coord_flip()

### Compare anthro conditions at parasite sampling locations to conditions across the host range.
anthro.host <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("host", "long", "lat", "ghsl", "hii"))

for(i in names(dat)) {
  if(length(dat[[i]])>1) {
    if(nrow(dat[[i]]$host.anthro)>0) {
      
      anthro <- as.data.frame(lapply(dat[[i]]$host.anthro[c("decimalLatitude", "decimalLongitude", "ghsl", "hii")],function(x) {as.numeric(as.character(x))})) ## For some reason numeric columns were saved as factors by code 03, and I can't change it - may be the saveRDS function.
      colnames(anthro)[1:2] <- c("Latitude", "Longitude")
      
      x <- aggregate(ghsl ~ Longitude + Latitude, anthro, function(x) min(x))
      y <- aggregate(hii ~ Longitude + Latitude, anthro, function(x) min(x))
      
      anthro <- merge(x,y)
      anthro$host <- i
      anthro.host <- rbind(anthro.host, anthro)
      
    }
  }
}  

### hii
median.host <- aggregate(hii ~ host, anthro.host, function(x) median(x))
min.host <- aggregate(hii ~ host, anthro.host, function(x) min(x))
max.host <- aggregate(hii ~ host, anthro.host, function(x) max(x))

median.par <- aggregate(hii ~ hostpar, anthro.par, function(x) median(x))
min.par <- aggregate(hii ~ hostpar, anthro.par, function(x) min(x))
max.par <- aggregate(hii ~ hostpar, anthro.par, function(x) max(x))

median.par$host <- as.character(unlist(as.data.frame(strsplit(median.par$hostpar, "_"))[2,]))
min.par$host <- as.character(unlist(as.data.frame(strsplit(min.par$hostpar, "_"))[2,]))
max.par$host <- as.character(unlist(as.data.frame(strsplit(max.par$hostpar, "_"))[2,]))

median.par.host <- merge(median.par, median.host, by.x="host", by.y="host")
min.par.host <- merge(min.par, min.host, by.x="host", by.y="host")
max.par.host <- merge(max.par, max.host, by.x="host", by.y="host")

median.par.host$prop <- median.par.host$hii.x - median.par.host$hii.y
hist(median.par.host$prop, main="Median hii of parasite sample\nlocations above or below median\nhii of host")

min.par.host$prop <- (min.par.host$hii.x-min.par.host$hii.y) / (max.par.host$hii.y - min.par.host$hii.y)
hist(min.par.host$prop, main="Proximity of min hii at parasite\nsample locations to min\nhii for host.\n0 is at min")

max.par.host$prop <- max.par.host$hii.x/max.par.host$hii.y
hist(max.par.host$prop, main="Proximity of max hii at parasite\nsample locations to max\nhii for host.\n1 is at max")

### ghsl
median.host <- aggregate(ghsl ~ host, anthro.host, function(x) median(x))
min.host <- aggregate(ghsl ~ host, anthro.host, function(x) min(x))
max.host <- aggregate(ghsl ~ host, anthro.host, function(x) max(x))

median.par <- aggregate(ghsl ~ hostpar, anthro.par, function(x) median(x))
min.par <- aggregate(ghsl ~ hostpar, anthro.par, function(x) min(x))
max.par <- aggregate(ghsl ~ hostpar, anthro.par, function(x) max(x))

median.par$host <- as.character(unlist(as.data.frame(strsplit(median.par$hostpar, "_"))[2,]))
min.par$host <- as.character(unlist(as.data.frame(strsplit(min.par$hostpar, "_"))[2,]))
max.par$host <- as.character(unlist(as.data.frame(strsplit(max.par$hostpar, "_"))[2,]))

median.par.host <- merge(median.par, median.host, by.x="host", by.y="host")
min.par.host <- merge(min.par, min.host, by.x="host", by.y="host")
max.par.host <- merge(max.par, max.host, by.x="host", by.y="host")

median.par.host$prop <- median.par.host$ghsl.x - median.par.host$ghsl.y
hist(median.par.host$prop, main="Median ghsl of parasite sample\nlocations above or below median\nghsl of host")

min.par.host$prop <- (min.par.host$ghsl.x-min.par.host$ghsl.y) / (max.par.host$ghsl.y - min.par.host$ghsl.y)
hist(min.par.host$prop, main="Proximity of min ghsl at parasite\nsample locations to min\nghsl for host.\n0 is at min")

max.par.host$prop <- max.par.host$ghsl.x/max.par.host$ghsl.y
hist(max.par.host$prop, main="Proximity of max ghsl at parasite\nsample locations to max\nghsl for host.\n1 is at max")

##### Niche position - raw #####
nichepos.far <- nichepos.cent <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("par", "pos", "hostpar"))

for(i in names(dat)) {
  if(length(dat[[i]])>1) {
    dist <- dat[[i]]$pardist.cor[,c("ParasiteCorrectedName", "distance")]
    
    cent <- aggregate(distance ~ ParasiteCorrectedName, dist, function(x) min(x))
    cent[,"host-par"] <- paste0(cent$ParasiteCorrectedName, "-", i)
    colnames(cent) <- colnames(nichepos.cent)
    nichepos.cent <- rbind(nichepos.cent, cent)
    
    far <- aggregate(distance ~ ParasiteCorrectedName, dist, function(x) max(x))
    far[,"host-par"] <- paste0(far$ParasiteCorrectedName, "-", i)
    colnames(far) <- colnames(nichepos.far)
    nichepos.far <- rbind(nichepos.far, far)
  }
}  
nichepos.farcent <- merge(nichepos.far, nichepos.cent, by.x=c("par", "hostpar"), by.y=c("par", "hostpar"))
nichepos.farcent <- nichepos.farcent[order(nichepos.farcent$pos.x),]
nichepos.farcent$y <- c(1:nrow(nichepos.farcent))

ggplot() +
  geom_segment(data = nichepos.farcent, aes(x=pos.x, y=y, xend=pos.y, yend=y), size = 1) +
  labs(x="Parasite sampling distance from host centroid", y="", title="Ordered by closest to centroid")

nichepos.farcent$pos.rng <- nichepos.farcent$pos.y - nichepos.farcent$pos.x
nichepos.farcent <- nichepos.farcent[order(nichepos.farcent$pos.rng),]
nichepos.farcent$y <- c(1:nrow(nichepos.farcent))

ggplot() +
  geom_segment(data = nichepos.farcent, aes(x=pos.x, y=y, xend=pos.y, yend=y), size = 1) +
  labs(x="Parasite sampling distance from host centroid", y="", title="Ordered by range of sampling positions")

##### Niche position - relative to host niche breadth #####
nichepos.far <- nichepos.cent <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("par", "pos", "hostpar"))

for(i in names(dat)) {
  if(length(dat[[i]])>1) {
    dist <- dat[[i]]$pardist.cor[,c("ParasiteCorrectedName", "distance")]
    
    cent <- aggregate(distance ~ ParasiteCorrectedName, dist, function(x) min(x))
    cent[,"host-par"] <- paste0(cent$ParasiteCorrectedName, "-", i)
    colnames(cent) <- colnames(nichepos.cent)
    nichepos.cent <- rbind(nichepos.cent, cent)
    
    far <- aggregate(distance ~ ParasiteCorrectedName, dist, function(x) max(x))
    far[,"host-par"] <- paste0(far$ParasiteCorrectedName, "-", i)
    colnames(far) <- colnames(nichepos.far)
    nichepos.far <- rbind(nichepos.far, far)
  }
}  
nichepos.farcent <- merge(nichepos.far, nichepos.cent, by.x=c("par", "hostpar"), by.y=c("par", "hostpar"))
nichepos.farcent <- nichepos.farcent[order(nichepos.farcent$pos.x),]
nichepos.farcent$y <- c(1:nrow(nichepos.farcent))

ggplot() +
  geom_segment(data = nichepos.farcent, aes(x=pos.x, y=y, xend=pos.y, yend=y), size = 1) +
  labs(x="Parasite sampling distance from host centroid", y="", title="Ordered by closest to centroid")

nichepos.farcent$pos.rng <- nichepos.farcent$pos.y - nichepos.farcent$pos.x
nichepos.farcent <- nichepos.farcent[order(nichepos.farcent$pos.rng),]
nichepos.farcent$y <- c(1:nrow(nichepos.farcent))

ggplot() +
  geom_segment(data = nichepos.farcent, aes(x=pos.x, y=y, xend=pos.y, yend=y), size = 1) +
  labs(x="Parasite sampling distance from host centroid", y="", title="Ordered by range of sampling positions")

