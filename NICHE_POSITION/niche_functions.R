## Written by Olivier Broennimann and Blaise Petitpierre. Departement of Ecology and Evolution (DEE). 
## University of Lausanne. Switzerland. April 2012.

##### Modified by Regan Early
##### Modified by Margaret Bolton
##### Probably modified by Henry Hakkinen at some point too!

##
## DESCRIPTION
##
## functions to perform measures of niche overlap and niche equivalency/similarity tests as described in Broennimann et al. (submitted)
## 
## list of functions:
##
## grid.clim(glob,glob1,sp,R,th.sp,th.env) 
## use the scores of an ordination (or SDM predictions) and create a grid z of RxR pixels 
## (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities
## Only scores of one, or two dimensions can be used 
## sp= scores for the occurrences of the species in the ordination, glob = scores for the whole studies areas, glob 1 = scores for the range of sp 
## R= resolution of the grid, th.sp=quantile of species densitie at species occurences used as a threshold to exclude low species density values, 
## th.env=quantile of environmental densitie at all study sites used as a threshold to exclude low environmental density values
##
## niche.overlap(z1,z2,cor)
## calculate the overlap metrics D and I (see Warren et al 2008) based on two species occurrence density grids z1 and z2 created by grid.clim
## cor=T correct occurrence densities of each species by the prevalence of the environments in their range
##
## niche.equivalency.test(z1,z2,rep)
## runs niche equivalency test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim and z2.sim.
## z1.sim and z2.sim are built from random reallocations of occurences of z1 and z2
## rep is the number of iterations. This is a one-sided version of the test of Broennimann et al. 2012 (alternative hypothese = greater)
##
## niche.similarity.test(z1,z2,rep)
## runs niche similarity test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches (z2.sim) in available in the range of z2 (z2$Z) 
## z2.sim have the same patterns as z2 but their center are randomly translatated in the availabe z2$Z space and weighted by z2$Z densities
## rep is the number of iterations. This is a one-sided version of the test of Broennimann et al. 2012 (alternative hypothese = greater)
##
## plot.niche(z,title,name.axis1,name.axis2)
## plot a niche z created by grid.clim. title,name.axis1 and name.axis2 are strings for the legend of the plot
##
## plot.contrib(contrib,eigen)
## plot the contribution of the initial variables to the analysis. Typically the eigen vectors and eigen values in ordinations
##
## plot.overlap.test(x,type,title)
## plot an histogram of observed and randomly simulated overlaps, with p-values of equivalency and similarity tests. 
## x must be an object created by niche.similarity.test or niche.equivalency.test.
## type is either "D" or "I". title is the title of the plot
##
## dynamic.index(z1,z2,intersection=NA)
## calculate niche expansion, stability and unfilling
## z1 : gridclim object for the native distribution
## z2 : gridclim object for the invaded range
## intersection : quantile of the environmental density used to remove marginal climates. 
## If intersection = NA, analysis is performed on the whole environmental extent (native and invaded)
## If intersection = 0, analysis is performed at the intersection between native and invaded range
## If intersection = 0.05, analysis is performed at the intersection of the 5th quantile of both native and invaded environmental densities 
## etc...
##
## plot.niche.dyn(z1,z2,quant,title,interest,colz1,colz2,colinter,colZ1=,colZ2=)
## plot niche categories and species density
## z1 : gridclim object for the native distribution
## z2 : gridclim object for the invaded range
## quant : quantile of the environmental density used to remove marginal climates.
## title : title of the figure
## interest : choose which density to plot. If interest=1 plot native density, if interest=2 plot invasive density
## colz1 : color used to depict unfilling area
## colz2 : color used to depict expansion area
## colinter : color used to depict overlap area
## colZ1 : color used to delimit the native extent
## colZ2 : color used to delimit the invaded extent
##
## pts2img <- function(pts,extent)
## convert plots into image
## pts = points coordinates (2 columns)
## extent : grid intervals (2 columns)
##
## fun.arrows(sp1,sp2,clim1,clim2)
## draw arrows linking the centroid of the native and inasive distribution (continuous line) and between native and invaded extent (dashed line)
## sp1 : scores of the species native distribution along the the 2 first axes of the PCA
## sp2 : scores of the species invasive distribution along the the 2 first axes of the PCA
## clim1 : scores of the entire native extent along the the 2 first axes of the PCA
## clim2 : scores of the entire invaded extent along the the 2 first axes of the PCA


##################################################################################################

## Examples for test run
# glob <- scores.env
# sp <- scores.occ
# gmpd <- scores.gmpd
# R <- 100
# th.sp <- 0
# th.env <- 0

grid.clim <- function(scores.env, scores.occ, scores.gmpd, R, th.sp, th.env){
  l <- list()
  glob <- as.matrix(scores.env); sp <- as.matrix(scores.occ)
  if(ncol(glob) == 2){ #if scores in two dimensions (e.g. PCA)
    
    library(adehabitatMA)
    library(adehabitatHR)
    xmin <- min(glob[, 1]); xmax <- max(glob[, 1]); ymin <- min(glob[, 2]); ymax <- max(glob[, 2]) ## data preparation
    globr <- data.frame(cbind((glob[, 1] - xmin)/abs(xmax - xmin), (glob[, 2] - ymin)/abs(ymax - ymin))) ## Standardise the PCA scores of each env grid-cell between 0 and 1 
    spr <- data.frame(cbind((sp[, 1] - xmin)/abs(xmax - xmin), (sp[, 2] - ymin)/abs(ymax - ymin))) 	## Standardise the PCA scores at species locations to match the env standardisation
    
    ### Make a spatial pixels data frame (grid) that encompasses the standardised environmental and species location scores
    mask.xy <- expand.grid(x = seq(0.01, 1, by = 0.01), y = seq(0.01, 1, by = 0.01))
    coordinates(mask.xy) <- ~ x + y
    gridded(mask.xy) <- TRUE
    
    ### Calculate the density of occurrences in a grid of RxR pixels along the score gradients
    ## The Utilization Distribution (UD) is the bivariate function giving the probability density that a species is found within a grid-cell, based on it's occupancy of environmental conditions.
    sp.dens <- kernelUD(SpatialPoints(spr[,1:2]), h = "href", grid = mask.xy, kern = "bivnorm")	## s4 object, length is the number of cells in the grid.				
    # using a gaussian kernel density function, with RxR bins.
    #    sp.dens <- asc2spixdf(sp.dens[[1]]$UD)							# data manipulation, no longer need this 
    #sp.dens$var[sp.dens$var>0 & sp.dens$var<1] <- 0
    
    ### Calculate the density of environmental conditions in a grid of RxR pixels along the score gradients
    glob.dens <- kernelUD(SpatialPoints(globr[, 1:2]), h = "href", grid = mask.xy, kern = "bivnorm") ## s4 object, length is the number of cells in the grid.
    #    glob.dens <- asc2spixdf(glob.dens[[1]]$UD)
    #glob.dens$var[glob.dens$var < 1 & glob.dens$var > 0] <- 0
    
    x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), length.out = R) ## The locations of the grid-cell breaks on PCA axis 1
    y <- seq(from = min(glob[, 2]), to = max(glob[, 2]), length.out = R) ## The locations of the grid-cell breaks on PCA axis 2
    
    ### Make scaled matrices of species and environmental densities.
    z <- rotate2(matrix((sp.dens@data*nrow(sp)/sum(sp.dens@data))[, 1], nrow = R, ncol = R, byrow = F))	## Rescale species' occurrence density to the number of grid-cells in which the sp occurs. Multiplies the densities of species' occurrences in each grid-cell by the number of grid-cells the species occupies, then divides by the total species density summed across all grid-cells. rotate2 is a function in this source file (below)
    Z <- rotate2(matrix((glob.dens@data*nrow(glob)/sum(glob.dens@data))[, 1], nrow = R, ncol = R, byrow = F))	## Rescale density of environmental conditions to the number of all sites with environmental data. Multiplies the densities of sites with conditions that match each grid-cell by the number of sites with env data, then divides by the total density of all env conditions summed across all grid-cells. rotate2 is a function in this source file (below)
    
    ### Remove infinitesimally small number generated by kernel density function
    spr <- pts2img(sp, cbind(x, y))
    globr <- pts2img(glob, cbind(x, y))
    
    z.th <- quantile(as.vector(z[which(spr == 1)]), th.sp) ## Identifies the proportion of the data with density below the specified threshold
    Z.th <- quantile(as.vector(Z[which(globr[] == 1)]), th.env)  
    z[z < z.th] <- 0 					#z[z<max(z)/(nrow(sp)/6)] <- 0 or z[z<(nrow(sp)^(1/6)*0.005)] <- 0 				
    Z[Z < Z.th] <- 0     # [Z<(nrow(glob)^(1/6)*0.005)]
    
    ### Correct species density for environmental prevalence, and scale for comparison between species.
    z.uncor <- z/max(z)	## Rescale uncorrected species densities between [0:1] for comparison with other species  
    w <- z.uncor 
    w[w > 0] <- 1 ## Identify the grid-cells where the species has a non-trivial density
    z <- z/Z	## Correct for environment prevalence
    z[is.na(z)] <- 0	## Remove situations where density of both species and environment is 0
    z[z == "Inf"] <- 0 ## Remove n/0 situations
    z.cor <- z/max(z)	## Rescale corrected densities between [0:1] for comparison with other species
    
    ### Calculate the distance between parasite sampling locations and the centre of the host species' distribution in environmental space.
    gmpdr <- vals2img(scores.gmpd, cbind(x, y)) 
    gmpdr.xy <- as.matrix(gmpdr[, c(1, 2)])
    
    gmpd.uncor <- distangles(gmpdr.xy, which(z.uncor == 1, arr.ind = T)) ## Distance and angle of parasite locations from the grid-cell where the uncorrected species density is highest, i.e. 1. Distangles is a function in this source file (below)
    gmpd.uncor <- cbind(gmpd.uncor, z.uncor[gmpdr.xy], gmpdr[, 3:ncol(gmpdr)]) ## Makes a nice dataframe containing the key information
    names(gmpd.uncor)[which(names(gmpd.uncor) == "z.uncor[gmpdr.xy]")] <- "zp.uncor"
    
    gmpd.cor <- distangles(gmpdr.xy, which(z.cor == 1, arr.ind = T)) ## Distance and angle of parasite locations from the grid-cell where the environmentally-corrected species density is highest, i.e. 1. Distangles is a function in this source file (below)
    gmpd.cor <- cbind(gmpd.cor, z.cor[gmpdr.xy], gmpdr[, 3:ncol(gmpdr)]) ## Makes a nice dataframe containing the key information
    names(gmpd.cor)[which(names(gmpd.cor) == "z.cor[gmpdr.xy]")] <- "zp.cor"
    
    ## Positions of the host range centroids
    zp.uncor <- c(x[which((z.uncor) == 1,arr.ind = T)[, 1]], y[which((z.uncor) == 1, arr.ind = T)[, 2]]) ## The xy locations of the centre of the host species' uncorrected distribution
    zp.cor <- c(x[which((z.cor) == 1, arr.ind = T)[, 1]], y[which((z.cor) == 1, arr.ind = T)[, 2]]) ## The xy locations of the centre of the host species' environmentally-corrected distribution
    
    ## Distance between the centroid and the most distant of the host species' range
    gbif <- vals2img(scores.occ, cbind(x, y)) 
    gbif.xy <- as.matrix(gbif[, c(1, 2)])
    
    gbif.uncor <- distangles(gbif.xy, which(z.uncor == 1, arr.ind = T)) ## Distance and angle of host locations from the grid-cell where the environmentally-corrected host density is highest, i.e. 1. Distangles is a function in this source file (below)
    gbif.uncor <- cbind(gbif.uncor, z.uncor[gbif.xy], gbif[, 3:ncol(gbif)]) ## Makes a nice dataframe containing the key information
    names(gbif.uncor)[which(names(gbif.uncor) == "z.uncor[gbif.xy]")] <- "zp.uncor"
    
    gbif.cor <- distangles(gbif.xy, which(z.cor == 1, arr.ind = T)) ## Distance and angle of host locations from the grid-cell where the environmentally-corrected host density is highest, i.e. 1. Distangles is a function in this source file (below)
    gbif.cor <- cbind(gbif.cor, z.cor[gbif.xy], gbif[, 3:ncol(gbif)]) ## Makes a nice dataframe containing the key information
    names(gbif.cor)[which(names(gbif.cor) == "z.cor[gbif.xy]")] <- "zp.cor"

    l$x <- x; l$y <- y; l$hostdens.uncor <- z.uncor; l$hostdens.cor <- z.cor; l$env.dens <- Z; l$scores.env <- glob; l$scores.occ <- sp; l$pardist.uncor <- gmpd.uncor; l$pardist.cor <- gmpd.cor; l$occdist.uncor <- gbif.uncor; l$occdist.cor <- gbif.cor; l$hostocc <- w; l$hostcent.uncor <- zp.uncor; l$hostcent.cor <- zp.cor
  }
  return(l)
}

dists <- function(locs, locus){
  locus <- matrix(locus, ncol = 2)
  locs <- as.matrix(locs)
  D <- sqrt((locs[, 1] - locus[, 1])^2+(locs[, 2] - locus[, 2])^2) ## Pythagoras
  return(D)
}

distangles <- function(locs, locus){
  D <- dists(locs, locus)
  locus <- matrix(locus, ncol = 2)
  dot.prods <- locs[, 1]*locus[, 1] + locs[, 2]*locus[, 2]
  norms.x <- dists(locs = matrix(rep(0, length(locs)), nrow = nrow(locs)), locus = locs)
  norms.y <- dists(locs = matrix(rep(0, length(locus)), nrow = nrow(locus)), locus = locus)
  thetas <- acos(dot.prods / (norms.x * norms.y))
  as.numeric(thetas)
  DA <- cbind(locs, D, thetas)
  colnames(DA) <- c("x", "y", "distance", "angle")
  return(DA)
}
#distangles(matrix(1:8,4),c(1,2))

rotate2 <- function(mat) t(mat[nrow(mat):1,,drop=FALSE])

twiddle <- function(mat,way){
  if (way == "rotateclock") t(mat[nrow(mat):1,, drop = FALSE])
  if (way == "rotateanti") t(mat[,ncol(mat):1, drop = FALSE])
  if (way == "fliphor") mat[nrow(mat):1,, drop = FALSE]
  if (way == "flipver") mat[,ncol(mat):1, drop = FALSE]
  if (way == "transpose") t(mat)
}

pts2img <- function(pts, extent){
  img <- matrix(0,nrow = nrow(extent), ncol = nrow(extent))
  x <- findInterval(pts[, 1], extent[, 1])
  y <- findInterval(pts[, 2], extent[, 2])
  xy <- cbind(x, y)
  img[xy] <- 1
  return(img)
}


vals2img <- function(pts, extent){
  x <- findInterval(pts[, 1], extent[, 1])
  y <- findInterval(pts[, 2], extent[, 2])
  img.vals <- cbind(x, y, pts)
  return(img.vals)
}

##################################################################################################

niche.overlap <- function(z1, z2, cor){ 
  
  # z1 = species 1 occurrence density grid created by grid.clim
  # z2 = species 2 occurrence density grid created by grid.clim
  # cor=T correct occurrence densities of each species by the prevalence of the environments in their range
  
  l <- list()
  
  if(cor == F){		
    p1 <- z1$z.uncor/sum(z1$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
    p2 <- z2$z.uncor/sum(z2$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
  }
  
  if(cor == T){
    p1 <- z1$z.cor/sum(z1$z.cor)	# rescale occurence densities so that the sum of densities is the same for both species
    p2 <- z2$z.cor/sum(z2$z.cor)	# rescale occurence densities so that the sum of densities is the same for both species
  }
  
  D <- 1 - (0.5*(sum(abs(p1 - p2))))				# overlap metric D
  I <- 1 - (0.5*(sqrt(sum((sqrt(p1) - sqrt(p2))^2))))	# overlap metric I
  l$D <- D
  l$I <- I
  return(l)
}

##################################################################################################

niche.equivalency.test <- function(z1, z2, rep){
  
  R <- length(z1$x)
  l <- list()
  x11(2, 2, pointsize = 12); par(mar = c(0,0,0,0));
  
  obs.o <- niche.overlap(z1, z2, cor = T)									#observed niche overlap
  sim.o <- data.frame(matrix(nrow = rep, ncol = 2))							#empty list of random niche overlap
  names(sim.o) <- c("D", "I")
  for (i in 1:rep){
    plot.new(); text(0.5, 0.5, paste("runs to go:", rep - i + 1))
    
    if(is.null(z1$y)){ #overlap on one axis
      
      occ1.sim <- sample(z1$x, size = nrow(z1$sp), replace = T, prob = z1$z.cor) 		#random sampling of occurrences following the corrected densities distribution
      occ2.sim <- sample(z2$x, size = nrow(z2$sp), replace = T, prob = z2$z.cor)
      
      occ.pool <- c(occ1.sim, occ2.sim)   # pool of random occurrences
      rand.row <- sample(1:length(occ.pool), length(occ1.sim), replace = T) 		# random reallocation of occurrences to datasets
      sp1.sim <- occ.pool[rand.row]
      sp2.sim <- occ.pool[-rand.row]
      
      z1.sim <- grid.clim(z1$glob, z1$glob1, data.frame(sp1.sim), R) # gridding
      z2.sim <- grid.clim(z2$glob, z2$glob1, data.frame(sp2.sim), R)	
    }
    
    if(!is.null(z1$y)){ 										#overlap on two axes
      coordinates <- which(z1$z.cor > 0, arr.ind = T)						# array of cell coordinates ((1,1),(1,2)...)
      weight <- z1$z.cor[z1$z.cor > 0] 								#densities in the same format as cells
      coordinates.sim1 <- coordinates[sample(1:nrow(coordinates), size = nrow(z1$sp), replace = T, prob = weight),] #random sampling of coordinates following z1$z.cor distribution
      occ1.sim <- cbind(z1$x[coordinates.sim1[, 1]], z1$y[coordinates.sim1[, 2]]) 	# random occurrences following the corrected densities distribution
      
      coordinates <- which(z2$z.cor > 0, arr.ind = T)						# array of cell coordinates ((1,1),(1,2)...)
      weight <- z2$z.cor[z2$z.cor > 0] 								#densities in the same format as cells
      coordinates.sim2 <- coordinates[sample(1:nrow(coordinates), size = nrow(z2$sp), replace = T, prob = weight),] #random sampling of coordinates following z1$z.cor distribution
      occ2.sim <- cbind(z2$x[coordinates.sim2[, 1]], z2$y[coordinates.sim2[, 2]]) 	# random occurrences following the corrected densities distribution
      
      occ.pool <- rbind(occ1.sim, occ2.sim) # pool of random occurrences
      rand.row <- sample(1:nrow(occ.pool), nrow(occ1.sim), replace = T) 			# random reallocation of occurrences to datasets
      sp1.sim <- occ.pool[rand.row,]
      sp2.sim <- occ.pool[-rand.row,]
      
      z1.sim <- grid.clim(z1$glob, z1$glob1, data.frame(sp1.sim), R)
      z2.sim <- grid.clim(z2$glob, z2$glob1, data.frame(sp2.sim), R)	
    }
    
    o.i <- niche.overlap(z1.sim, z2.sim, cor = F)							# overlap between random and observed niches
    sim.o$D[i] <- o.i$D											# storage of overlaps
    sim.o$I[i] <- o.i$I
  }
  
  dev.off()
  l$sim <- sim.o	# storage
  l$obs <- obs.o	# storage
  l$p.D <- (sum(sim.o$D >= obs.o$D ) + 1)/(length(sim.o$D) + 1)	# storage of p-values one sided test
  l$p.I <- (sum(sim.o$I >= obs.o$D ) + 1)/(length(sim.o$I) + 1)	# storage of p-values one sided test
  return(l)
}

##################################################################################################

niche.similarity.test <- function(z1,z2,rep){
  
  R <- length(z1$x)
  x11(2, 2, pointsize = 12); par(mar=c(0, 0, 0, 0));						# countdown window
  l <- list()
  obs.o <- niche.overlap(z1, z2, cor = T) 								#observed niche overlap
  sim.o <- data.frame(matrix(nrow = rep, ncol = 2))						#empty list of random niche overlap
  names(sim.o) <- c("D", "I")
  
  for (k in 1:rep){
    plot.new(); text(0.5, 0.5, paste("similarity tests:", "\n", "runs to go:", rep - k + 1))	# countdown
    
    if(is.null(z2$y)){
      center <- which(z2$z.cor == 1, arr.ind = T) 					# define the centroid of the observed niche
      Z <- z2$Z/max(z2$Z)
      rand.center <- sample(1:R, size = 1, replace = F, prob = Z)				# randomly (weighted by environment prevalence) define the new centroid for the niche
      
      xshift <- rand.center - center							# shift on x axis
      z2.sim <- z2
      z2.sim$z.cor <- rep(0, R)								# set intial densities to 0
      for(i in 1:R){
        i.trans <- i + xshift
        if(i.trans > R|i.trans < 0)next()						# densities falling out of the env space are not considered
        z2.sim$z.cor[i.trans] <- z2$z.cor[i]					# shift of pixels
      }
      z2.sim$z.cor <- (z2$Z != 0)*1*z2.sim$z.cor 					# remove densities out of existing environments
    }
    
    if(!is.null(z2$y)){
      centroid <- which(z2$z.cor == 1, arr.ind = T)[1,] 				# define the centroid of the observed niche
      Z <- z2$Z/max(z2$Z)
      rand.centroids <- which(Z > 0, arr.ind = T)					# all pixels with existing environments in the study area
      weight <- Z[Z > 0]
      rand.centroid <- rand.centroids[sample(1:nrow(rand.centroids)
                                           , size = 1, replace = F, prob = weight),]					# randomly (weighted by environment prevalence) define the new centroid for the niche
      xshift <- rand.centroid[1] - centroid[1]					# shift on x axis
      yshift <- rand.centroid[2] - centroid[2]					# shift on y axis
      z2.sim <- z2
      z2.sim$z.cor <- matrix(rep(0, R*R), ncol = R, nrow = R)				# set intial densities to 0
      for(i in 1:R){
        for(j in 1:R){
          i.trans <- i+xshift
          j.trans <- j+yshift
          if(i.trans > R | i.trans < 0)next()					# densities falling out of the env space are not considered
          if(j.trans > R | j.trans < 0)next()
          z2.sim$z.cor[i.trans, j.trans] <- z2$z.cor[i, j]		# shift of pixels
        }
      }
      z2.sim$z.cor <- (z2$Z != 0)*1*z2.sim$z.cor 					# remove densities out of existing environments
    }
    o.i <- niche.overlap(z1, z2.sim, cor = T)							# overlap between random and observed niches
    sim.o$D[k] <- o.i$D										# storage of overlaps
    sim.o$I[k] <- o.i$I
  }
  dev.off()
  l$sim <- sim.o											# storage
  l$obs <- obs.o											# storage
  l$p.D <- (sum(sim.o$D >= obs.o$D ) + 1)/(length(sim.o$D) + 1)	# storage of p-values one sided test
  l$p.I <- (sum(sim.o$I >= obs.o$D ) + 1)/(length(sim.o$I) + 1)	# storage of p-values one sided test
  
  return(l)
}

##################################################################################################

plot.niche <- function(z,title,name.axis1="PC1",name.axis2="PC2",cor=F){
  
  if(is.null(z$y)){
    R <- length(z$x)
    x <- z$x
    xx <- sort(rep(1:length(x), 2))
    
    
    if(cor == F)y1 <- z$z.uncor/max(z$z.uncor)
    if(cor == T)y1 <- z$z.cor/max(z$z.cor)
    
    Y1 <- z$Z/max(z$Z)
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1)*2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1)*2)]
    
    plot(x, y1,type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = "grey")
    lines(x[xx], c(0, Y1[YY1], 0, 0))
  }
  
  if(!is.null(z$y)){
    if(cor==F)image(z$x,z$y, z$z.uncor, col = gray(100:0 / 100), zlim = c(0.000001, max(z$z.uncor)), xlab = name.axis1, ylab = name.axis2)
    if(cor==T)image(z$x,z$y, z$z.cor, col = gray(100:0 / 100), zlim = c(0.000001, max(z$z.cor)), xlab = name.axis1, ylab = name.axis2)
    contour(z$x, z$y, z$Z, add = T, levels = quantile(z$Z[z$Z > 0], c(0, 0.25)), drawlabels = F, lty = c(1, 2))
  }
  title(title)
}

plot.contrib <- function(contrib, eigen){
  
  if(ncol(contrib) == 1){
    h <- c(unlist(contrib))
    n <- row.names(contrib)
    barplot(h, space = 0, names.arg = n)
    title(main = "variable contribution")
  }
  if(ncol(contrib) == 2){
    s.corcircle(contrib[, 1:2]/max(abs(contrib[, 1:2])), grid = F)
    title(main = "Correlation circle", line = 2)
    title(sub = paste("Axis 1 = ", round(eigen[1]/sum(eigen)*100, 2), "%", ", axis 2 = ", round(eigen[2]/sum(eigen)*100, 2), "%"))
  }
}

plot.overlap.test <- function (x, type, title) {
  if(type == "D") {
    obs <- x$obs$D
    sim <- x$sim$D
    p <- x$p.D
  }
  if(type == "I") {
    obs <- x$obs$I
    sim <- x$sim$I
    p <- x$p.I
  }
  r0 <- c(sim, obs)
  l0 <- max(sim) - min(sim)
  w0 <- l0/(log(length(sim), base = 2) + 1)
  xlim0 <- range(r0) + c(-w0, w0)
  h0 <- hist(sim, plot = FALSE, nclass = 10)
  y0 <- max(h0$counts)
  hist(sim, plot = TRUE, nclass = 10, xlim = xlim0, col = grey(0.8), main = title, xlab = type, sub = paste("p.value = ",round(p,5)))
  lines(c(obs, obs), c(y0/2, 0), col = "red")
  points(obs, y0/2, pch = 18, cex = 2, col = "red")
  invisible()
}

# plot.niche.dyn <- function(z1,z2,quant=0,title,interest=1,colz1="#00FF0050",colz2="#FF000050",colinter="#0000FF50",colZ1="green3",colZ2="red3",name.axis1="PC1",name.axis2="PC2"){ ## Green is native and red is invasive, in my normal scheme of things

plot.niche.dyn <- function(z1, z2, quant = 0, title, interest = 1, colz1 = ntv.col, colz2 = usa.col, colinter = overlap.col, colZ1 = ntv.col.cont, colZ2 = usa.col.cont, name.axis1 = "PC1", name.axis2 = "PC2", yl = c(min(z1$y, z2$y), max(z1$y, z2$y))){ ## Blue is native and Yellow is invasive, in my normal scheme of things
  
  if(is.null(z1$y)){
    R <- length(z1$x)
    x <- z1$x
    xx <- sort(rep(1:length(x), 2))
    
    y1 <- z1$z.uncor/max(z1$z.uncor)
    #	if(cor == T)y1 <- z$z.cor/max(z$z.cor)
    Y1 <- z1$Z/max(z1$Z)
    if (quant>0){
      Y1.quant <- quantile(z1$Z[which(z1$Z > 0)], probs = seq(0, 1, quant))[2] / max(z1$Z)
    }else{
      Y1.quant <- 0
    }
    Y1.quant <- Y1-Y1.quant; Y1.quant[Y1.quant < 0] <- 0
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2,length(y1)*2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2,length(Y1)*2)]
    
    y2 <- z2$z.uncor/max(z2$z.uncor)
    #	if(cor==T)y1 <- z$z.cor/max(z$z.cor)
    Y2 <- z2$Z/max(z2$Z)
    if (quant > 0){
      Y2.quant <- quantile(z2$Z[which(z2$Z > 0)], probs = seq(0, 1, quant))[2] / max(z2$Z)
    }else{
      Y2.quant = 0
    }
    Y2.quant <- Y2 - Y2.quant; Y2.quant[Y2.quant < 0] <- 0
    yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2)*2)]
    YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2)*2)]
    
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = colz1, border = 0)
    polygon(x[xx], c(0, y2[yy2], 0, 0), col = colz2, border = 0)
    polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, min, na.exclude = T), 0, 0), col = colinter, border = 0)
    lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, lty = "dashed")
    lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, lty = "dashed")
    lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
    lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
    segments(x0 = 0, y0 = 0, x1 = max(x[xx]), y1 = 0, col = "white")
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, col = "white")
    
    seg.cat <- function(inter, cat, col.unf, col.exp, col.stab){
      if (inter[3] == 0){my.col = 0}
      if (inter[3] == 1){my.col = col.unf}
      if (inter[3] == 2){my.col = col.stab}
      if (inter[3] == -1){my.col = col.exp}
      segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2], col = my.col, lwd = 4, lty = 2)
    }
    cat <- dynamic.index(z1, z2, intersection = quant)$dyn
    inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1])
    apply(inter, 1, seg.cat, col.unf = "#00FF0050", col.exp = "#FF000050", col.stab = "#0000FF50")
    
  }
  if(!is.null(z1$y)){
    
    z <- z1$w+2*z2$w
    if(interest == 1){
      image(z1$x, z1$y, z1$z.uncor, col = gray(100:0 / 100), zlim = c(0.00001, max(z1$z.uncor)), xlab = name.axis1, ylab = name.axis2, ylim = yl)
      image(z1$x, z1$y, z, col = c("#FFFFFF00", colz1, colz2, colinter), add = T)
    }
    if(interest == 2){
      image(z2$x, z2$y, z2$z.uncor, col = gray(100:0 / 100), zlim = c(0.00001, max(z2$z.uncor)), xlab = name.axis1, ylab = name.axis2, ylim = yl)
      image(z2$x, z2$y, z, col = c("#FFFFFF00", colz1, colz2, colinter), add = T )
    }
    title(title)
    contour(z1$x, z1$y, z1$Z, add = T, levels = quantile(z1$Z[z1$Z > 0], c(0, quant)), drawlabels = F, lty = c(1, 2), col = colZ1)
    contour(z2$x, z2$y, z2$Z, add = T, levels = quantile(z2$Z[z2$Z > 0], c(0, quant)), drawlabels = F, lty = c(1, 2), col = colZ2)
  }
}



fun.arrows <- function(sp1, sp2, clim1, clim2, col = "red"){
  
  if (ncol(as.matrix(sp1)) == 2){
    arrows(median(sp1[, 1]), median(sp1[, 2]), median(sp2[, 1]), median(sp2[, 2]), col = "red", lwd = 2, length = 0.1)
    arrows(median(clim1[, 1]), median(clim1[, 2]), median(clim2[, 1]), median(clim2[, 2]), lty = "11", col = col, lwd = 2, length = 0.1)
  }else{
    arrows(median(sp1), 0.025, median(sp2), 0.025, col = "red", lwd = 2, length = 0.1)
    arrows(median(clim1), -0.025, median(clim2), -0.025, lty = "11", col = col, lwd = 2, length = 0.1)
  }
}

count <- function(interval, values){
  count <- c()
  for (i in 1:(length(interval)-1)){
    count <- c(count, length(which(values >= interval[i] & values < interval[i + 1])))
    if (i == length(interval)-1 ){count[i] <- length(which(values >= interval[i] & values <= interval[i + 1]))}
  }
  return(count)
}

##################################################################################################
### Rewritten by Regan
dynamic.index <- function(z1, z2, thresh) {
  glob1 <- z1$Z              # Native environmental extent densities
  glob2 <- z2$Z              # Invaded environmental extent densities
  if (thresh == 'xxx'){
    w1 <- z1$w         # Environmental native distribution in all climate space
    w2 <- z2$w         # Environmental invasive distribution in all climate space
    glob.pot <- glob2 # Climate space available in USA
    glob.pot[glob.pot[] > 0] <- 1 
    
  } else {
    quant.val <- quantile(glob1[glob1 > 0], probs = seq(0, 1, 0.01))[thresh + 1]     # threshold do delimit native environmental mask. Exclude climate below the given percentile
    glob1[glob1[] < quant.val] <- 0; glob1[glob1[] >= quant.val] <- 1                #  native environmental mask
    quant.val <- quantile(glob2[glob2 > 0], probs = seq(0, 1, 0.01))[thresh + 1]      # threshold do delimit invaded environmental mask
    glob2[glob2[] < quant.val] <- 0; glob2[glob2[] >= quant.val] <- 1                #  invaded environmental mask                                                                 
    glob <- glob.pot <- glob1*glob2     # delimitation of the intersection between the native and invaded extents    
    w1 <- z1$w*glob         # Environmental native distribution at the intersection 
    w2 <- z2$w*glob         # Environmental invasive distribution at the intersection 
  }  
  
  z.exp.cat <- (w1 + 2*w2)/2; z.exp.cat[z.exp.cat != 1] <- 0             #categorizing expansion pixels. Treats cells as occupied or not. A value of 1 from the first equation means that the cell is not occupied in Europe (w1+2=2), and occupied in the USA (w2=1), therefore (w1+2*w2)/2 = 2
  z.stable.cat <- (w1 + 2*w2)/3; z.stable.cat[z.stable.cat != 1] <- 0    #categorizing stable pixels
  z.res.cat <- w1 + 2*w2; z.res.cat[z.res.cat != 1] <- 0              #categorizing restriction pixels
  obs.exp <- z2$z.uncor*z.exp.cat                             #Correct for density of USA occupancy points
  obs.stab <- z2$z.uncor*z.stable.cat                         #Correct for density of USA occupancy points
  obs.res <- z1$z.uncor*z.res.cat						  #Correct for density of European occupancy points
  
  dyn <- (-1*z.exp.cat) + (2*z.stable.cat) + z.res.cat;    # draw matrix (100*100) with 3 categories of niche dynamic
  expansion.index.w <- sum(obs.exp)/sum(obs.stab + obs.exp); # expansion as a proportion of the invasive range
  stability.index.w <- sum(obs.stab)/sum(obs.stab + obs.exp) # stability as a proportion of the invasive range
  restriction.index.w <- sum(obs.res)/sum(obs.res + (z.stable.cat*z1$z.uncor)) #unfilling as a proportion of the native range
  total.occupancy <- sum(z2$z.uncor*glob.pot)
  
  ##### Potential Expansion:
  z.exp.cat.pot.x <- z.stable.cat + z.res.cat #z.exp.cat + 
  z.exp.cat.pot <- z.exp.cat.pot.x
  z.exp.cat.pot[z.exp.cat.pot.x == 1] <- 0     
  z.exp.cat.pot[z.exp.cat.pot.x == 0] <- 1  ## Probably a much more straightforward way to do this! But basically identifies all cells with no expansion, stability or unfilling 
  z.exp.cat.pot <- z.exp.cat.pot*glob.pot ## Control for analogue climate space if necessary
  exp.pot.cells <- sum(z.exp.cat.pot)
  exp.obs.cells <- sum(z.exp.cat)
  
  part <- list();part$dyn <- dyn
  part$dynamic.index.w <- c(round(c(expansion.index.w, stability.index.w, restriction.index.w, total.occupancy), 5), exp.obs.cells, exp.pot.cells)
  names(part$dynamic.index.w) <- c("expansion", "stability", "unfilling", "total occupancy", "expansion observed cells", "expansion potential cells")
  return(part)
}

#################################################################################################
### Rewritten by Regan
dynamic.index.nnd <- function(z1, z2, thresh) {
  glob1 <- z1$Z              # Invaded environmental extent densities
  glob2 <- z2$Z              # Native environmental extent densities
  if (thresh == 'xxx'){
    w1 <- z1$w         # Environmental invasive distribution in all climate space
    w2 <- z2$w         # Environmental native distribution in all climate space
    glob.pot <- glob2 # Climate space available in Europe
    glob.pot[glob.pot[] > 0] <- 1 
    
  } else {
    quant.val <- quantile(glob1[glob1 > 0], probs = seq(0, 1, 0.01))[thresh + 1]     # threshold do delimit native environmental mask. Exclude climate below the given percentile
    glob1[glob1[]<quant.val] <- 0; glob1[glob1[] >= quant.val] <- 1                #  native environmental mask
    quant.val <- quantile(glob2[glob2 > 0], probs = seq(0, 1, 0.01))[thresh + 1]      # threshold do delimit invaded environmental mask
    glob2[glob2[]<quant.val] <- 0; glob2[glob2[] >= quant.val] <- 1                #  invaded environmental mask                                                                 
    glob <- glob.pot <- glob1*glob2 # delimitation of the intersection between the native and invaded extents    
    w1 <- z1$w*glob         # Environmental invasive distribution at the intersection 
    w2 <- z2$w*glob         # Environmental native distribution at the intersection 
  }  
  
  z.exp.cat <- (w1 + 2*w2)/2; z.exp.cat[z.exp.cat != 1] <- 0             #categorizing expansion pixels. Treats cells as occupied or not. A value of 1 from the first equation means that the cell is not occupied in Europe (w1+2=2), and occupied in the USA (w2=1), therefore (w1+2*w2)/2 = 2
  z.stable.cat <- (w1 + 2*w2)/3; z.stable.cat[z.stable.cat != 1] <- 0    #categorizing stable pixels
  z.res.cat <- w1 + 2*w2; z.res.cat[z.res.cat != 1] <- 0              #categorizing restriction pixels
  
  obs.exp <- z2$z.uncor*z.exp.cat                             #Correct for density of European occupancy points
  obs.stab <- z2$z.uncor*z.stable.cat                         #Correct for density of European occupancy points
  obs.res <- z1$z.uncor*z.res.cat    					  #Correct for density of USA occupancy points
  
  dyn <- (-1*z.exp.cat) + (2*z.stable.cat) + z.res.cat;    # draw matrix (100*100) with 3 categories of niche dynamic
  nnd.propn <- sum(obs.res)/sum(obs.stab + obs.exp); # expansion as a proportion of the native range
  nnd <- sum(obs.res)
  
  exp.obs.cells <- sum(z.res.cat)
  exp.obs.cells.propn <- sum(z.res.cat) / (sum(z.stable.cat) + sum(z.exp.cat))# expansion as a proportion of the native range
  
  part <- list(); part$dyn <- dyn
  part$nnd <- round(c(nnd, nnd.propn, exp.obs.cells, exp.obs.cells.propn), 3)
  names(part$nnd) <- c('nnd', 'nnd.propn', 'nnd.observed.cells', 'nnd.observed.cells.propn')
  return(part)
}