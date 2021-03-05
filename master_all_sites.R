#master for all forested sites

# load packages
library(lidR)
library(gstat)
library(neondiversity)
library(neonUtilities)
library(geoNEON)
library(dplyr, quietly=T)
library(downloader)
library(ggplot2)
library(tidyr)
library(doBy)
library(sf)
library(sp)
library(devtools)



wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

#read in data for all sites
#NIWO <- readLAS(paste0(wd,"NEON_D13_NIWO_DP1_454000_4425000_classified_point_cloud_colorized.laz"))
NIWO <- readLAS(paste0(wd,"NEON_D13_NIWO_DP1_454000_4425000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 2460 -drop_z_above 2947")

#WREF <- readLAS(paste0(wd,"NEON_D16_WREF_DP1_578000_5072000_classified_point_cloud_colorized.laz"))
WREF <- readLAS(paste0(wd,"NEON_D16_WREF_DP1_578000_5072000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 425 -drop_z_above 1082")


#loop through sites
sites <- list(NIWO, WREF)
sites <- list(NIWO, WREF, KONZ, YELL, ONAQ, SJER, CLBJ, OSBS, TALL, HARV, SCBI, UNDE, PUUM, SRER, BONA, GUAN)
#sites <- as.factor(c('NIWO', 'WREF', 'KONZ', 'YELL', 'ONAQ', 'SJER', 'CLBJ', 'OSBS', 'TALL', 'HARV', 'SCBI', 'UNDE', 'PUUM', 'SRER', 'BONA', 'GUAN'))

#dfarray <- list(df1, df2, df3)
#lapply(dfarray, function(i) plot(i$x, i$y))

for (i in sites) {

#Let's correct for elevation and measure structural diversity for SITE
x <- ((max(i$X) - min(i$X))/2)+ min(i$X)
y <- ((max(i$Y) - min(i$Y))/2)+ min(i$Y)


data.200m <- lasclipRectangle(i, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)


#Zip up all the code we previously used and write function to 
#run all 13 metrics in a single function. 
structural_diversity_metrics <- function(data.40m) {
  chm <- grid_canopy(data.40m, res = 1, dsmtin()) 
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  rumple <- rumple_index(chm) 
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) 
  deepgap.fraction <- deepgaps/cells 
  cover.fraction <- 1 - deepgap.fraction 
  vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 
  sd.1m2 <- grid_metrics(data.40m, sd(Z), 1) 
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  Zs <- data.40m@data$Z
  Zs <- Zs[!is.na(Zs)]
  entro <- entropy(Zs, by = 1) 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3)
  GFP.AOP <- mean(gap_frac$gf) 
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  out.plot <- data.frame(
    matrix(c(x, y, mean.max.canopy.ht,max.canopy.ht, 
             rumple,deepgaps, deepgap.fraction, 
             cover.fraction, top.rugosity, vert.sd, 
             sd.sd, entro, GFP.AOP, VAI.AOP,VCI.AOP),
           ncol = 15)) 
  colnames(out.plot) <- 
    c("easting", "northing", "mean.max.canopy.ht.aop",
      "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
      "deepgap.fraction.aop", "cover.fraction.aop",
      "top.rugosity.aop","vert.sd.aop","sd.sd.aop", 
      "entropy.aop", "GFP.AOP.aop",
      "VAI.AOP.aop", "VCI.AOP.aop") 
  print(out.plot)
}

i_structural_diversity <- structural_diversity_metrics(data.40m)





#####################################
#diversity data
cover <- loadByProduct (dpID = "DP1.10058.001", site = i, check.size= FALSE)

coverDiv <- cover[[2]]

unique(coverDiv$divDataType)

cover2 <- coverDiv %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2$scientificName))

summary(cover2$nativeStatusCode)


#subset of invasive only
inv <- cover2 %>%
  filter(nativeStatusCode=="I")

#total SR of exotics across all plots
exotic_SR <-length(unique(inv$scientificName))

#mean plot percent cover of exotics
exotic_cover <- inv %>%
  group_by(plotID) %>%
  summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
  summarize(exotic_cov = mean(sumz))

i_table <- cbind(i_structural_diversity, all_SR, exotic_SR, exotic_cover)

i_table <- i_table %>%
  mutate(Site.ID = i)

i_table <- i_table %>%
  select(-easting, -northing)

veg_types <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/field-sites.csv') %>%
  select(Site.ID, Dominant.NLCD.Classes)

i_table <- i_table %>%
  left_join(veg_types)

i_table

}

#############################################
#calculate spectral reflectance as CV
#as defined here: https://www.mdpi.com/2072-4292/8/3/214/htm
f <- paste0(wd,"NEON_D13_NIWO_DP3_453000_4427000_reflectance.h5")


###
#for each of the 426 bands, I need to calculate the mean reflectance and the SD reflectance across all pixels 

myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)

dat <- data.frame()

for (i in 1:426){
  #extract one band
  b <- h5read(f,"/i/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 
  
  # set all values equal to -9999 to NA
  b[b == myNoDataValue] <- NA
  
  #calculate mean and sd
  meanref <- mean(b, na.rm = TRUE)
  SDref <- sd(b, na.rm = TRUE)
  
  rowz <- cbind(i, meanref, SDref)
  
  dat <- rbind(dat, rowz)
}

dat

dat$calc <- dat$SDref/dat$meanref

CV <- sum(dat$calc)/426


i_table$specCV <- CV
