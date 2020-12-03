#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)
library(neondiversity)

############### Set working directory ######
#set the working of the downloaded data
wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

############ Read in LiDAR data ###########
#UNDE <- readLAS(paste0(wd,"NEON_D05_UNDE_DP1_310000_5116000_classified_point_cloud_colorized.laz"))
UNDE <- readLAS(paste0(wd,"NEON_D05_UNDE_DP1_310000_5116000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 164 -drop_z_above 1463")

summary(UNDE)
#plot(UNDE)


#set center of plot based on extent
x <- ((max(UNDE$X) - min(UNDE$X))/2)+ min(UNDE$X)
y <- ((max(UNDE$Y) - min(UNDE$Y))/2)+ min(UNDE$Y)

data.200m <- lasclipRectangle(UNDE, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)
#looks a little weird

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

UNDE_structural_diversity <- structural_diversity_metrics(data.40m)



#####################################
#diversity data
devtools::insUNDE_github("admahood/neondiversity")

# load packages
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
library(neondiversity)

#no cover data at this site
coverU <- loadByProduct (dpID = "DP1.10058.001", site = 'UNDE')

coverDivU <- coverU[[2]]

unique(coverDivU$divDataType)

cover2U <- coverDivU %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2U$scientificName))

summary(cover2U$nativeStatusCode)

#subset of invasive only
inv <- cover2U %>%
  filter(nativeStatusCode=="I")

exotic_SR <-length(unique(inv$scientificName))


UNDE_table <- cbind(UNDE_structural_diversity, all_SR, exotic_SR)

UNDE_table <- UNDE_table %>%
  mutate(Site.ID = "UNDE")

UNDE_table <- UNDE_table %>%
  select(-easting, -northing)

UNDE_table <- UNDE_table %>%
  left_join(veg_types)


combo12 <- rbind(combo11, UNDE_table)
combo12

write.table(combo12, file = "prelim_results.csv", sep = ",", row.names = FALSE)

library(ggplot2)
ggplot(combo12, aes(x = mean.max.canopy.ht.aop, y = exotic_SR))+
  geom_point()

ggplot(combo12, aes(x = max.canopy.ht.aop, y = exotic_SR))+
  geom_point()

ggplot(combo12, aes(x = rumple.aop, y = exotic_SR))+
  geom_point()
