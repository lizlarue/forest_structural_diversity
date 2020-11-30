#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)

#######################################
#example data only
wd <- "/Users/rana7082/Dropbox/forest_structure_disturbance/"
setwd(wd)

TEAK <- readLAS(paste0(wd,"NEON_D17_TEAK_DP1_316000_4091000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 1694 -drop_z_above 2500")
summary(TEAK)

HARV <- readLAS(paste0(wd,"NEON_D01_HARV_DP1_727000_4702000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 150 -drop_z_above 325")

summary(HARV)
#######################################

############### Set working directory ######
#set the working of the downloaded data
wd <- "/Users/rana7082/Documents/research/forest_structural_disturbance/data"
setwd(wd)

############ Read in LiDAR data ###########

#Watch out for outlier Z points - this function also allows for the
#ability to filter outlier points well above or below the landscape
#(-drop_z_blow and -drop_z_above). See how we have done this here 
#for you.

#how to choose good z-values???
ONAQ <- readLAS(paste0(wd,"NEON_D15_ONAQ_DP1_367000_4454000_classified_point_cloud_colorized.laz"))
ONAQ <- readLAS(paste0(wd,"NEON_D15_ONAQ_DP1_367000_4454000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 1960 -drop_z_above 2650")

summary(ONAQ)
plot(ONAQ)


#set center of plot based on extent
x <- 367500 #easting 
y <- 4454500 #northing

#Cut out a 200 x 200 m buffer by adding 100 m to easting and 
#northing coordinates (x,y).
data.200m <- 
  lasclipRectangle(ONAQ,
                   xleft = (x - 100), ybottom = (y - 100),
                   xright = (x + 100), ytop = (y + 100))

#Correct for ground height using a kriging function to interpolate 
#elevation from ground points in the .laz file.
#If the function will not run, then you may need to checkfor outliers
#by adjusting the 'drop_z_' arguments when reading in the .laz files.
dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))


data.200m <- lasnormalize(data.200m, dtm)

#Will often give a warning if not all points could be corrected, 
#but visually check to see if it corrected for ground height. 
plot(data.200m)
#There's only a few uncorrected points and we'll fix these in 
#the next step. 

#Clip 20 m out from each side of the easting and northing 
#coordinates (x,y).
data.40m <- 
  lasclipRectangle(data.200m, 
                   xleft = (x - 20), ybottom = (y - 20),
                   xright = (x + 20), ytop = (y + 20))

data.40m@data$Z[data.40m@data$Z < 0] <- NA  
#This line filters out all z_vals below .5 m as we are less 
#interested in shrubs/trees. 
#You could change it to zero or another height depending on interests. 

#visualize the clipped plot point cloud
plot(data.40m) 
summary(data.40m)

############# Structural diversity metrics  ######

#GENERATE CANOPY HEIGHT MODEL (CHM) (i.e. a 1 m2 raster grid of 
#vegetations heights)
#res argument specifies pixel size in meters and dsmtin is 
#for raster interpolation
chm <- grid_canopy(data.40m, res = 1, dsmtin())  

#visualize CHM
plot(chm)


#MEAN OUTER CANOPY HEIGHT (MOCH)
#calculate MOCH, the mean CHM height value
mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 

#MAX CANOPY HEIGHT
#calculate HMAX, the maximum CHM height value
max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 

#RUMPLE
#calculate rumple, a ratio of outer canopy surface area to 
#ground surface area (1600 m^2)
rumple <- rumple_index(chm) 

#TOP RUGOSITY
#top rugosity, the standard deviation of pixel values in chm and 
#is a measure of outer canopy roughness
top.rugosity <- sd(chm@data@values, na.rm = TRUE) 

#DEEP GAPS & DEEP GAP FRACTION
#number of cells in raster (also area in m2)
cells <- length(chm@data@values) 
chm.0 <- chm
chm.0[is.na(chm.0)] <- 0 #replace NAs with zeros in CHM
#create variable for the number of deep gaps, 1 m^2 canopy gaps
zeros <- which(chm.0@data@values == 0) 
deepgaps <- length(zeros) #number of deep gaps
#deep gap fraction, the number of deep gaps in the chm relative 
#to total number of chm pixels
deepgap.fraction <- deepgaps/cells 

#COVER FRACTION
#cover fraction, the inverse of deep gap fraction
cover.fraction <- 1 - deepgap.fraction 

#HEIGHT SD
#height SD, the standard deviation of height values for all points
#in the plot point cloud
vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 

## Warning: 'lasmetrics' is deprecated.
## Use 'cloud_metrics' instead.
## See help("Deprecated")

#SD of VERTICAL SD of HEIGHT
#rasterize plot point cloud and calculate the standard deviation 
#of height values at a resolution of 1 m^2
sd.1m2 <- grid_metrics(data.40m, sd(Z), 1)
#standard deviation of the calculated standard deviations 
#from previous line 
#This is a measure of internal and external canopy complexity
sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 


#some of the next few functions won't handle NAs, so we need 
#to filter these out of a vector of Z points
Zs <- data.40m@data$Z
Zs <- Zs[!is.na(Zs)]

#ENTROPY 
#entropy, quantifies diversity & evenness of point cloud heights 
#by = 1 partitions point cloud in 1 m tall horizontal slices 
#ranges from 0-1, with 1 being more evenly distributed points 
#across the 1 m tall slices 
entro <- entropy(Zs, by = 1) 

#GAP FRACTION PROFILE 
#gap fraction profile, assesses the distribution of gaps in the 
#canopy volume 
#dz = 1 partitions point cloud in 1 m horizontal slices 
#z0 is set to a reasonable height based on the age and height of 
#the study sites 
gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3) 
#defines gap fraction profile as the average gap fraction in each 
#1 m horizontal slice assessed in the previous line
GFP.AOP <- mean(gap_frac$gf) 

#VAI
#leaf area density, assesses leaf area in the canopy volume 
#k = 0.5 is a standard extinction coefficient for foliage 
#dz = 1 partitions point cloud in 1 m horizontal slices 
#z0 is set to the same height as gap fraction profile above
LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
#vegetation area index, sum of leaf area density values for 
#all horizontal slices assessed in previous line
VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 

#VCI
#vertical complexity index, fixed normalization of entropy 
#metric calculated above
#set zmax comofortably above maximum canopy height
#by = 1 assesses the metric based on 1 m horizontal slices in 
#the canopy
VCI.AOP <- VCI(Zs, by = 1, zmax=100)


#OUTPUT CALCULATED METRICS INTO A TABLE
#creates a dataframe row, out.plot, containing plot descriptors 
#and calculated metrics
ONAQ_structural_diversity <- 
  data.frame(matrix(c(x, y, mean.max.canopy.ht, max.canopy.ht, 
                      rumple, deepgaps,deepgap.fraction,
                      cover.fraction,top.rugosity, vert.sd, 
                      sd.sd, entro,GFP.AOP, VAI.AOP, VCI.AOP),
                    ncol = 15)) 

#provide descriptive names for the calculated metrics
colnames(ONAQ_structural_diversity) <- 
  c("easting", "northing", "mean.max.canopy.ht.aop",
    "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
    "deepgap.fraction.aop","cover.fraction.aop", 
    "top.rugosity.aop", "vert.sd.aop", "sd.sd.aop", 
    "entropy.aop", "GFP.AOP.aop", "VAI.AOP.aop", "VCI.AOP.aop") 

#View the results
ONAQ_structural_diversity



#####################################
#diversity data
devtools::install_github("admahood/neondiversity")

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
coverO <- loadByProduct (dpID = "DP1.10058.001", site = 'ONAQ', check.size= TRUE)

coverDivO <- coverO[[2]]

unique(coverDivO$divDataType)

cover2O <- coverDivO %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2O$scientificName))

summary(cover2O$nativeStatusCode)

#subset of invasive only
inv <- cover2O %>%
  filter(nativeStatusCode=="I")

exotic_SR <-length(unique(inv$scientificName))


ONAQ_table <- cbind(ONAQ_structural_diversity, all_SR, exotic_SR)

ONAQ_table <- ONAQ_table %>%
  mutate(Site.ID = "ONAQ")


ONAQ_table <- ONAQ_table %>%
  select(-easting, -northing)

ONAQ_table <- ONAQ_table %>%
  left_join(veg_types)


combo2 <- rbind(combo, ONAQ_table)
combo2


######################
#soil chem for NIWO
ONAQ_soil_chem_0 <- read.csv(file = 'NEON.D15.ONAQ.DP1.00096.001.mgp_perbiogeosample.2014-06.basic.20201006T215549Z.csv')


ONAQ_soil_chem <- read.csv(file = 'NEON.D15.ONAQ.DP1.00096.001.mgp_perbiogeosample.2014-06.basic.20201006T215549Z.csv') %>%
  select(siteID, carbonTot, horizonName) %>%
  rename(Site.ID = siteID) %>%
  mutate(horizon = ifelse(horizonName == "A1" | horizonName == "A2", "A", ifelse(horizonName == "Bw" | horizonName == "Bk"| horizonName == "2Bkm", "B", "C")))

