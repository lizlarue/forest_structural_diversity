#structural diversity tutorial
#https://www.neonscience.org/structural-diversity-discrete-return

library(lidR)
library(gstat)
library(neondiversity)
library(rhdf5)
library(raster)
library(rgdal)
library(maps)
library(plyr)
library(reshape2)
library(ggplot2)

############### Set working directory ######
#set the working of the downloaded data
wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

############ Read in LiDAR data ###########
#SJER <- readLAS(paste0(wd,"NEON_D17_SJER_DP1_255000_4110000_classified_point_cloud_colorized.laz"))
SJER <- readLAS(paste0(wd,"NEON_D17_SJER_DP1_255000_4110000_classified_point_cloud_colorized.laz"),
                filter = "-drop_z_below 187 -drop_z_above 616")
summary(SJER)
#plot(SJER)


#set center of plot based on extent
x <- ((max(SJER$X) - min(SJER$X))/2)+ min(SJER$X)
y <- ((max(SJER$Y) - min(SJER$Y))/2)+ min(SJER$Y)

#x <- 255500 #easting 
#y <- 4110500 #northing

data.200m <- lasclipRectangle(SJER, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))

data.200m <- lasnormalize(data.200m, dtm)

data.40m <- lasclipRectangle(data.200m, 
                             xleft = (x - 20), ybottom = (y - 20),
                             xright = (x + 20), ytop = (y + 20))
data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
plot(data.40m)
#not sure if this looks right


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

SJER_structural_diversity <- structural_diversity_metrics(data.40m)



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
coverS <- loadByProduct (dpID = "DP1.10058.001", site = 'SJER', check.size = FALSE)

coverDivS <- coverS[[2]]

unique(coverDivS$divDataType)

cover2S <- coverDivS %>%
  filter(divDataType=="plantSpecies")

all_SR <-length(unique(cover2S$scientificName))

summary(cover2S$nativeStatusCode)

#subset of invasive only
inv <- cover2S %>%
  filter(nativeStatusCode=="I")

#total SR of exotics across all plots
exotic_SR <-length(unique(inv$scientificName))

#mean plot percent cover of exotics
exotic_cover <- inv %>%
  group_by(plotID) %>%
  summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
  summarize(exotic_cov = mean(sumz))


SJER_table <- cbind(SJER_structural_diversity, all_SR, exotic_SR, exotic_cover)

SJER_table <- SJER_table %>%
  mutate(Site.ID = "SJER")


SJER_table <- SJER_table %>%
  select(-easting, -northing)

SJER_table <- SJER_table %>%
  left_join(veg_types)


#combo6 <- rbind(combo5, SJER_table)
#combo6



#############################################
#calculate spectral reflectance as CV
#as defined here: https://www.mdpi.com/2072-4292/8/3/214/htm
f <- paste0(wd,"NEON_D17_SJER_DP3_257000_4112000_reflectance.h5")


###
#for each of the 426 bands, I need to calculate the mean reflectance and the SD reflectance across all pixels 

myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)

dat <- data.frame()

for (i in 1:426){
#extract one band
b <- h5read(f,"/SJER/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 

# set all values equal to -9999 to NA
b[b == myNoDataValue] <- NA

#calculate mean and sd
meanref <- mean(b)
SDref <- sd(b)

rowz <- cbind(i, meanref, SDref)

dat <- rbind(dat, rowz)
}

dat

dat$calc <- dat$SDref/dat$meanref

CV <- sum(dat$calc)/426


SJER_table$specCV <- CV

combo6 <- rbind(combo5, SJER_table)
combo6











##############################################
#other exploration with spectral diversity
######
# Next, we read the different dimensions
nRows <- reflInfo$Dimensions[1]
nCols <- reflInfo$Dimensions[2]
nBands <- reflInfo$Dimensions[3]

nRows
nCols
nBands
#426 bands


View(h5ls(f,all=T))



# get information about the wavelengths of this dataset
wavelengthInfo <- h5readAttributes(f,"/SJER/Reflectance/Metadata/Spectral_Data/Wavelength")
wavelengthInfo


# read in the wavelength information from the HDF5 file
wavelengths <- h5read(f,"/SJER/Reflectance/Metadata/Spectral_Data/Wavelength")
head(wavelengths)
tail(wavelengths)

# First, we need to extract the reflectance metadata:
reflInfo <- h5readAttributes(f, "/SJER/Reflectance/Reflectance_Data")
reflInfo




# convert from array to matrix
b9 <- b9[1,,]

# Extract the EPSG from the h5 dataset
myEPSG <- h5read(f, "/SJER/Reflectance/Metadata/Coordinate_System/EPSG Code")
# convert the EPSG code to a CRS string
myCRS <- crs(paste0("+init=epsg:",myEPSG))
# if there is a no data value in our raster - let's define it



# We need to transpose x and y values in order for our 
# final image to plot properly
b9 <- t(b9)

# define final raster with projection info 
b9r <- raster(b9, 
              crs=myCRS)










# Calculate NDVI
# select bands to use in calculation (red, NIR)
ndvi_bands <- c(58,90) #bands c(58,90) in full NEON hyperspectral dataset

# create raster list and then a stack using those two bands
ndvi_rast <- lapply(ndvi_bands,FUN=band2Raster, file = f,
                    noDataValue=myNoDataValue, 
                    extent=rasExt, CRS=myCRS)
ndvi_stack <- stack(ndvi_rast)

# make the names pretty
bandNDVINames <- paste("Band_",unlist(ndvi_bands),sep="")
names(ndvi_stack) <- bandNDVINames

# view the properties of the new raster stack
ndvi_stack

#calculate NDVI
NDVI <- function(x) {
  (x[,2]-x[,1])/(x[,2]+x[,1])
}
ndvi_calc <- calc(ndvi_stack,NDVI)
plot(ndvi_calc, main="NDVI for the NEON SJER Field Site")


# Now, play with breaks and colors to create a meaningful map
# add a color map with 4 colors
myCol <- rev(terrain.colors(4)) # use the 'rev()' function to put green as the highest NDVI value
# add breaks to the colormap, including lowest and highest values (4 breaks = 3 segments)
brk <- c(0, .25, .5, .75, 1)

# plot the image using breaks
plot(ndvi_calc, main="NDVI for the NEON SJER Field Site", col=myCol, breaks=brk)

