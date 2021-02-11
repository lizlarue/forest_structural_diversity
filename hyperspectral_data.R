#hyperspectral data

#install.packages("BiocManager")
#BiocManager::install("rhdf5")

library(rhdf5)
library(raster)
library(rgdal)
library(maps)


#When we refer to bands in this tutorial, we will note the band numbers for this example dataset, which are different from NEON production data. To convert a band number (b) from this example data subset to the equivalent band in a full NEON hyperspectral file (b'), use the following equation: b' = 1+4*(b-1).

wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

f <- paste0(wd,"NEON_hyperspectral_tutorial_example_subset.h5")

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

# Next, we read the different dimensions
nRows <- reflInfo$Dimensions[1]
nCols <- reflInfo$Dimensions[2]
nBands <- reflInfo$Dimensions[3]

nRows
nCols
nBands


# Extract or "slice" data for band 9 from the HDF5 file
b9 <- h5read(f,"/SJER/Reflectance/Reflectance_Data",index=list(9,1:nCols,1:nRows)) 

# what type of object is b9?
class(b9)


# convert from array to matrix by selecting only the first band
b9 <- b9[1,,]

# check it
class(b9)


# look at the metadata for the reflectance dataset
h5readAttributes(f,"/SJER/Reflectance/Reflectance_Data")

# plot the image
image(b9)

# oh, that is hard to visually interpret.
# what happens if we plot a log of the data?
image(log(b9))

# Plot range of reflectance values as a histogram to view range
# and distribution of values.
hist(b9,breaks=40,col="darkmagenta")

# View values between 0 and 5000
hist(b9,breaks=40,col="darkmagenta",xlim = c(0, 5000))

# View higher values
hist(b9, breaks=40,col="darkmagenta",xlim = c(5000, 15000),ylim=c(0,100))

#The data scale factor in the metadata tells us to divide all reflectance values by 10,000. Thus, a value of 5,000 equates to a reflectance value of 0.50. 

# there is a no data value in our raster - let's define it
myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)
myNoDataValue

# set all values equal to -9999 to NA
b9[b9 == myNoDataValue] <- NA

# plot the image now
image(b9)

image(log(b9))

# We need to transpose x and y values in order for our 
# final image to plot properly
b9 <- t(b9)
image(log(b9), main="Transposed Image")



# Extract the EPSG from the h5 dataset
myEPSG <- h5read(f, "/SJER/Reflectance/Metadata/Coordinate_System/EPSG Code")

# convert the EPSG code to a CRS string
myCRS <- crs(paste0("+init=epsg:",myEPSG))

# define final raster with projection info 
# note that capitalization will throw errors on a MAC.
# if UTM is all caps it might cause an error!
b9r <- raster(b9, 
              crs=myCRS)

# view the raster attributes
b9r

# let's have a look at our properly oriented raster. Take note of the 
# coordinates on the x and y axis.

image(log(b9r), 
      xlab = "UTM Easting", 
      ylab = "UTM Northing",
      main = "Properly Oriented Raster")


# Grab the UTM coordinates of the spatial extent
xMin <- reflInfo$Spatial_Extent_meters[1]
xMax <- reflInfo$Spatial_Extent_meters[2]
yMin <- reflInfo$Spatial_Extent_meters[3]
yMax <- reflInfo$Spatial_Extent_meters[4]

# define the extent (left, right, top, bottom)
rasExt <- extent(xMin,xMax,yMin,yMax)
rasExt


# assign the spatial extent to the raster
extent(b9r) <- rasExt

# look at raster attributes
b9r

# let's change the colors of our raster and adjust the zlims 
col <- terrain.colors(25)

image(b9r,  
      xlab = "UTM Easting", 
      ylab = "UTM Northing",
      main= "Raster w Custom Colors",
      col=col, 
      zlim=c(0,3000))

#We've now created a raster from band 9 reflectance data. We can export the data as a raster, using the writeRaster command.
# write out the raster as a geotiff
writeRaster(b9r,
            file=paste0(wd,"band9.tif"),
            format="GTiff",
            overwrite=TRUE)

# It's always good practice to close the H5 connection before moving on!
# close the H5 file
H5close()





####################################################################################

#The function band2Rast slices a band of data from the HDF5 file, and extracts the reflectance. It them converts the data to a matrix, converts it to a raster and returns a spatially corrected raster for the specified band.


# file: the hdf file
# band: the band you want to process
# returns: a matrix containing the reflectance data for the specific band

band2Raster <- function(file, band, noDataValue, extent, CRS){
  # first, read in the raster
  out <- h5read(file,"/SJER/Reflectance/Reflectance_Data",index=list(band,NULL,NULL))
  # Convert from array to matrix
  out <- (out[1,,])
  # transpose data to fix flipped row and column order 
  # depending upon how your data are formatted you might not have to perform this
  # step.
  out <- t(out)
  # assign data ignore values to NA
  # note, you might chose to assign values of 15000 to NA
  out[out == myNoDataValue] <- NA
  
  # turn the out object into a raster
  outr <- raster(out,crs=CRS)
  
  # assign the extents to the raster
  extent(outr) <- extent
  
  # return the raster object
  return(outr)
}


#Let's start with a typical RGB (red, green, blue) combination. We will use bands 14, 9, and 4 (bands 58, 34, and 19 in a full NEON hyperspectral dataset).

# create a list of the bands we want in our stack
rgb <- list(14,9,4) #list(58,34,19) when using full NEON hyperspectral dataset

# lapply tells R to apply the function to each element in the list
rgb_rast <- lapply(rgb,FUN=band2Raster, file = f,
                   noDataValue=myNoDataValue, 
                   extent=rasExt,
                   CRS=myCRS)

# check out the properties or rgb_rast
# note that it displays properties of 3 rasters.
rgb_rast


# finally, create a raster stack from our list of rasters
rgbStack <- stack(rgb_rast)

# Create a list of band names
bandNames <- paste("Band_",unlist(rgb),sep="")

# set the rasterStack's names equal to the list of bandNames created above
names(rgbStack) <- bandNames

# check properties of the raster list - note the band names
rgbStack

# scale the data as specified in the reflInfo$Scale Factor
rgbStack <- rgbStack/as.integer(reflInfo$Scale_Factor)

# plot one raster in the stack to make sure things look OK.
plot(rgbStack$Band_14, main="Band 14")

# change the colors of our raster 
myCol <- terrain.colors(25)
image(rgbStack$Band_14, main="Band 14", col=myCol)

# adjust the zlims or the stretch of the image
myCol <- terrain.colors(25)
image(rgbStack$Band_14, main="Band 14", col=myCol, zlim = c(0,.5))

# try a different color palette
myCol <- topo.colors(15, alpha = 1)
image(rgbStack$Band_14, main="Band 14", col=myCol, zlim=c(0,.5))

# create a 3 band RGB image
plotRGB(rgbStack,
        r=1,g=2,b=3,
        stretch = "lin")

# write out final raster    
# note: if you set overwrite to TRUE, then you will overwite or lose the older
# version of the tif file! Keep this in mind.
writeRaster(rgbStack, file=paste0(wd,"NEON_hyperspectral_tutorial_example_RGB_stack_image.tif"), format="GTiff", overwrite=TRUE)


#Use different band combinations to create other "RGB" images. Suggested band combinations are below for use with the full NEON hyperspectral reflectance datasets (for this example dataset, divide the band number by 4 and round to the nearest whole number):
  
#Color Infrared/False Color: rgb (90,34,19)
#SWIR, NIR, Red Band: rgb (152,90,58)
#False Color: rgb (363,246,55)


# Calculate NDVI
# select bands to use in calculation (red, NIR)
ndvi_bands <- c(16,24) #bands c(58,90) in full NEON hyperspectral dataset

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

###
#Calculate EVI using the following formula : EVI<- 2.5 * ((b4-b3) / (b4 + 6 * b3- 7.5*b1 + 1))

#Calculate Normalized Difference Nitrogen Index (NDNI) using the following equation: log(1/p1510)-log(1/p1680)/ log(1/p1510)+log(1/p1680)

#Explore the bands in the hyperspectral data. What happens if you average reflectance values across multiple red and NIR bands and then calculate NDVI?
###


####################################################################################



