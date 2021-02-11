#hyperspectral data

#install.packages("BiocManager")
#BiocManager::install("rhdf5")

library(rhdf5)
library(raster)
library(rgdal)


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

