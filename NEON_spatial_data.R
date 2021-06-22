#working with NEON spatial data to spatially align field sampling plots with AOP data

library(sp)
library(sf)
library(rgdal)
library(broom)
library(ggplot2)
library(tidyverse)
library(neonUtilities)
library(devtools)

devtools::install_github("NEONScience/NEON-geolocation/geoNEON")

wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

neonDomains <- readOGR(dsn = "." , layer="NEON_Domains")


# First, add a new column termed "id" composed of the row names of the data
neonDomains@data$id <- rownames(neonDomains@data)

# Now, use tidy() to convert to a dataframe
# if you previously used fortify(), this does the same thing. 
neonDomains_points<- tidy(neonDomains, region="id")

# Finally, merge the new data with the data from our spatial object
neonDomainsDF <- merge(neonDomains_points, neonDomains@data, by = "id")

# view data structure for each variable
str(neonDomainsDF)

# plot domains
domainMap <- ggplot(neonDomainsDF) + 
  geom_map(map = neonDomainsDF,
           aes(x = long, y = lat, map_id = id),
           fill="white", color="black", size=0.3)
domainMap


#read in the data
neonSites <- read.delim("field-sites.csv", sep=",", header=T)

# view data structure for each variable
str(neonSites)

# plot the sites
neonMap <- domainMap + 
  geom_point(data=neonSites, 
             aes(x=Longitude, y=Latitude))

neonMap



# color is determined by the order that the unique values show up. Check order
unique(neonSites$Site.Type)
## [1] "Relocatable Aquatic"     "Core Terrestrial"        "Core Aquatic"            "Relocatable Terrestrial"

# add color
sitesMap <- neonMap + 
  geom_point(data=neonSites, 
             aes(x=Longitude, y=Latitude, color=Site.Type)) + 
  scale_color_manual(values=c("lightskyblue", "forest green", 
                              "blue4", "light green"),
                     name="",
                     breaks=unique(neonSites$Site.Type))
sitesMap



## load TOS plot readme
rdme <- read.delim('readme.csv', sep=',', header=T)

## View the variables
rdme[,1]
#this has lat, long as well as easting, northing







plot_locations <- readOGR(dsn=".", layer = "NEON_TOS_Plot_Centroids")

ecoregions <- readOGR(dsn=".", layer = "us_eco_l3_state_boundaries")

tot_table_plots <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/cover_by_plot.csv')




#initial plot
plot(site_locations)
#this looks like it works

#plot(ecoregions[2])

ggplot() + 
  geom_polygon(data = site_locations, aes(x = lat, y = long))
#Regions defined for each Polygons

#does not have the site polygon, obviously, but avoids the Regions problem from geom_polygon
ggplot() + 
  geom_point(data = site_locations, aes(x = lat, y = long))

summary(site_locations)
head(site_locations)


#bring in shapefile of ecoregions
ggplot() + 
  geom_polygon(data = site_locations, aes(x = lat, y = long))







###
#bring in kmz of AOP flight boundaries

file <- "data/Burrows_et_al_Nature_traj_ocean_NH1.kmz"
SST_start = readOGR(dsn = file, layer = "SST_start") 




###plot plot cover data
ggplot() + 
  geom_point(data = tot_table_plots, aes(x = decimalLatitude, y = decimalLongitude, color = Dominant.NLCD.Classes))

ggplot() + 
  geom_polygon(data = ecoregions[2], aes(x = coords, y = coords))




#extract the LIDAR AOP 1 km2 tiles that match the plot centroids specified but can only do 1 year and site at a time
#will ask if you want to download - indicate y

#filter to two sites of interest for most recent year

NIWO <- tot_table_plots %>%
  filter(siteID == "NIWO", year == 2020)
#33 obs

SOAP <- tot_table_plots %>%
  filter(siteID == "SOAP", year == 2019)
#32 obs

cord.decN <- SpatialPoints(cbind(NIWO$decimalLongitude, -NIWO$decimalLatitude), proj4string=CRS("+proj=longlat"))
cord.decS <- SpatialPoints(cbind(SOAP$decimalLongitude, -SOAP$decimalLatitude), proj4string=CRS("+proj=longlat"))

cord.UTMN <- spTransform(cord.decN, CRS("+init=epsg:32611"))
cord.UTMS <- spTransform(cord.decS, CRS("+init=epsg:32611"))


cord.UTMN <- spTransform(cord.decN, CRS('+proj=utm +zone=13 ellps=WGS84'))
cord.UTMS <- spTransform(cord.decS, CRS('+proj=utm +zone=10 ellps=WGS84'))

byTileAOP(dpID = "DP1.30003.001", site = "NIWO",  year = 2020, 
          easting = cord.UTMN$coords.x1, northing = cord.UTMN$coords.x2,
          check.size = T, buffer = 500, savepath = ".data/")
