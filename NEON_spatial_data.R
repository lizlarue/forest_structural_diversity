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




plot_centroids <- read.delim('All_NEON_TOS_Plot_Centroids_V8.csv', sep=',', header=T) %>%
  rename(
    decimalLatitude = latitude,
    decimalLongitude = longitude
  )
#this has easting and northing

#plot_points <- read.delim('All_NEON_TOS_Plot_Points_V8.csv', sep=',', header=T)
#this also has easting and northing

#bring in my data
tot_table_plots <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/cover_by_plot.csv')

#left join to add the easting and northing variables
tot_table_plots_en <- tot_table_plots %>%
  left_join(plot_centroids)

write.table(tot_table_plots_en, file = "data/tot_table_plots_en.csv", sep = ",", row.names = FALSE)


#####################################
#plot_locations <- readOGR(dsn=".", layer = "NEON_TOS_Plot_Centroids")

ecoregions <- readOGR(dsn=".", layer = "us_eco_l3_state_boundaries")
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

NIWO <- tot_table_plots_en %>%
  filter(siteID == "NIWO", year == 2020)
#33 obs

SOAP <- tot_table_plots_en %>%
  filter(siteID == "SOAP", year == 2019)
#32 obs


byTileAOP(dpID = "DP1.30003.001", site = "NIWO",  year = 2020, 
          easting = NIWO$easting, northing = NIWO$northing,
          check.size = T, buffer = 500, savepath = ".")

byTileAOP(dpID = "DP1.30003.001", site = "SOAP",  year = 2019, 
          easting = SOAP$easting, northing = SOAP$northing,
          check.size = T, buffer = 500, savepath = ".")
