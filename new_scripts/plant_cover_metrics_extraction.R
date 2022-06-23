# This code extracts plant cover data on a site and a plot level.
# Then combines plot-level plant cover data with plot-level LiDAR data.

library(lidR)
library(gstat)
library(neonUtilities)
library(dplyr)

coverdata <- read.csv("prelim_cover.csv")
veg_types <- read.csv('field-sites.csv')
plot_data_table <- read.csv("plot_data_table.csv")
lidar_data <- read.csv("lidar_structural_metrics.csv")

veg_types <- veg_types %>% 
  select(Site.ID, Dominant.NLCD.Classes) %>% 
  rename(siteID = Site.ID, Dominant_NLCD_Classes = Dominant.NLCD.Classes)

#---------------------------------------------------------------
# SITE-LEVEL DATA

unique_sitemonthyear <- unique(coverdata$sitemonthyear)

#will contain site level data
tot_table = data.frame()

for (sitemonthyear in unique_sitemonthyear) {
  siteID <- substr(sitemonthyear,1,4) 
  monthyear <- substr(sitemonthyear,5,11)
  year <- substr(monthyear, 1, 4)
  
  #filter for site/date combinations of interest
  sub <- coverdata[coverdata$sitemonthyear == sitemonthyear,]
  
  #calculating number of plots within a site
  numplots <- length(unique(sub$plotID))
  
  #total species richness
  species_richness <-length(unique(sub$scientificName))
  
  #subset of invasive only
  inv <- sub[sub$nativeStatusCode == 'I',]
  
  #total species richness of exotics across all plots
  exotic_richness <-length(unique(inv$scientificName))
  
  #calculating mean percent cover of exotic species per plot
  exotic_cover <- mean(inv$percentCover, na.rm = TRUE)
  
  i_table <- cbind(sitemonthyear, siteID, monthyear, year, numplots, species_richness, exotic_richness, exotic_cover)
  tot_table <- rbind(tot_table,i_table)
}

#join veg table data with tot_table
tot_table <- tot_table %>% left_join(veg_types)

#export table of plant cover by site 
write.table(tot_table, file = "cover_by_site.csv", sep = ",", row.names = FALSE)

#---------------------------------------------------------------
# PLOT-LEVEL DATA

#calculates total species richness for each sitemonthyear, plot combination
all_species_richness <- coverdata %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(species_richness = n_distinct(scientificName)) 

#calculates exotic species richness for each sitemonthyear, plot combo
all_exotic_richness <- coverdata %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_richness =  n_distinct(scientificName))

#calculates exotic cover for each sitemonthyear, plot combo
all_exotic_cover <- coverdata %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_cover = sum(percentCover))

#joins these tables together
closer <- left_join(all_species_richness, all_exotic_richness, by=c("sitemonthyear", "plotID")) %>%
  left_join(., all_exotic_cover, by=c("sitemonthyear", "plotID"))

plot_data_table <- plot_data_table %>% 
  select(sitemonthyear, plotID, latitude, longitude, easting, northing)

#joins these tables together to make plot level total table
tot_table_plots <- left_join(closer, plot_data_table, by=c("sitemonthyear", "plotID"))

#get rows where latitude is NA
missing_lat_long <- tot_table_plots[is.na(tot_table_plots$latitude),]

#get latitude longitude values from coverdata
lat_long_table <- unique(coverdata[c("sitemonthyear","plotID","decimalLatitude", "decimalLongitude")]) %>% 
  rename(latitude = decimalLatitude, longitude = decimalLongitude)

#fill latitude longitude values in tot_table_plots where NA
for (i in 1:nrow(missing_lat_long)){
  row <- missing_lat_long[i,]
  lat_long_row <- lat_long_table[lat_long_table$sitemonthyear == row$sitemonthyear & 
                          lat_long_table$plotID == row$plotID,]
  tot_table_plots[tot_table_plots$sitemonthyear == row$sitemonthyear & 
                    tot_table_plots$plotID == row$plotID,]$latitude <- lat_long_row$latitude
  tot_table_plots[tot_table_plots$sitemonthyear == row$sitemonthyear & 
                    tot_table_plots$plotID == row$plotID,]$longitude <- lat_long_row$longitude
}

#creating variables to line up with other datasets
tot_table_plots$siteID <- substr(tot_table_plots$sitemonthyear,1,4) 
tot_table_plots$monthyear <- substr(tot_table_plots$sitemonthyear,5,11) 
tot_table_plots$year <- substr(tot_table_plots$monthyear, 1, 4)

#join veg table data with tot_table_plots
tot_table_plots <- tot_table_plots %>%
  left_join(veg_types)

#re-arranging the columns
tot_table_plots <- tot_table_plots[, c(1, 10, 11, 12 , 2:9, 13)]

#export table of plant cover by plot 
write.table(tot_table_plots, file = "cover_by_plot.csv", sep = ",", row.names = FALSE)

structural_metrics <- left_join(tot_table_plots, lidar_data, by=c("sitemonthyear", "siteID", "monthyear", "easting", "northing", "plotID"))

#export combined table of plant cover and lidar data by plot 
write.table(structural_metrics, file = "structural_metrics_by_plot.csv", sep = ",", row.names = FALSE)
