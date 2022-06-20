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
library(gdata)
library(foreach)
library(tidyverse)
library(raster)
library(rhdf5)
library(rgdal)



wd <- "/Users/rana7082/Documents/research/forest_structural_diversity/data/"
setwd(wd)

getwd()


#####################################
#loop through sites
sites <- list('BONA', 'CLBJ', 'HARV', 'KONZ', 'NIWO', 'ONAQ', 'OSBS', 'SCBI', 'SJER', 'TALL', 'UNDE', 'WREF', 'YELL', 'ABBY', 'MOAB', 'SOAP', 'TEAK', 'HEAL', 'DEJU', 'TREE', 'JERC', 'BART', 'GRSM', 'STEI', 'LENO', 'MLBS', 'UKFS', 'BLAN', 'SERC', 'RMNP', 'DELA')


#prep dates to subset cover data to dates of interest; create list of site/date combos of interest
###
dates <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/NEON_sites_dates_for_cover.csv')

datezz <- dates %>%
  gather(key = "datez", value = "date", -siteID) %>%
  drop_na(.)

#create sitemonthyear as variable to sort the data to the site/date combos of interest
datezz$sitemonthyear <- stringr::str_c(datezz$siteID, datezz$date)

#create a list of dates
datezzz <- as.list(datezz$sitemonthyear)




#creating a dataframe to fill with the plant cover data
tot_cover = data.frame()

#create diversity dataframe of all 31 sites for dates of interest from plot level plant cover data; 
#this contains all dates at each site that will be split out later
for (i in sites) {
  #loading plant cover data product
  cover <- loadByProduct (dpID = "DP1.10058.001", site = i, check.size= FALSE)
  
  #selecting the nested plant diversity information at XX m2 plots
  coverDiv <- cover[[3]]
  
  cover2 <- coverDiv %>%
  
  filter(divDataType=="plantSpecies")
  
  #creating sitemonthyear part 1 to match up with dates of interest
  cover2$monthyear <- substr(cover2$endDate,1,7) 
  
  #part 2
  cover2$sitemonthyear <- stringr::str_c(cover2$siteID, cover2$monthyear)
  
  tot_cover <- rbind(tot_cover, cover2)
}



#filter for the site/date combos of interest
tot_cover <- tot_cover %>%
  filter(sitemonthyear %in% datezzz)

#write .csv for just plant cover diversity data at sites/dates of interest
write.table(tot_cover, file = "prelim_cover.csv", sep = ",", row.names = FALSE)






#############################################
#creating metrics of interest from plant cover data including total species richness (allSR), 
#exotic species richness (exoticSR), exotic plant cover AT SITE LEVEL for site, date combos

#creating a dataframe
tot_table = data.frame()


#looping through dates
for (i in datezzz) {
  
  #filter for site/date combinations of interest
  sub <- tot_cover %>%
    filter (sitemonthyear == i) 
  
  #calculating number of plots within a site
  numplots <-length(unique(sub$plotID))
  
  #total species richness
  all_SR <-length(unique(sub$scientificName))

  #subset of invasive only
  inv <- sub %>%
    filter(nativeStatusCode=="I")

  #total species richness of exotics across all plots
  exotic_SR <-length(unique(inv$scientificName))

  #calculating mean percent cover of exotic species per plot
  exotic_cover <- inv %>%
    summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
    summarize(exotic_cov = mean(sumz))

  i_table <- cbind(numplots, all_SR, exotic_SR, exotic_cover, i)
  
  tot_table <- rbind(tot_table,i_table)
}

#check: 82 site/date combos

#creating variables to line up with other datasets
tot_table$siteID <- substr(tot_table$i,1,4) 
tot_table$monthyear <- substr(tot_table$i,5,11) 
tot_table$year <- substr(tot_table$monthyear, 1, 4)

#check number of sites
count <-unique(tot_table$siteID)
#31



#############################################
#BRING IN NLCD DATA
#read .csv file
veg_types <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/field-sites.csv') %>%
  dplyr::select(Site.ID, Dominant.NLCD.Classes) %>%
  rename(siteID = Site.ID)


#join NLCD information with site level plot cover data
tot_table <- tot_table %>%
  left_join(veg_types)
#############################################




#############################################
#to create allSR, exoticSR, exotic cover AT PLOT LEVEL for site, date combos

#check: this tells us how many combinations of sitemonthyear and plot we should have in our dataframe
test <- unique(tot_cover[c("sitemonthyear","plotID")])
#995

#calculates total species richness for each sitemonthyear, plot combination
all_SR <- tot_cover %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(all_SR = n_distinct(scientificName)) 

#calculates exotic species richness for each sitemonthyear, plot combo
exotic_SR <- tot_cover %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_SR =  n_distinct(scientificName))

#calculates exotic cover for each sitemonthyear, plot combo
exotic_cover <- tot_cover %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_cov = sum(percentCover))

#joins these tables together
closer <- left_join(all_SR, exotic_cover, by=c("sitemonthyear", "plotID")) %>%
  left_join(., exotic_SR, by=c("sitemonthyear", "plotID"))




###
#pulling the plot location information from the plant cover data 

#extracting plot locations
latitude <- unique(tot_cover[c("sitemonthyear","plotID","decimalLatitude")])
longitude <- unique(tot_cover[c("sitemonthyear","plotID","decimalLongitude")])

#adding just longitude to the plant cover data by plot (assuming this will be unique enough to match up with easting)
tot_table_plots <- left_join(closer, latitude, by=c("sitemonthyear", "plotID")) %>%
  left_join(., longitude, by=c("sitemonthyear", "plotID")) 


#creating variables to join by later
tot_table_plots$siteID <- substr(tot_table_plots$sitemonthyear,1,4) 
tot_table_plots$monthyear <- substr(tot_table_plots$sitemonthyear,5,11) 
tot_table_plots$year <- substr(tot_table_plots$monthyear, 1, 4)


######################
#data brought in above
#joining with forest type info from NLCD data
tot_table_plots <- tot_table_plots %>%
  left_join(veg_types)
#######################


#export table of plant cover by plot with longitude
write.table(tot_table_plots, file = "cover_by_plot.csv", sep = ",", row.names = FALSE)



#bring in NEON metadata on plot locations which has both latitude/longitude and easting/northing
plot_centroids <- read.delim('All_NEON_TOS_Plot_Centroids_V8.csv', sep=',', header=T) %>%
  rename(
    decimalLatitude = latitude,
    decimalLongitude = longitude
  )

#left join to add the easting and northing variables from NEON metadata to our plant cover dataset
tot_table_plots_en <- tot_table_plots %>%
  left_join(plot_centroids)

#writing the plant cover plot level data with easting northing that we need to match with lidar data
write.table(tot_table_plots_en, file = "tot_table_plots_en.csv", sep = ",", row.names = FALSE)


#setting NAs for zeros 
tot_table_plots_en$exotic_cov[is.na(tot_table_plots_en$exotic_cov)] <- 0
tot_table_plots_en$exotic_SR[is.na(tot_table_plots_en$exotic_SR)] <- 0

#creating dataset of plot locations (easting/northing only) for each plot that we have plant cover data
plots <- as.data.frame(unique(tot_table_plots_en[c("easting","northing")]))
#537
#############################################









#############################################
#calculate structural metrics

str_table = data.frame()

files <- list.files(path="DP1.30003.001", pattern="*.laz", full.names=TRUE, recursive=TRUE)

tot_table_plots_en <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/tot_table_plots_en.csv')
plots <- as.data.frame(unique(tot_table_plots_en[c("easting","northing")]))

files 


for (q in files)  {
  i <- readLAS(q) # load each file
  # apply function to remove noise
  i <- lidR::classify_noise(i, sor(15,3))
  i <- filter_poi(i, Classification != LASNOISE)
  
  #find extent of tile
  minx <- min(i$X) 
  maxx <- max(i$X)
  miny <- min(i$Y)
  maxy <- max(i$Y)
  
  #select only the plot centroids that are found in this tile (plus a buffer?)
  matches<-filter(plots, plots$easting <= (maxx) & plots$easting >= (minx) & plots$northing <= (maxy) & plots$northing >= (miny)) 
  
  #loop through the matches
  foreach(x = matches$easting, y = matches$northing, .packages="sp") %do% {
 
  #subset to 200 m x 200 m (40,000m2) subtile
      data.200m <- clip_rectangle(i, 
                              xleft = (x - 100), ybottom = (y - 100),
                              xright = (x + 100), ytop = (y + 100))

  
  #creates rasterized digital terrain model
      dtm <- grid_terrain(data.200m, 1, kriging(k = 10L))
  
  #normalize the data based on the digital terrain model (correct for elevation)
      data.200m <- normalize_height(data.200m, dtm)
  
  #replaces really short veg with zeros
      data.200m@data$Z[data.200m@data$Z <= .5] <- 0  

  #Calculate 13 structural metrics in a single function 
  structural_diversity_metrics <- function(data.200m) {
    chm <- grid_canopy(data.200m, res = 1, dsmtin()) 
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
    vert.sd <- cloud_metrics(data.200m, sd(Z, na.rm = TRUE)) 
    sd.1m2 <- grid_metrics(data.200m, sd(Z), 1) 
    sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
    Zs <- data.200m@data$Z
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
  
  #create table that contains 1 row for each plot centroid
  new_structural_diversity <- structural_diversity_metrics(data.200m)
  str_table <- rbind(str_table, new_structural_diversity)
  }  
}


str_table

write.table(str_table, file = "str_table.csv", sep = ",", row.names = FALSE)


######note: 
#will need to match these up both based on easting, northing and by monthyear

str_table <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/str_table.csv')


#join cover table with structure table by easting and northing
#NOTE: this is not joining properly, more than the 995 original observations....
tot_table_plots_en_str <- tot_table_plots_en %>%
  left_join(str_table)
write.table(tot_table_plots_en_str, file = "tot_table_plots_en_str.csv", sep = ",", row.names = FALSE)

################################




