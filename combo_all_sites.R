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


###############################
###if decide to do recent only
#recent <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/NEON_sites_recent_dates.csv')
#recent$sitemonthyear <- stringr::str_c(recent$siteID, recent$recent)

#recentdatezzz <- as.list(recent$sitemonthyear)
###############################

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







#creating metrics of interest from plant cover data including total species richness (allSR), 
#exotic species richness (exoticSR), exotic plant cover at site level for site, date combos

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







###
#looking at summary statistics of percent cover and diversity information based on forest types AT SITE LEVEL

#bring in NLCD information for each site

#read .csv file
veg_types <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/field-sites.csv') %>%
  dplyr::select(Site.ID, Dominant.NLCD.Classes) %>%
  rename(siteID = Site.ID)

#join NLCD information with plot cover data
tot_table <- tot_table %>%
  left_join(veg_types)

#look at distribution of total species richness and exotic species richness across sites 
hist(tot_table$all_SR, breaks = 20, xlab="Total species richness", ylab= "Number of sites", main="")
hist(tot_table$exotic_SR, breaks = 20, xlab="Non-native species richness", ylab= "Number of sites", main="")

#calculate mean exotic species richness across sites
meanz <- tot_table %>%
  filter (exotic_SR > 0) %>%
  summarize(meaninv = mean(exotic_SR))
#8.74

#calculate maximum exotic species richness across sites
maxinv <- tot_table %>%
  filter (exotic_SR > 0) %>%
  summarize(maxinv = max(exotic_SR))
#53; SCBI

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

#check to see how many NAs
test3 <- exotic_cover %>%
  filter(is.na(exotic_cov))
#49 already are NAs

#joins these tables together
closer <- left_join(all_SR, exotic_cover, by=c("sitemonthyear", "plotID")) %>%
  left_join(., exotic_SR, by=c("sitemonthyear", "plotID")) 

#check on how many plots do not have any exotic species
test2 <- tot_table_plots %>%
  filter(is.na(exotic_SR))
#2661; this is correct (2661 + 2102 = 4763)

#check on how many plots do not have any exotic plant cover
test3 <- tot_table_plots %>%
  filter(is.na(exotic_cov))
#2710; this is from 2661 in the 


###
#pulling the plot location information from the plant cover data 
#note: (this is in decimaldegrees and will need to be converted to easting/northing to join with lidar data)

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

#joining with forest type info from NLCD data
tot_table_plots <- tot_table_plots %>%
  left_join(veg_types)

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


#myDataframe[is.na(myDataframe)] = 0
#tot_table_plots_en[is.na(tot_table_plots_en)] = 0

#setting NAs for zeros 
tot_table_plots_en$exotic_cov[is.na(tot_table_plots_en$exotic_cov)] <- 0
tot_table_plots_en$exotic_SR[is.na(tot_table_plots_en$exotic_SR)] <- 0

#rm(plots)
#creating dataset of plot locations (easting/northing only) for each plot that we have plant cover data
plots <- as.data.frame(unique(tot_table_plots_en[c("easting","northing")]))
#537

plots

head(plots)
tail(plots)

#easting <- plots[1]
#northing <- plots[2]

#names(plots)[1] <- "x"
#names(plots)[2] <- "y"

#looking at distribution of total species richness, exotic species richness, and exotic plant cover by plot
hist(tot_table_plots_en$all_SR, breaks = 20, xlab="Total species richness", ylab= "Number of plots", main="")
hist(tot_table_plots_en$exotic_SR, breaks = 20, xlab="Non-native species richness", ylab= "Number of plots", main="")
hist(tot_table_plots_en$exotic_cov, breaks = 40, xlab="Non-native species percent cover", ylab= "Number of plots", main="")

#take out zeros and show just those invaded
onlyinv <- tot_table_plots_en %>%
  filter(exotic_SR > 0)

#looking at only plots that are invaded
hist(onlyinv$exotic_cov, breaks = 40, xlab="Non-native species percent cover", ylab= "Number of plots", main="")


#how many plots are invaded?

#first look at total number of invaded sites/dates
numinv <- tot_table_plots_en %>%
  filter(exotic_SR > 0) %>%
  summarize(plots_w_nonnatives = n_distinct(plotID)) 
#63 site/date combos

#then look at total number of sites/dates
numplots <- tot_table_plots_en %>%
  summarize(totplots = n_distinct(plotID)) 
#82 site/date combos

#calculate the percent of invaded plots across site/date combinations
percent_invaded <- numplots %>%
  left_join(numinv) %>%
  replace_na(list(plots_w_nonnatives = 0)) %>%
  mutate(perc_inv = (plots_w_nonnatives / totplots)*100)


#merge percent_invaded with tot_table

#rm(tot_table_expanded)
tot_table_expanded <- tot_table %>%
  dplyr::rename(sitemonthyear = i) %>%
  left_join(percent_invaded) %>%
  dplyr::select(-numplots)

#creating date field 
tot_table_expanded$date <- lubridate::as_date(tot_table_expanded$monthyear, format = '%Y-%m')
#as_date(x, tz = NULL, format = NULL)


#tot_table_expanded$monthyear<-as.factor(tot_table_expanded$monthyear)
#tot_table_expanded$abis<-strptime(tot_table_expanded$monthyear,format="%Y-%m") #defining what is the original format of your date
#tot_table_expanded$dated<-as.Date(tot_table_expanded$abis,format="%Y-%m")

#don't need this code
recent <- tot_table_expanded %>%
  group_by(siteID) %>%
  slice_max(year) %>%
  filter(sitemonthyear != "SERC2017-07")

is.character(tot_table_expanded$monthyear) #TRUE


tabforpres <- recent %>%
  select(siteID, all_SR, exotic_SR, exotic_cov) %>%
  arrange(all_SR) %>%
  dplyr::rename(nonnative_SR = exotic_SR) %>%
  dplyr::rename(nonnative_cov = exotic_cov) %>%
  dplyr::rename(total_SR = all_SR)


write.table(tabforpres, file = "tabforpres.csv", sep = ",", row.names = FALSE)









#############################################
#calculate structural metrics

str_table = data.frame()

#files <- list.files(path="DP1.30003.001/2019/FullSite/D17/2019_SOAP_4/L1/DiscreteLidar/ClassifiedPointCloud", pattern="*.laz", full.names=TRUE, recursive=FALSE)
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

  #subset to 40 m x 40 m (1600m2) subtile
    #data.40m <- clip_rectangle(data.200m, 
                                  #xleft = (x - 20), ybottom = (y - 20),
                                  #xright = (x + 20), ytop = (y + 20))
  #when subsetting further, get the following error message: Error in { : 
      #task 1 failed - "Interpolation failed (NAs everywhere). Input parameters might be wrong."
  
  
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



sub <- tot_table_plots_en %>%
  filter(siteID == "SOAP" | siteID == "NIWO" | siteID == "DELA")


subby <- sub %>%
  left_join(str_table)

subby[is.na(subby)] = 0

means <- subby %>%
  filter(sitemonthyear == "SOAP2019-06" | sitemonthyear == "NIWO2020-08" | sitemonthyear == "DELA2017-05") %>%
  summarise(meangf = mean(deepgap.fraction.aop), meanoutcanht = mean(mean.max.canopy.ht.aop), meaninthet = mean(sd.sd.aop), meanextht = mean(top.rugosity.aop), meanvertsd = mean(vert.sd.aop), meanentropy = mean(entropy.aop))


#mean outer canopy height
ggplot(data = subby) +
  geom_point(aes(x = mean.max.canopy.ht.aop, y = exotic_cov, color = siteID)) + 
  xlab("Mean Outer Canopy Height (m)") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = mean.max.canopy.ht.aop, y = exotic_SR, color = siteID)) + 
  xlab("Mean Outer Canopy Height (m)") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#gap fraction
ggplot(data = subby) +
  geom_point(aes(x = deepgap.fraction.aop, y = exotic_cov, color = siteID)) +
  xlab("Gap Fraction") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = deepgap.fraction.aop, y = exotic_SR, color = siteID)) +
  xlab("Gap Fraction") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#vert (height) sd
ggplot(data = subby) +
  geom_point(aes(x = vert.sd.aop, y = exotic_cov, color = siteID)) +
  xlab("Height SD") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = vert.sd.aop, y = exotic_SR, color = siteID)) +
  xlab("Height SD") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#sd sd
ggplot(data = subby) +
  geom_point(aes(x = sd.sd.aop, y = exotic_cov, color = siteID)) +
  xlab("SD of Height SD") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = sd.sd.aop, y = exotic_SR, color = siteID)) +
  xlab("SD of Height SD") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#entropy
ggplot(data = subby) +
  geom_point(aes(x = entropy.aop, y = exotic_cov, color = siteID)) +
  xlab("Entropy") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = entropy.aop, y = exotic_SR, color = siteID)) +
  xlab("Entropy") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Top Rugosity, Outer Canopy Roughness
ggplot(data = subby) +
  geom_point(aes(x = top.rugosity.aop, y = exotic_cov, color = siteID)) +
  xlab("Outer Canopy Roughness") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = subby) +
  geom_point(aes(x = top.rugosity.aop, y = exotic_SR, color = siteID)) +
  xlab("Outer Canopy Roughness") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))




###############################################################
#calculate spectral reflectance as CV
#as defined here: https://www.mdpi.com/2072-4292/8/3/214/htm
#f <- paste0(wd,"NEON_D13_NIWO_DP3_453000_4427000_reflectance.h5")


#########################
#based on plot centroid, not on tile


files2 <- list.files(path="DP3.30006.001/2019", pattern="*.h5", full.names=TRUE, recursive=TRUE)
files3 <- list.files(path="DP3.30006.001/2020", pattern="*.h5", full.names=TRUE, recursive=TRUE)
files4 <- list.files(path="DP3.30006.001/2017", pattern="*.h5", full.names=TRUE, recursive=TRUE)

#try with one tile first
#files2 <- paste0(wd,"DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4100000_reflectance.h5")

#t <- paste0(wd,"DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4100000_reflectance.h5")


myNoDataValue <- -9999
spec_table <- data.frame()

for (t in files4) {
  
  reflInfo <- h5readAttributes(t, "/DELA/Reflectance/Reflectance_Data")
  
  #find extent of tile 
  minx <- reflInfo$Spatial_Extent_meters[1] 
  maxx <- reflInfo$Spatial_Extent_meters[2]
  miny <- reflInfo$Spatial_Extent_meters[3]
  maxy <- reflInfo$Spatial_Extent_meters[4]
  
  #select only the plot centroids that are found in this tile (plus a buffer?)
  matches<-filter(plots, plots$easting <= (maxx) & plots$easting >= (minx) & plots$northing <= (maxy) & plots$northing >= (miny)) 
  
  #loop through the matches
  foreach(x = matches$easting, y = matches$northing, .packages="sp") %do% {
    
    for (i in 1:426){
    #read metadata
    
    
    nRows <- reflInfo$Dimensions[1]
    nCols <- reflInfo$Dimensions[2]
    nBands <- reflInfo$Dimensions[3]
    
    #extract one band
    b <- h5read(t,"/DELA/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 
    
    # set all values equal to -9999 to NA
    b[b == myNoDataValue] <- NA
    
    #calculate mean and sd
    meanref <- mean(b, na.rm = TRUE)
    SDref <- sd(b, na.rm = TRUE)
    
    rowz <- cbind(i, meanref, SDref)
    dat <- data.frame()
    dat <- rbind(dat, rowz)
    }
  
  #calculate CV
  dat$calc <- dat$SDref/dat$meanref
  
  CV <- sum(dat$calc)/426
  
  out.plot <- data.frame(
    matrix(c(t, CV, x, y),
           ncol = 4)) 
  colnames(out.plot) <- 
    c("tile", "CV", "easting", "northing") 
  print(out.plot)
  
  #create table that contains 1 row for each plot centroid
  #newspec <- out.plot
  spec_table <- rbind(spec_table, out.plot)
  #spec_table[t] <- out.plot
  
  #spec_table[t, ] <- out.plot
  
  
  h5closeAll()
  }
}
write.table(spec_table, file = "spec_table.csv", sep = ",", row.names = FALSE)

spec_table <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/spec_table.csv')


tot_table_plots_en_str <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/tot_table_plots_en_str.csv')

is.factor(spec_table$easting)

spec_table$easting <- as.numeric(spec_table$easting)
spec_table$northing <- as.numeric(spec_table$northing)

tot_table_plots_en_str_spec <- tot_table_plots_en_str %>%
  left_join(spec_table)


write.table(tot_table_plots_en_str_spec, file = "tot_table_plots_en_str_spec.csv", sep = ",", row.names = FALSE)



sub3 <- tot_table_plots_en_str_spec %>%
  filter(siteID == "SOAP" | siteID == "NIWO" | siteID == "DELA")

#df %>% drop_na(col1)
sub3 <- sub3 %>% drop_na(CV)

#once add in easting and northing, can line up with plots to map as a function of SR and cover
ggplot(data = sub3) +
  geom_point(aes(x = CV, y = exotic_SR, color = siteID)) + 
  xlab("Spectral Diversity: CV of Reflectance") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = sub3) +
  geom_point(aes(x = CV, y = exotic_cov, color = siteID)) + 
  xlab("Spectral Diversity: CV of Reflectance") + 
  ylab("Percent Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#######################################################
#####################################################
#linear regression with structural diversity and cover data

#external heterogeneity
ggplot(sub3, aes(x = top.rugosity.aop, y = log(exotic_cov +1), color = siteID))+
  geom_point()+
  #geom_text(aes(label = SiteID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")+
  xlab("Outer Canopy Roughness") + 
  ylab("Log (% Cover of Non-Native Species) + 1") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ top.rugosity.aop, data = sub3)
summary(fit)
#r2 = 0.07058, p=0.00329


#vert (height) SD
ggplot(sub3, aes(x = vert.sd.aop, y = log(exotic_cov +1), color = siteID))+
  geom_point()+
  #geom_text(aes(label = SiteID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")+
  xlab("SD of Canopy Height") + 
  ylab("Log (% Cover of Non-Native Species) + 1") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#expect to see negative relationship here

fit <- lm(log(exotic_cov +1) ~ vert.sd.aop, data = sub3)
summary(fit)
#r2 = 0.05677, p=0.007714


#entropy
ggplot(sub3, aes(x = entropy.aop, y = log(exotic_cov +1), color = siteID))+
  geom_point()+
  #geom_text(aes(label = SiteID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")+
  xlab("Entropy") + 
  ylab("Log (% Cover of Non-Native Species) + 1") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#expect to see positive relationship here

fit <- lm(log(exotic_cov +1) ~ entropy.aop, data = sub3)
summary(fit)
#r2 = 0.01438, p=0.1136







#CV of reflectance
ggplot(sub3, aes(x = CV, y = log(exotic_cov +1), color = siteID))+
  geom_point()+
  #geom_text(aes(label = SiteID),hjust=0, vjust=0) +
  geom_smooth(method = "lm")+
  xlab("CV of Reflectance") + 
  ylab("Log (% Cover of Non-Native Species) + 1") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#expect to see positive relationship here

fit <- lm(log(exotic_cov +1) ~ CV, data = sub3)
summary(fit)
#r2 = 0.07258, p=0.002906