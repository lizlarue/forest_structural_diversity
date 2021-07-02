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

datezz$sitemonthyear <- stringr::str_c(datezz$siteID, datezz$date)

datezzz <- as.list(datezz$sitemonthyear)


###############################
###if decide to do recent only
#recent <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/NEON_sites_recent_dates.csv')
#recent$sitemonthyear <- stringr::str_c(recent$siteID, recent$recent)

#recentdatezzz <- as.list(recent$sitemonthyear)
###############################

tot_cover = data.frame()

#create diversity dataframe of all 31 sites for dates of interest
for (i in sites) {
  cover <- loadByProduct (dpID = "DP1.10058.001", site = i, check.size= FALSE)
  coverDiv <- cover[[3]]
  
  cover2 <- coverDiv %>%
  filter(divDataType=="plantSpecies")
  
  cover2$monthyear <- substr(cover2$endDate,1,7) 

  cover2$sitemonthyear <- stringr::str_c(cover2$siteID, cover2$monthyear)
  
  #cover2[cover2$sitemonthyear %in% datezzz ,]
  
  tot_cover <- rbind(tot_cover, cover2)
}


#tot_cover[tot_cover$sitemonthyear %in% datezzz ,]

tot_cover <- tot_cover %>%
  filter(sitemonthyear %in% datezzz)

write.table(tot_cover, file = "prelim_cover.csv", sep = ",", row.names = FALSE)







#to create allSR, exoticSR, exotic cover at site level for site, date combos
tot_table = data.frame()

for (i in datezzz) {
  
  sub <- tot_cover %>%
    filter (sitemonthyear == i) 
  
  numplots <-length(unique(sub$plotID))
  
  all_SR <-length(unique(sub$scientificName))

  #subset of invasive only
  inv <- sub %>%
    filter(nativeStatusCode=="I")

  #total SR of exotics across all plots
  exotic_SR <-length(unique(inv$scientificName))

  #mean percent cover of exotics
  exotic_cover <- inv %>%
    summarize(sumz = sum(percentCover, na.rm = TRUE)) %>%
    summarize(exotic_cov = mean(sumz))

  i_table <- cbind(numplots, all_SR, exotic_SR, exotic_cover, i)
  
  tot_table <- rbind(tot_table,i_table)
}

#82 site, date combos

tot_table$siteID <- substr(tot_table$i,1,4) 
tot_table$monthyear <- substr(tot_table$i,5,11) 
tot_table$year <- substr(tot_table$monthyear, 1, 4)

count <-unique(tot_table$siteID)
#31



#remove rows where numplots = 0
#tot_table_align <- tot_table %>%
  #filter(numplots != 0)
#57 site date combos that align of where they sampled 

#othercount <- unique(tot_table_align$siteID)
#27; only 27 sites with 0 numplots

#choose the most recent of these???
#one_date <- tot_table %>%
  #group_by(siteID) %>%
  #top_n(n=1)



###
veg_types <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/field-sites.csv') %>%
  dplyr::select(Site.ID, Dominant.NLCD.Classes) %>%
  rename(siteID = Site.ID)

tot_table <- tot_table %>%
  left_join(veg_types)


hist(tot_table$all_SR, breaks = 20, xlab="Total species richness", ylab= "Number of sites", main="")
hist(tot_table$exotic_SR, breaks = 20, xlab="Non-native species richness", ylab= "Number of sites", main="")

meanz <- tot_table %>%
  filter (exotic_SR > 0) %>%
  summarize(meaninv = mean(exotic_SR))
#8.74

maxinv <- tot_table %>%
  filter (exotic_SR > 0) %>%
  summarize(maxinv = max(exotic_SR))
#53; SCBI

#############################################
#to create allSR, exoticSR, exotic cover at plot level for site, date combos

#this tells us how many combinations of sitemonthyear and plot we should have in our dataframe
test <- unique(tot_cover[c("sitemonthyear","plotID")])
#995

#calculates total SR for each sitemonthyear, plot combo
all_SR <- tot_cover %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(all_SR = n_distinct(scientificName)) 

#calculates invasive SR for each sitemonthyear, plot combo
exotic_SR <- tot_cover %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_SR =  n_distinct(scientificName))

#calculates invasive cover for each sitemonthyear, plot combo
exotic_cover <- tot_cover %>%
  filter(nativeStatusCode =="I") %>%
  group_by(sitemonthyear, plotID) %>%
  summarize(exotic_cov = sum(percentCover))

test3 <- exotic_cover %>%
  filter(is.na(exotic_cov))
#49 already are NAs

#joins these tables together
#close <- left_join(all_SR, exotic_SR, by=c("sitemonthyear", "plotID"))
closer <- left_join(all_SR, exotic_cover, by=c("sitemonthyear", "plotID")) %>%
  left_join(., exotic_SR, by=c("sitemonthyear", "plotID")) 

test2 <- tot_table_plots %>%
  filter(is.na(exotic_SR))
#2661; this is correct (2661 + 2102 = 4763)

test3 <- tot_table_plots %>%
  filter(is.na(exotic_cov))
#2710; this is from 2661 in the 



#plot locations
latitude <- unique(tot_cover[c("sitemonthyear","plotID","decimalLatitude")])
longitude <- unique(tot_cover[c("sitemonthyear","plotID","decimalLongitude")])

tot_table_plots <- left_join(closer, latitude, by=c("sitemonthyear", "plotID")) %>%
  left_join(., longitude, by=c("sitemonthyear", "plotID")) 



tot_table_plots$siteID <- substr(tot_table_plots$sitemonthyear,1,4) 
tot_table_plots$monthyear <- substr(tot_table_plots$sitemonthyear,5,11) 
tot_table_plots$year <- substr(tot_table_plots$monthyear, 1, 4)

tot_table_plots <- tot_table_plots %>%
  left_join(veg_types)

write.table(tot_table_plots, file = "cover_by_plot.csv", sep = ",", row.names = FALSE)



#bring in NEON metadata on plot locations
plot_centroids <- read.delim('All_NEON_TOS_Plot_Centroids_V8.csv', sep=',', header=T) %>%
  rename(
    decimalLatitude = latitude,
    decimalLongitude = longitude
  )

#left join to add the easting and northing variables
tot_table_plots_en <- tot_table_plots %>%
  left_join(plot_centroids)

write.table(tot_table_plots_en, file = "tot_table_plots_en.csv", sep = ",", row.names = FALSE)


#myDataframe[is.na(myDataframe)] = 0
tot_table_plots_en[is.na(tot_table_plots_en)] = 0

rm(plots)
plots <- as.data.frame(unique(tot_table_plots_en[c("easting","northing")]))
#766

plots

head(plots)
tail(plots)

#easting <- plots[1]
#northing <- plots[2]

#names(plots)[1] <- "x"
#names(plots)[2] <- "y"


hist(tot_table_plots_en$all_SR, breaks = 20, xlab="Total species richness", ylab= "Number of plots", main="")
hist(tot_table_plots_en$exotic_SR, breaks = 20, xlab="Non-native species richness", ylab= "Number of plots", main="")
hist(tot_table_plots_en$exotic_cov, breaks = 40, xlab="Non-native species percent cover", ylab= "Number of plots", main="")

numinv <- tot_table_plots_en %>%
  filter(exotic_SR > 0) %>%
  summarize(plots_w_nonnatives = n_distinct(plotID)) 
#63 site/date combos

numplots <- tot_table_plots_en %>%
  summarize(totplots = n_distinct(plotID)) 
#82 site/date combos

percent_invaded <- numplots %>%
  left_join(numinv) %>%
  replace_na(list(plots_w_nonnatives = 0)) %>%
  mutate(perc_inv = (plots_w_nonnatives / totplots)*100)


#merge percent_invaded with tot_table


tot_table_expanded <- tot_table %>%
  rename(sitemonthyear = i) %>%
  left_join(percent_invaded) %>%
  dplyr::select(-numplots)


#not sure what this tells us?
meaninv <- tot_table_plots_en %>%
  summarize(meaninv = mean(exotic_cov))

check <- tot_table_plots_en %>%
  summarize(n_distinct(sitemonthyear))
#82






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


######note: 
#will need to match these up both based on easting, northing and by monthyear

#join cover table with structure table by easting and northing
tot_table_plots_en_str <- tot_table_plots_en %>%
  left_join(str_table)
write.table(tot_table_plots_en_str, file = "tot_table_plots_en_str.csv", sep = ",", row.names = FALSE)

################################

sub <- tot_table_plots_en %>%
  filter(siteID == "SOAP" | siteID == "NIWO")


subby <- sub %>%
  left_join(str_table)


means <- subby %>%
  filter(sitemonthyear == "SOAP2019-06") %>%
  summarise(meangf = mean(deepgap.fraction.aop), meanoutcanht = mean(mean.max.canopy.ht.aop), meaninthet = mean(sd.sd.aop), meanextht = mean(top.rugosity.aop))

ggplot(data = subby) +
  geom_point(aes(x = deepgap.fraction.aop, y = exotic_cov, color = siteID)) +
  xlab("Gap Fraction") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(data = subby) +
  geom_point(aes(x = mean.max.canopy.ht.aop, y = exotic_cov, color = siteID)) + 
  xlab("Mean Outer Canopy Height (m)") + 
  ylab("% Cover of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(data = subby) +
  geom_point(aes(x = deepgap.fraction.aop, y = exotic_SR, color = siteID)) +
  xlab("Gap Fraction") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(data = subby) +
  geom_point(aes(x = mean.max.canopy.ht.aop, y = exotic_SR, color = siteID)) + 
  xlab("Mean Outer Canopy Height (m)") + 
  ylab("Species Richness of Non-Native Species") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))



###############################################################
#copy using centroid of tile rather than plot centroid

#############################################
#calculate structural metrics

str_table = data.frame()

#files <- list.files(path="DP1.30003.001/2019/FullSite/D17/2019_SOAP_4/L1/DiscreteLidar/ClassifiedPointCloud", pattern="*.laz", full.names=TRUE, recursive=FALSE)
files <- list.files(path="DP1.30003.001", pattern="*.laz", full.names=TRUE, recursive=TRUE)


files 


for (q in files)  {
  i <- readLAS(q) # load each file
  # apply function to remove noise
  i <- lidR::classify_noise(i, sor(15,3))
  i <- filter_poi(i, Classification != LASNOISE)
  
  #find center point of tile
  x <- ((max(i$X) - min(i$X))/2)+ min(i$X)
  y <- ((max(i$Y) - min(i$Y))/2)+ min(i$Y)
  
  
    #subset to 40 m x 40 m (1600m2) subtile
    data.40m <- clip_rectangle(i, 
                               xleft = (x - 20), ybottom = (y - 20),
                               xright = (x + 20), ytop = (y + 20))
    
    
    #creates rasterized digital terrain model
    dtm <- grid_terrain(data.40m, 1, kriging(k = 10L))
    
    #normalize the data based on the digital terrain model (correct for elevation)
    data.40m <- normalize_height(data.40m, dtm)
    
    #might not need this step
    #subset to a 40 m x 40 m (1600m2) subtile
    #data.40m <- clip_rectangle(data.200m, 
    #xleft = (x - 20), ybottom = (y - 20),
    #xright = (x + 20), ytop = (y + 20))
    
    
    #replaces really short veg with zeros
    data.40m@data$Z[data.40m@data$Z <= .5] <- 0  
    #plot(data.40m)
    
    
    #Calculate 13 structural metrics in a single function 
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
    
    #create table that contains 1 row for each tile
    new_structural_diversity <- structural_diversity_metrics(data.40m)
    str_table <- rbind(str_table, new_structural_diversity)
}


str_table



#join cover table with structure table by easting and northing
tot_table_plots_en_str <- tot_table_plots_en %>%
  left_join(str_table)
#no points in common












#############################################
#calculate spectral reflectance as CV
#as defined here: https://www.mdpi.com/2072-4292/8/3/214/htm
#f <- paste0(wd,"NEON_D13_NIWO_DP3_453000_4427000_reflectance.h5")

files2 <- list.files(path="DP3.30006.001/2019", pattern="*.h5", full.names=TRUE, recursive=TRUE)
files3 <- list.files(path="DP3.30006.001/2020", pattern="*.h5", full.names=TRUE, recursive=TRUE)


files2[1]

###
#for each of the 426 bands, I need to calculate the mean reflectance and the SD reflectance across all pixels in the tile
#reflInfo <- h5readAttributes(f, "/SJER/Reflectance/Reflectance_Data")
#myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)

myNoDataValue <- -9999

dat <- data.frame()
spec_table <- data.frame()


for (t in files2) {

  for (i in 1:426){
  #read metadata
    reflInfo <- h5readAttributes(t, "/SOAP/Reflectance/Reflectance_Data")
    
    nRows <- reflInfo$Dimensions[1]
    nCols <- reflInfo$Dimensions[2]
    nBands <- reflInfo$Dimensions[3]
    
  #extract one band
    b <- h5read(t,"/SOAP/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 
  
  # set all values equal to -9999 to NA
    b[b == myNoDataValue] <- NA
  
  #calculate mean and sd
    meanref <- mean(b, na.rm = TRUE)
    SDref <- sd(b, na.rm = TRUE)
  
    rowz <- cbind(i, meanref, SDref)
  
    dat <- rbind(dat, rowz)
  }

dat$calc <- dat$SDref/dat$meanref

CV <- sum(dat$calc)/426


newspec <- CV
spec_table <- rbind(spec_table, newspec)

}

h5closeAll()
