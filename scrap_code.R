#scrap code



#copy using centroid of tile rather than plot centroid


files2 <- list.files(path="DP3.30006.001/2019", pattern="*.h5", full.names=TRUE, recursive=TRUE)
files3 <- list.files(path="DP3.30006.001/2020", pattern="*.h5", full.names=TRUE, recursive=TRUE)


files2[1]

###
#for each of the 426 bands, I need to calculate the mean reflectance and the SD reflectance across all pixels in the tile
#reflInfo <- h5readAttributes(f, "/SJER/Reflectance/Reflectance_Data")
#myNoDataValue <- as.numeric(reflInfo$Data_Ignore_Value)

myNoDataValue <- -9999

spec_table <- data.frame()
#spec_table <-  data.frame(matrix(ncol = 2, nrow = length(files2))) 
#colnames(spec_table) <- c("tile", "CV") 

#files2 <- c("DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4100000_reflectance.h5", "DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4101000_reflectance.h5")

#files2 <- paste0(wd,"DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4100000_reflectance.h5")
#t <- paste0(wd,"DP3.30006.001/2019/FullSite/D17/2019_SOAP_4/L3/Spectrometer/Reflectance/NEON_D17_SOAP_DP3_296000_4101000_reflectance.h5")

for (t in files3) {
  
  for (i in 1:426){
    #read metadata
    reflInfo <- h5readAttributes(t, "/NIWO/Reflectance/Reflectance_Data")
    
    nRows <- reflInfo$Dimensions[1]
    nCols <- reflInfo$Dimensions[2]
    nBands <- reflInfo$Dimensions[3]
    
    #extract one band
    b <- h5read(t,"/NIWO/Reflectance/Reflectance_Data",index=list(i,1:nCols,1:nRows)) 
    
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
    matrix(c(t, CV),
           ncol = 2)) 
  colnames(out.plot) <- 
    c("tile", "CV") 
  print(out.plot)
  
  #create table that contains 1 row for each tile
  #newspec <- out.plot
  spec_table <- rbind(spec_table, out.plot)
  #spec_table[t] <- out.plot
  
  #spec_table[t, ] <- out.plot
  
  
  h5closeAll()
}

write.table(spec_table, file = "spec_table.csv", sep = ",", row.names = FALSE)


###won't need this 
spec_table <- spec_table %>%
  mutate(siteID = str_sub(tile, start = 38, end = 41))

spec_table$CV <- as.numeric(spec_table$CV)
#once add in easting and northing, can line up with plots to map as a function of SR and cover
ggplot(data = spec_table) +
  geom_point(aes(x = tile, y = CV, color = siteID)) + 
  ylab("Spectral Diversity: CV of Reflectance") + 
  xlab("Tile") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

is.numeric(spec_table$CV)
###










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