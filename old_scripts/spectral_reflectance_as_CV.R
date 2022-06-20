#spectral reflectance code from Chelsea


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
