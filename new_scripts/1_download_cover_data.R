# This script will download the cover data for all the date-site 
# combinations and will put all of it in a csv.

library(lidR)
library(gstat)
library(neonUtilities)

root <- paste0(getwd(), "/")

data <- read.csv(paste0(root, "NEON_sites_dates_for_cover.csv"))
#print(data)

tot_cover = data.frame()

for (i in 1:nrow(data)) {
  site <- data[i,][1]
  print(site)
  for (j in 2:ncol(data))
  {
    cell <- data[i,][j]
    if (!is.na(cell)){
      cell <- toString(cell[[1]])
      print(cell)
      
      cover <- loadByProduct (dpID = "DP1.10058.001", site = site, 
                              startdate = cell, enddate = cell,
                              check.size= F)#, savepath = "~/plantcover/")
      coverDiv <- cover$div_1m2Data
      coverDiv2 <- coverDiv[coverDiv$divDataType == 'plantSpecies',]
      coverDiv2$monthyear <- substr(coverDiv2$endDate,1,7) 
      coverDiv2$sitemonthyear <- stringr::str_c(coverDiv2$siteID, coverDiv2$monthyear)
      tot_cover <- rbind(tot_cover, coverDiv2)
    }
  }
}

write.table(tot_cover, file = "prelim_cover.csv", sep = ",", row.names = FALSE)
