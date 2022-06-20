#plant cover summary statistics code


###
#SITE LEVEL SUMMARY STATISTICS
#looking at summary statistics of percent cover and diversity information based on forest types AT SITE LEVEL



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
###
#END SITE LEVEL SUMMARY STATISTICS


#####
#check on how many plots do not have any exotic species
test2 <- tot_table_plots %>%
  filter(is.na(exotic_SR))
#2661; this is correct (2661 + 2102 = 4763)

#check on how many plots do not have any exotic plant cover
test3 <- tot_table_plots %>%
  filter(is.na(exotic_cov))
#2710; 





######
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



#tot_table_expanded$monthyear<-as.factor(tot_table_expanded$monthyear)
#tot_table_expanded$abis<-strptime(tot_table_expanded$monthyear,format="%Y-%m") #defining what is the original format of your date
#tot_table_expanded$dated<-as.Date(tot_table_expanded$abis,format="%Y-%m")

#this is to remove a second date for a site if only concerned with the most recent data; won't need this once we have all dates/sites lined up
#recent <- tot_table_expanded %>%
#group_by(siteID) %>%
#slice_max(year) %>%
#filter(sitemonthyear != "SERC2017-07")

#check
is.character(tot_table_expanded$monthyear) #TRUE

#clean table for presentation with diversity metrics of interest
tabforpres <- tot_table_expanded %>%
  select(siteID, all_SR, exotic_SR, exotic_cov) %>%
  arrange(all_SR) %>%
  dplyr::rename(nonnative_SR = exotic_SR) %>%
  dplyr::rename(nonnative_cov = exotic_cov) %>%
  dplyr::rename(total_SR = all_SR)

#writing table
write.table(tabforpres, file = "tabforpres.csv", sep = ",", row.names = FALSE)


