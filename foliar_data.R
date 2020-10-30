#playing with foliar data from ORNL

library(ggplot2)
library(dplyr)
library(sp)

wd <- "/Users/rana7082/Dropbox/forest_structure_disturbance/"

#data from: https://daac.ornl.gov/VEGETATION/guides/Leaf_carbon_nutrients.html
foliar <- read.csv(paste0(wd,"LEAF_CARBON_NUTRIENTS_1106/data/Leaf_Carbon_Nutrients_data.csv"))

foliar_US <- foliar %>%
  filter(Country == "USA")

unique(foliar_US$Region)

is.factor(foliar_US$Region)

foliar2 <- foliar_US %>% 
  mutate(Region = replace(Region, Region == 'Alaska', 'AK'))

unique(foliar2$Region)

foliar2 <- foliar2 %>%
  mutate(N_P_green_leaf = N_green_leaf/P_green_leaf)

summary(foliar2$N_P_green_leaf)  

foliar2 <- foliar2 %>%
  rename(State = Region)

foliar2 <- foliar2 %>%
  mutate(Region = ifelse(State == "HI"| State == "CO"| State == "AK" | State == "WA" | State == "CA"| State == "NM" | State == "UT" | State == "NV", 'west', 'east'))

foliar2NP <- foliar2[!is.na(foliar2$N_P_green_leaf),]

sumzNP <- foliar2NP %>%
  group_by(Vegetation) %>%
  summarize(mean = mean(N_P_green_leaf))

sumzNPgrowth <- foliar2NP %>%
  group_by(Growth_habit) %>%
  summarize(mean = mean(N_P_green_leaf))
