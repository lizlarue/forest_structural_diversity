#FIA data

#Chelsea Nagy, University of Colorado Boulder
#October 30, 2020

invsub <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_INVASIVE_SUBPLOT_SPP.csv') 
vegstr <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_P2VEG_SUBP_STRUCTURE.csv') 
vegsub <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_P2VEG_SUBPLOT_SPP.csv') 
soils <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_SOILS_LAB.csv') 
treebio <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_TREE_REGIONAL_BIOMASS.csv') 
species <- read.csv(file = '/Users/rana7082/Documents/research/forest_structural_diversity/data/CO_VEG_PLOT_SPECIES.csv') 

#this gives invasive species richness by plot
invSRplot <- invsub %>%
  group_by(PLOT) %>%
  summarise(SR = length(unique(VEG_FLDSPCD)))
