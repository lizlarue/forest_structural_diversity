#plotting structural metrics




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

