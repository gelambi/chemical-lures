#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### script 1: BAT COMMUNITIES WITH AND WITHOUT LURES ###

rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(effects)
library(glmmTMB)
library(emmeans)
library(lme4)
library(viridis)
library(viridisLite)
library(MetBrewer)
library(vegan)
library(MetBrewer)
library(ggpubr)
library(parameters)
library(MASS)
library(DHARMa)
library(performance)
library(stats)
library(car)

### DATA

setwd("~/Desktop/lures/objective1") # set directory to objective 1 folder
data <- read.csv("data_nodates.csv")
head(data)
data <- data %>% 
  filter(!bat_species %in% c("bats", "fruit_bats")) # eliminate total numbers

### BAR GRAPHS COMPARING BAT COMMUNITIES CONTROL VERSUS CHEMICAL LURES ###
supp.labs <- c("Control", "Chemical\nlures")
names(supp.labs) <- c("control", "lures")

figure3_total <- ggplot(data, aes(x = site, fill = bat_species, y = values)) +
  theme_test(base_size = 15) +
  ggtitle("c")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  theme(strip.background = element_rect(colour = "NA", fill = "NA")) +
  geom_bar(stat = "identity", color = "NA") +
  scale_fill_viridis_d(option = "inferno", name = "Bat ID",
                       labels = c("artibeus_spp" = parse(text = "italic('Artibeus') ~ 'sp.'"),
                                  "carollia_spp" = parse(text = "italic('Carollia') ~ 'spp.'"),
                                  "desmodus_rotundus" = "Hematophagous\nbats",
                                  "ectophylla.alba" = expression(italic("Ectophylla alba")),
                                  "insectivorous_bats" = "Insectivorous\nbats",
                                  "nectarivorous_bats" = "Nectarivorous\nbats",
                                  "sturnira_spp" = parse(text = "italic('Sturnira') ~ 'spp.'"),
                                  "uroderma_spp" = parse(text = "italic('Uroderma') ~ 'spp.'"))) +
  scale_x_discrete(labels = c("site 1" = "A",
                              "site 2" = "B",
                              "site 3" = "C",
                              "site 4" = "D",
                              "site 5" = "E",
                              "site 6" = "F")) +
  ylab("Total number of bats captured") +
  xlab("Sites") +
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

figure3_total 
ggsave(file="figure3_total.jpg", 
       plot= figure3_total,
       width=18,height=16,units="cm",dpi=300)

### NMDS ###
data <- read.csv("data.csv")
head(data)
mmatrix <- data[, c(8:11, 12)]
head(mmatrix) # just 5 genera of fruit bats
matrix <- mmatrix + 1  #add a small constant to deal with the excess of zeros. Without the constant the NMDS does not run. 
matrix <- as.matrix(matrix) # turn data frame into matrix

nmds_results <- metaMDS(comm = matrix,
                        autotransform = FALSE,
                        distance = "bray",
                        trymax = 100)

nmds_results$stress
plot(nmds_results, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$site <- data$site
data.scores$treatment <- data$treatment
data.scores

#Graph
# PERMANOVA
adonis_site <- adonis2(matrix ~ data.scores$site, distance = "bray", perm=999)
adonis_site
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data$site
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

nmdsgraph_site <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2, color = site)) +
  theme_test(base_size = 15) +
  ggtitle("a")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(size = 2) + 
  scale_color_viridis_d(option="D", name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  annotate("text", x = 1, y = 1, size = 3, label = paste("PERMDISP2, P < 0.01 ***\nPERMANOVA, P = 0.001 ***")) +
  stat_ellipse(level = 0.95, aes(color = site)) + 
  ylim(-1.8, 1.5) +
  xlim(-1.5, 3)
nmdsgraph_site 

# PERMANOVA
adonis_t <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=999)
adonis_t
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)

nmdsgraph_treatment <- ggplot(data.scores, aes(x = MDS1, y = MDS2, color = treatment)) +
  theme_test(base_size = 15) +
  ggtitle("b")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(size = 2) + 
  scale_color_viridis_d(option="D", name= "Treatment", labels = c("control" = "Control", "lures" = "Chemical\nlures")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  annotate("text", x = 1, y = 1, size = 3,label = paste("PERMDISP2, P = 0.553\nPERMANOVA, P = 0.608")) +
  stat_ellipse(level = 0.95) + 
  ylim(-1.8, 1.5) +
  xlim(-1.5, 3)

nmdsgraph_treatment

### Group all community figures together
nmds <- ggarrange(nmdsgraph_site,
                  nmdsgraph_treatment,
                  ncol = 1,
                  nrow = 2,
                  align = "hv")

bat_community <- ggarrange(nmds,
                           figure3_total,
                           ncol = 2,
                           nrow = 1)
bat_community

ggsave(file="bat_community.jpg", 
       plot= bat_community,
       width=11,height=6,units="in",dpi=600)




#######################################
#######################################
#######################################
#######################################
