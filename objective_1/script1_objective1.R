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

#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### Bats captured using mist nets 
setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_1") # set directory to objective 1 folder
data <- read.csv("data_nodates.csv")
head(data)
data <- data  %>% filter(!bat_species %in% c("bats", "fruit_bats"))
### Bar graphs comparing bat communities control versus chemical lures 
supp.labs <- c("Control", "Chemical lures")
names(supp.labs) <- c("control", "lures")

figure_1 <-ggplot(data, aes(x = site, fill = bat_species, y = values)) +
  theme_pubclean(base_size = 15) +
  theme(panel.border = element_blank(), strip.background = element_rect(colour="NA", fill="NA"))+
  geom_bar(stat = "identity", color = "NA") +
  scale_fill_viridis_d(option = "inferno", name= "Bat ID", labels = c("artibeus_spp" = "Artibeus sp.",
                                                      "carollia_spp" = "Carollia spp.",
                                                      "desmodus_rotundus" = "Hematophagous bats",
                                                      "ectophylla.alba" = "Ectophylla alba",
                                                      "insectivorous_bats" = "Insectivorous bats",
                                                      "nectarivorous_bats" = "Nectarivorous bats",
                                                      "sturnira_spp" = "Sturnira spp.",
                                                      "uroderma_spp" = "Uroderma spp.")) + 
  scale_x_discrete(labels = c("site 1" = "A",
                              "site 2" = "B",
                              "site 3" = "C",
                              "site 4" = "D",
                              "site 5" = "E",
                              "site 6" = "F")) + 
  ylab ("Total number of bats captured per site") +
  xlab (" ") + 
  theme(legend.position = "right") + 
  facet_wrap(~treatment,
             labeller = labeller(treatment = supp.labs))

figure_1
ggsave(file="figure_1.jpg", 
       plot= figure_1,
       width=16,height=14,units="cm",dpi=300)

### all bats together in one GLMM
fruitbatsp_data <- data  %>% filter(!bat_species %in% c("bats", 
                                             "fruit_bats",
                                             "nectarivorous_bats",
                                             "insectivorous_bats",
                                             "desmodus_rotundus" )) # delete categories that are not fruit bat genera/species

fruitbatsp_data

fruitbatsp_glmm <- glmmTMB(values ~ treatment*bat_species + (1|site), data = fruitbatsp_data, family = poisson)
Anova(fruitbatsp_glmm)

summary(fruitbatsp_glmm)

# models 
fruitbatsp_glmm <- glmmTMB(values ~ treatment*bat_species + (1|site), data = fruitbatsp_data, family = poisson)
fruitbatsp_glmm_nb1 <- glmmTMB(values ~ treatment*bat_species + (1|site), data = fruitbatsp_data, family = nbinom1(link = "log"))
fruitbatsp_glmm_nb2 <- glmmTMB(values ~ treatment*bat_species + (1|site), data = fruitbatsp_data, family = nbinom2(link = "log"))
fruitbatsp_glmm_zi <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family = poisson)
fruitbatsp_glmm_nb1_zi <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family = nbinom1(link = "log"))
fruitbatsp_glmm_nb2_zi <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family = nbinom2(link = "log"))
fruitbatsp_glmm_hurdle_poisson_zi <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family = truncated_poisson)
fruitbatsp_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family=truncated_nbinom2)
fruitbatsp_glmm_zi_site <- glmmTMB(values ~ treatment*bat_species + (1|site), zi = ~1, data = fruitbatsp_data, family = poisson)

# models with zi have convergence problem; extreme or very small eigenvalues detected.
# compare all models
check_distribution(fruitbatsp_glmm) #  neg. binomial (zero-infl.) 62%

AICtab(fruitbatsp_glmm, # poisson
       fruitbatsp_glmm_nb1, fruitbatsp_glmm_nb2, # negative binomial
       fruitbatsp_glmm_zi, # zero inflated poisson
       fruitbatsp_glmm_nb1_zi, fruitbatsp_glmm_nb2_zi,  # zero inflated negative binomial
       fruitbatsp_glmm_hurdle_poisson_zi, # hurdle model, poisson
       fruitbatsp_glmm_hurdle_bn2_zi, # hurdle model, negative binomial
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(fruitbatsp_glmm_nb2) ### tab best model: # negative binomial 1
Anova(fruitbatsp_glmm_nb2) # overall effects

# Check for overdispesion and zero inflation 
summary(fruitbatsp_glmm_nb2)
check_overdispersion(fruitbatsp_glmm_nb2) # No overdispersion detected.
check_zeroinflation(fruitbatsp_glmm_nb2) # Model is overfitting zeros.
# Check for overdispesion and zero inflation 
summary(fruitbatsp_glmm_nb2_zi)
check_overdispersion(fruitbatsp_glmm_nb2_zi) # No overdispersion detected.
check_zeroinflation(fruitbatsp_glmm_nb2_zi) # Model is overfitting zeros.
# Check for overdispesion and zero inflation 
summary(fruitbatsp_glmm_nb1)
check_overdispersion(fruitbatsp_glmm_nb1) # No overdispersion detected.
check_zeroinflation(fruitbatsp_glmm_nb1) # Model is underfitting zeros (probable zero-inflation).
# Check for overdispesion and zero inflation 
summary(fruitbatsp_glmm)
check_overdispersion(fruitbatsp_glmm) # Overdispersion detected.
check_zeroinflation(fruitbatsp_glmm_nb1) # Model is underfitting zeros (probable zero-inflation).

# calculate effect size
glmm_emmeans <-emmeans(fruitbatsp_glmm_nb2,~treatment*bat_species, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
fruitbatsp_data$prediction_fruitbatsp <- predict(fruitbatsp_glmm_nb1, fruitbatsp_data, type="response")
glmm_graph_allbats <- ggplot(fruitbatsp_data, aes(x = treatment, y = prediction_fruitbatsp)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = fruitbatsp_data, aes(x = treatment, y = values, group = bat_species, position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = fruitbatsp_data, aes(x = treatment, y = values, color = bat_species, shape = site), width = 0.1, size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="inferno", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of fruit bat species") +
  xlab ("") +
  facet_wrap(~bat_species)
  
glmm_graph_allbats
ggsave(file="fruitbatsp_a.jpg", 
       plot= glmm_graph_allbats,
       width=15,height=15,units="cm",dpi=300)

# The differences in species abundances might be causing biased in the model!!!!

### Standardize counts
fruitbatsp_data <- fruitbatsp_data %>%
  group_by(bat_species) %>%
  mutate(values_std = (values - mean(values))/sd(values))

# models 
fruitbatsp_glmm <- glmmTMB(values_std ~ treatment*bat_species + (1|site), data = fruitbatsp_data)
check_distribution(fruitbatsp_glmm) #  neg. binomial (zero-infl.) 62%
tab_model(fruitbatsp_glmm) ### tab best model: # negative binomial 1
Anova(fruitbatsp_glmm) # overall effects
check_heteroscedasticity(fruitbatsp_glmm)

# calculate effect size
glmm_emmeans <-emmeans(fruitbatsp_glmm,~treatment*bat_species, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
supp.labs <- c("Artibeus sp.", "Carollia spp.", "Ectophylla alba", "Sturnira spp.", "Uroderma spp.")
names(supp.labs) <- c("artibeus_spp", "carollia_spp","ectophylla.alba","sturnira_spp", "uroderma_spp")

fruitbatsp_data$prediction_fruitbatsp <- predict(fruitbatsp_glmm, fruitbatsp_data, type="response")
glmm_graph_allbats <- ggplot(fruitbatsp_data, aes(x = treatment, y = prediction_fruitbatsp)) +
  theme_pubclean(base_size = 15) +  
  theme(panel.border = element_blank(), strip.background = element_rect(colour="NA", fill="NA")) +
  geom_line(data = fruitbatsp_data, aes(x = treatment, y = values_std, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") + 
  geom_jitter(data = fruitbatsp_data, aes(x = treatment, y = values_std, color = site), position = position_dodge(0.06), size = 3) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Standardized counts of fruit bat species") +
  xlab ("") +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~bat_species, labeller = labeller(bat_species= supp.labs))

glmm_graph_allbats
ggsave(file="fruitbatsp_a.jpg", 
       plot= glmm_graph_allbats,
       width=30,height=20,units="cm",dpi=300)

# all fruit bats (JUST total numbers)
data <- read.csv("data_nodates.csv")
head(data)
fruitbats_data <- data %>% filter(bat_species %in% c("fruit_bats"))

# models 
fruitbats_glmm <- glmmTMB(values ~ treatment + (1|site), data = fruitbats_data, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = fruitbats_data, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = fruitbats_data, family = nbinom2(link = "log"))
fruitbats_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = fruitbats_data, family = poisson)
fruitbats_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = fruitbats_data, family = nbinom1(link = "log"))
fruitbats_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = fruitbats_data, family = nbinom2(link = "log"))
fruitbats_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = fruitbats_data, family = truncated_poisson)
fruitbats_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = fruitbats_data, family=truncated_nbinom2)

# compare all models
check_distribution(fruitbats_glmm) # neg. binomial (zero-infl.) 84%

AICtab(fruitbats_glmm,
       fruitbats_glmm_nb1,
       fruitbats_glmm_nb2,
       fruitbats_glmm_zi,
       fruitbats_glmm_nb1_zi,
       fruitbats_glmm_nb2_zi,
       fruitbats_glmm_hurdle_poisson_zi,
       fruitbats_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(fruitbats_glmm_nb1) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(fruitbats_glmm_nb1)
check_overdispersion(fruitbats_glmm_nb1) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_nb1) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(fruitbats_glmm_nb1))
check_model(fruitbats_glmm_nb1)

# calculate effect size
glmm_emmeans <-emmeans(fruitbats_glmm_nb1,~treatment, type="response")
glmm_emmeans
  # Fold change
  13.4/11 # Sites with lures capture 1.2 more bats 

# compute predictions for graph using the function predict
fruitbats_data$prediction_fruitbats <- predict(fruitbats_glmm_nb1, fruitbats_data, type="response")
glmm_graph_allbats <- ggplot(fruitbats_data, aes(x = treatment, y = prediction_fruitbats)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = fruitbats_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = fruitbats_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of fruit bat species") +
  xlab ("") + 
  theme(legend.position = "bottom")
glmm_graph_allbats
ggsave(file="fruitbats_a.jpg", 
       plot= glmm_graph_allbats,
       width=15,height=15,units="cm",dpi=300)

# carollia
data
carollia_data <- data %>% filter(bat_species %in% c("carollia_spp"))

# models 
carollia_glmm <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = carollia_data, family = poisson)
carollia_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = carollia_data, family = truncated_poisson)
carollia_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~site, data = carollia_data, family=truncated_nbinom2)

# compare all models
check_distribution(carollia_glmm) # neg. binomial (zero-infl.) 56%

AICtab(carollia_glmm,
       carollia_glmm_nb1,
       carollia_glmm_nb2,
       carollia_glmm_zi,
       carollia_glmm_nb1_zi,
       carollia_glmm_nb2_zi,
       carollia_glmm_hurdle_poisson_zi,
       carollia_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(carollia_glmm) ### tab best model

# Check for overdispesion and zero inflation 
summary(carollia_glmm)
check_overdispersion(carollia_glmm) # No overdispersion detected.
check_zeroinflation(carollia_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(carollia_glmm))
check_model(carollia_glmm)

# calculate effect size
glmm_emmeans <-emmeans(carollia_glmm_nb1,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
carollia_data$prediction_carollia <- predict(carollia_glmm_nb1, carollia_data, type="response")
glmm_graph_carollia <- ggplot(carollia_data, aes(x = treatment, y = prediction_carollia)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = carollia_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06),color = "black") +
  geom_jitter(data = carollia_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Carollia spp.") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,14,by=2))
glmm_graph_carollia
ggsave(file="carollia_a.jpg", 
       plot= glmm_graph_carollia,
       width=15,height=15,units="cm",dpi=300)

# uroderma
data
uroderma_data <- data %>% filter(bat_species %in% c("uroderma_spp"))

# models 
uroderma_glmm <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = poisson)
uroderma_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = uroderma_data, family = poisson)
uroderma_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = uroderma_data, family = truncated_poisson)
uroderma_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = uroderma_data, family=truncated_nbinom2)

# compare all models
check_distribution(uroderma_glmm) # neg. binomial (zero-infl.) 56%

AICtab(uroderma_glmm,
       uroderma_glmm_nb1,
       uroderma_glmm_nb2,
       uroderma_glmm_zi,
       uroderma_glmm_nb1_zi,
       uroderma_glmm_nb2_zi,
       uroderma_glmm_hurdle_poisson_zi,
       uroderma_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(uroderma_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(uroderma_glmm)
check_overdispersion(uroderma_glmm) # No overdispersion detected.
check_zeroinflation(uroderma_glmm) # Model is overfitting zeros.
plot(allEffects(uroderma_glmm))
check_model(uroderma_glmm)

# calculate effect size
glmm_emmeans <-emmeans(uroderma_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
uroderma_data$prediction_uroderma <- predict(uroderma_glmm, uroderma_data, type="response")
glmm_graph_uroderma <- ggplot(uroderma_data, aes(x = treatment, y = prediction_uroderma)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = uroderma_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = uroderma_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Uroderma spp.") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,40,by=10))
glmm_graph_uroderma
ggsave(file="uroderma_a.jpg", 
       plot= glmm_graph_uroderma,
       width=15,height=15,units="cm",dpi=300)

# ectophylla alba
data
ectophylla_data <- data %>% filter(bat_species %in% c("ectophylla.alba"))

# models 
ectophylla_glmm <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = poisson)
ectophylla_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = nbinom1(link = "log"))
ectophylla_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = nbinom2(link = "log"))
ectophylla_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = ectophylla_data, family = poisson)
ectophylla_glmm_nbsite_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = ectophylla_data, family = nbinom1(link = "log"))
ectophylla_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = ectophylla_data, family = nbinom2(link = "log"))
ectophylla_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = ectophylla_data, family = truncated_poisson)
ectophylla_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = ectophylla_data, family=truncated_nbinom2)

# compare all models
check_distribution(ectophylla_glmm) # neg. binomial (zero-infl.) 53%

AICtab(ectophylla_glmm,
       ectophylla_glmm_nb1,
       ectophylla_glmm_nb2,
       ectophylla_glmm_zi,
       ectophylla_glmm_nb1_zi,
       ectophylla_glmm_nb2_zi,
       ectophylla_glmm_hurdle_poisson_zi,
       ectophylla_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(ectophylla_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(ectophylla_glmm)
check_overdispersion(ectophylla_glmm) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(ectophylla_glmm))
check_model(ectophylla_glmm)

summary(ectophylla_glmm_nb1)
check_overdispersion(ectophylla_glmm_nb1) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm_nb1) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(ectophylla_glmm_nb1))
check_model(ectophylla_glmm_nb1)

# calculate effect size
glmm_emmeans <-emmeans(ectophylla_glmm_nb1,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
ectophylla_data$prediction_ectophylla <- predict(ectophylla_glmm_nb1, ectophylla_data, type="response")
glmm_graph_ectophylla <- ggplot(ectophylla_data, aes(x = treatment, y = prediction_ectophylla)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = ectophylla_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = ectophylla_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Ectophylla alba") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,75,by=15))
glmm_graph_ectophylla
ggsave(file="ectophylla_a.jpg", 
       plot= glmm_graph_ectophylla,
       width=15,height=15,units="cm",dpi=300)

# artibeus
data
artibeus_data <- data %>% filter(bat_species %in% c("artibeus_spp"))

# models 
artibeus_glmm <- glmmTMB(values ~ treatment + (1|site), data = artibeus_data, family = poisson)
artibeus_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = artibeus_data, family = nbinom1(link = "log"))
artibeus_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = artibeus_data, family = nbinom2(link = "log"))
artibeus_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = artibeus_data, family = poisson)
artibeus_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = artibeus_data, family = nbinom1(link = "log"))
artibeus_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = artibeus_data, family = nbinom2(link = "log"))
artibeus_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = artibeus_data, family = truncated_poisson)
artibeus_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = artibeus_data, family=truncated_nbinom2)

# compare all models
check_distribution(artibeus_glmm) # neg. binomial (zero-infl.) 53%

AICtab(artibeus_glmm,
       artibeus_glmm_nb1,
       artibeus_glmm_nb2,
       artibeus_glmm_zi,
       artibeus_glmm_nb1_zi,
       artibeus_glmm_nb2_zi,
       artibeus_glmm_hurdle_poisson_zi,
       artibeus_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(artibeus_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(artibeus_glmm)
check_overdispersion(artibeus_glmm) # No overdispersion detected.
check_zeroinflation(artibeus_glmm) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(artibeus_glmm))
check_model(artibeus_glmm)

# calculate effect size
glmm_emmeans <-emmeans(artibeus_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
artibeus_data$prediction_artibeus <- predict(artibeus_glmm, artibeus_data, type="response")
glmm_graph_artibeus <- ggplot(artibeus_data, aes(x = treatment, y = prediction_artibeus)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = artibeus_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = artibeus_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Artibeus sp.") +
  xlab ("") +
  scale_y_continuous(breaks=seq(0,1,by=1))

glmm_graph_artibeus
ggsave(file="artibeus_a.jpg", 
       plot= glmm_graph_artibeus,
       width=15,height=15,units="cm",dpi=300)

# sturnira
data
sturnira_data <- data %>% filter(bat_species %in% c("sturnira_spp"))

# models 
sturnira_glmm <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = poisson)
sturnira_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom1(link = "log"))
sturnira_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = sturnira_data, family = poisson)
sturnira_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = sturnira_data, family = nbinom1(link = "log"))
sturnira_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_hurdle_poisson_zi  <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = sturnira_data, family = truncated_poisson)
sturnira_glmm_hurdle_bn2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~1, data = sturnira_data, family=truncated_nbinom2)

# compare all models
check_distribution(sturnira_glmm) # neg. binomial (zero-infl.) 53%

AICtab(sturnira_glmm,
       sturnira_glmm_nb1,
       sturnira_glmm_nb2,
       sturnira_glmm_zi,
       #sturnira_glmm_nb1_zi,
       sturnira_glmm_nb2_zi,
       sturnira_glmm_hurdle_poisson_zi,
       sturnira_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(sturnira_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(sturnira_glmm)
check_overdispersion(sturnira_glmm) # No overdispersion detected.
check_zeroinflation(sturnira_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm))
check_model(sturnira_glmm)

# calculate effect size
glmm_emmeans <-emmeans(sturnira_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
sturnira_data$prediction_sturnira <- predict(sturnira_glmm, sturnira_data, type="response")
glmm_graph_sturnira <- ggplot(sturnira_data, aes(x = treatment, y = prediction_sturnira)) +
  theme_pubclean(base_size = 15) +  
  #geom_line(data = sturnira_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = sturnira_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Sturnira spp.") +
  xlab ("")
glmm_graph_sturnira
ggsave(file="sturnira_a.jpg", 
       plot= glmm_graph_sturnira,
       width=15,height=15,units="cm",dpi=300)

sp_bats <- ggarrange(glmm_graph_carollia,
                     glmm_graph_ectophylla,
                     glmm_graph_uroderma,
                     glmm_graph_sturnira,
                     glmm_graph_artibeus,
                     ncol = 2,
                     nrow = 3,
                     align = "hv",
                     common.legend = TRUE,
                     legend="none")
sp_bats

figure_2 <- ggarrange(glmm_graph_allbats,
                      sp_bats,
                      ncol = 2,
                      nrow = 1)

figure_2
ggsave(file="figure_2_a.jpg", 
       plot= figure_2,
       width=35,height=20,units="cm",dpi=300)
