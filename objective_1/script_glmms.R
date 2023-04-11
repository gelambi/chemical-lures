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
library(bbmle) 
library(sjPlot)

#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### read data: bats captured using mist nets 
setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_1") # set directory to objective 1 folder
data <- read.csv("data.csv")
head(data)

################
### all data ###
################

ha_glmm <- glmmTMB(bats ~ treatment + (1|site), data = data, family = poisson)
ha_glmm_nb1 <- glmmTMB(bats ~ treatment + (1|site), data = data, family = nbinom1(link = "log"))
ha_glmm_nb2 <- glmmTMB(bats ~ treatment + (1|site), data = data, family = nbinom2(link = "log"))
ha_glmm_zi <- glmmTMB(bats ~ treatment + (1|site), data = data, zi =~ treatment, family = poisson)
ha_glmm_zi_nb2 <- glmmTMB(bats ~ treatment + (1|site), data = data, zi =~ treatment, family = nbinom2(link = "log"))
ha_glmm_hurdle <- glmmTMB(bats ~ treatment + (1|site), data = data, zi =~ treatment, family=truncated_poisson)

# poisson
summary(ha_glmm)
glmm1 <- tab_model(ha_glmm)
r2(ha_glmm)
check_overdispersion(ha_glmm) # Overdispersion detected
check_zeroinflation(ha_glmm) # Model is underfitting zeros (probable zero-inflation)
check_distribution(ha_glmm) #  neg. binomial (zero-infl.) 56%

# negative binomial 1 
summary(ha_glmm_nb1)
r2(ha_glmm_nb1)
check_overdispersion(ha_glmm_nb1) # No overdispersion detected
check_zeroinflation(ha_glmm_nb1) # Model is underfitting zeros (probable zero-inflation)

# negative binomial 2
summary(ha_glmm_nb2)
r2(ha_glmm_nb2)
check_overdispersion(ha_glmm_nb2) # No overdispersion detected.
check_zeroinflation(ha_glmm_nb2) # Model is overfitting zeros.

# zero inflation 
summary(ha_glmm_zi)
r2(ha_glmm_zi)
check_overdispersion(ha_glmm_zi) # Overdispersion detected
check_zeroinflation(ha_glmm_zi) # Model is underfitting zeros (probable zero-inflation).

# negative binomial 2, zero inflation 
summary(ha_glmm_zi_nb2)
r2(ha_glmm_zi_nb2)
check_overdispersion(ha_glmm_zi_nb2) # No overdispersion detected.
check_zeroinflation(ha_glmm_zi_nb2) # Model is overfitting zeros.

# hurdle
summary(ha_glmm_hurdle)
r2(ha_glmm_hurdle)
check_overdispersion(ha_glmm_hurdle) # Overdispersion detected
check_zeroinflation(ha_glmm_hurdle) # Model is underfitting zeros (probable zero-inflation).

### Model comparison 

AICtab(ha_glmm,
       ha_glmm_nb1,
       ha_glmm_nb2,
       ha_glmm_zi,
       ha_glmm_zi_nb2,
       ha_glmm_hurdle,
       base=TRUE, weights=TRUE, logLik=TRUE)

# Best model, ha_glmm_nb2 <- glmmTMB(bats ~ treatment + (1|site), data = data, family = nbinom2(link = "log")) 
# Now, lets include the interaction: 
ha_glmm_nb2_intec <- glmmTMB(bats ~ treatment*site, data = data, family = nbinom2(link = "log"))
summary(ha_glmm_nb2_intec)
r2(ha_glmm_nb2_intec)
check_overdispersion(ha_glmm_nb2_intec) # No overdispersion detected.
check_zeroinflation(ha_glmm_nb2_intec) # Model is overfitting zeros.

AICtab(ha_glmm,
       ha_glmm_nb1,
       ha_glmm_nb2,
       ha_glmm_zi,
       ha_glmm_zi_nb2,
       ha_glmm_hurdle,
       ha_glmm_nb2_intec,
       base=TRUE, weights=TRUE, logLik=TRUE)

### High activity sites, different bat species

# Carollia
carollia_glmm <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, zi =~ treatment, family = poisson)
carollia_glmm_zi_nb2 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, zi =~ treatment, family = nbinom2(link = "log"))
carollia_glmm_hurdle <- glmmTMB(carollia_spp ~ treatment + (1|site), data = data, zi =~ treatment, family=truncated_poisson)

# poisson
summary(carollia_glmm)
r2(carollia_glmm)
check_overdispersion(carollia_glmm) # Overdispersion detected
check_zeroinflation(carollia_glmm) # Model is underfitting zeros (probable zero-inflation)
check_distribution(carollia_glmm) # poisson (zero-infl.) 38%

# negative binomial 1 
summary(carollia_glmm_nb1)
r2(carollia_glmm_nb1)
check_overdispersion(carollia_glmm_nb1) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb1) # Model is overfitting zeros.

# negative binomial 2
summary(carollia_glmm_nb2)
r2(carollia_glmm_nb2)
check_overdispersion(carollia_glmm_nb2) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.

# zero inflation 
summary(carollia_glmm_zi)
r2(carollia_glmm_zi)
check_overdispersion(carollia_glmm_zi) # No overdispersion detected.
check_zeroinflation(carollia_glmm_zi) # Model is underfitting zeros (probable zero-inflation).

# negative binomial 2, zero inflation 
summary(carollia_glmm_zi_nb2)
r2(carollia_glmm_zi_nb2)
check_overdispersion(carollia_glmm_zi_nb2) # No overdispersion detected.
check_zeroinflation(carollia_glmm_zi_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.

# hurdle
summary(carollia_glmm_hurdle)
r2(carollia_glmm_hurdle)
check_overdispersion(carollia_glmm_hurdle) # No overdispersion detected.
check_zeroinflation(carollia_glmm_hurdle) # Model is underfitting zeros (probable zero-inflation).

### Model comparison 

AICtab(carollia_glmm,
       carollia_glmm_nb1,
       carollia_glmm_nb2,
       carollia_glmm_zi,
       carollia_glmm_zi_nb2,
       carollia_glmm_hurdle,
       base=TRUE, weights=TRUE, logLik=TRUE)

### carollia_glmm_nb2 seems the best model

carollia_glmm_nb2_int <- glmmTMB(carollia_spp ~ treatment*site, data = data, family = nbinom2(link = "log"))
summary(carollia_glmm_nb2_int)
r2(carollia_glmm_nb2_int)
check_overdispersion(carollia_glmm_nb2_int) # No overdispersion detected
check_zeroinflation(carollia_glmm_nb2_int) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.


################################################
### data per low and high bat activity sites ###
################################################

### subset data
low_act <- data  %>% filter(site %in% c("site 1", "site 5", "site 4"))
high_act <- data  %>% filter(site %in% c("site 2", "site 3", "site 6"))

############################
### high activity sites ###
############################

# glmm with different distributions 
ha_glmm <- glmmTMB(bats ~ treatment + (1|site), data = high_act, family = poisson)
ha_glmm_nb1 <- glmmTMB(bats ~ treatment + (1|site), data = high_act, family = nbinom1(link = "log"))
ha_glmm_nb2 <- glmmTMB(bats ~ treatment + (1|site), data = high_act, family = nbinom2(link = "log"))
ha_glmm_zi <- glmmTMB(bats ~ treatment + (1|site), data = high_act, zi =~ treatment, family = poisson)
ha_glmm_zi_nb2 <- glmmTMB(bats ~ treatment + (1|site), data = high_act, zi =~ treatment, family = nbinom2(link = "log"))
ha_glmm_hurdle <- glmmTMB(bats ~ treatment + (1|site), data = high_act, zi =~ treatment, family=truncated_poisson)

# poisson
summary(ha_glmm)
r2(ha_glmm)
check_overdispersion(ha_glmm) # Overdispersion detected
check_zeroinflation(ha_glmm) # Model is underfitting zeros (probable zero-inflation)
check_distribution(ha_glmm) #  neg. binomial (zero-infl.) 72%

# negative binomial 1 
summary(ha_glmm_nb1)
r2(ha_glmm_nb1)
check_overdispersion(ha_glmm_nb1) # Overdispersion detected
check_zeroinflation(ha_glmm_nb1) # Model is underfitting zeros (probable zero-inflation)

# negative binomial 2
summary(ha_glmm_nb2)
r2(ha_glmm_nb2)
check_overdispersion(ha_glmm_nb2) # No overdispersion detected.
check_zeroinflation(ha_glmm_nb2) # Model is overfitting zeros.

# zero inflation 
summary(ha_glmm_zi)
r2(ha_glmm_zi)
check_overdispersion(ha_glmm_zi) # Overdispersion detected
check_zeroinflation(ha_glmm_zi) # Model is underfitting zeros (probable zero-inflation).

# negative binomial 2, zero inflation 
summary(ha_glmm_zi_nb2)
r2(ha_glmm_zi_nb2)
check_overdispersion(ha_glmm_zi_nb2) # No overdispersion detected.
check_zeroinflation(ha_glmm_zi_nb2) # Model is overfitting zeros.

# hurdle
summary(ha_glmm_hurdle)
r2(ha_glmm_hurdle)
check_overdispersion(ha_glmm_hurdle) # Overdispersion detected
check_zeroinflation(ha_glmm_hurdle) # Model is underfitting zeros (probable zero-inflation).

### Model comparison 
AICtab(ha_glmm,
       ha_glmm_nb1,
       ha_glmm_nb2,
       ha_glmm_zi,
       ha_glmm_zi_nb2,
       ha_glmm_hurdle,
       base=TRUE, weights=TRUE, logLik=TRUE)
              

### High activity sites, different bat species

# Carollia
carollia_glmm <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, family = poisson)
carollia_glmm_nb1 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, zi =~ treatment, family = poisson)
carollia_glmm_zi_nb2 <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, zi =~ treatment, family = nbinom2(link = "log"))
carollia_glmm_hurdle <- glmmTMB(carollia_spp ~ treatment + (1|site), data = high_act, zi =~ treatment, family=truncated_poisson)

# poisson
summary(carollia_glmm)
r2(carollia_glmm)
check_overdispersion(carollia_glmm) # Overdispersion detected
check_zeroinflation(carollia_glmm) # Model is underfitting zeros (probable zero-inflation)
check_distribution(carollia_glmm) # negative binomial 56%

# negative binomial 1 
summary(carollia_glmm_nb1)
r2(carollia_glmm_nb1)
check_overdispersion(carollia_glmm_nb1) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb1) # Model is overfitting zeros.

# negative binomial 2
summary(carollia_glmm_nb2)
r2(carollia_glmm_nb2)
check_overdispersion(carollia_glmm_nb2) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.

# zero inflation 
summary(carollia_glmm_zi)
r2(carollia_glmm_zi)
check_overdispersion(carollia_glmm_zi) # No overdispersion detected.
check_zeroinflation(carollia_glmm_zi) # Model is underfitting zeros (probable zero-inflation).

# negative binomial 2, zero inflation 
summary(carollia_glmm_zi_nb2)
r2(carollia_glmm_zi_nb2)
check_overdispersion(carollia_glmm_zi_nb2) # No overdispersion detected.
check_zeroinflation(carollia_glmm_zi_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.

# hurdle
summary(carollia_glmm_hurdle)
r2(carollia_glmm_hurdle)
check_overdispersion(carollia_glmm_hurdle) # No overdispersion detected.
check_zeroinflation(carollia_glmm_hurdle) # Model is underfitting zeros (probable zero-inflation).

### Model comparison 
AICtab(carollia_glmm,
       carollia_glmm_nb1,
       carollia_glmm_nb2,
       carollia_glmm_zi,
       carollia_glmm_zi_nb2,
       carollia_glmm_hurdle,
       base=TRUE, weights=TRUE, logLik=TRUE)






