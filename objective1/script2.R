#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### script 2: BAT SPP. GLMMs ###

rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(effects)
library(glmmTMB)
library(sjPlot)
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
library(bbmle) 

### Load data: bats captured using mist nets 
setwd("~/Desktop/lures/objective1") # set directory to objective 1 folder
data <- read.csv("data.csv")
head(data)

### Describe data, get N per control/treatment
table(data$treatment) # control = 28 ; lures = 27
table(interaction(data$treatment, data$site))

### Group fruit bats by date and site
total_bats <- data %>%
  group_by(date, treatment, site) %>%
  summarise(total_bats = sum(across(matches("fruit_bats"))))
head(total_bats)

fruitbats_glmm <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = nbinom2(link = "log"))
fruitbats_glmm_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = poisson)
fruitbats_glmm_nb1_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = nbinom1(link = "log"))
fruitbats_glmm_nb2_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = nbinom2(link = "log"))
fruitbats_glmm_hurdle_poisson_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = truncated_poisson)
fruitbats_glmm_hurdle_bn2_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family=truncated_nbinom2)
fruitbats_glmm_zi_site <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = poisson)
fruitbats_glmm_null <- glmmTMB(total_bats ~ 1  + (1|site) + (1|date), data =  total_bats, family = poisson)
# Model selection
AIC_total_bats <- AICtab(fruitbats_glmm_null,
                         fruitbats_glmm,
                         fruitbats_glmm_nb1,
                         fruitbats_glmm_nb2,
                         fruitbats_glmm_zi,
                         fruitbats_glmm_nb1_zi,
                         fruitbats_glmm_nb2_zi,
                         fruitbats_glmm_hurdle_poisson_zi,
                         fruitbats_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_total_bats
write.csv(AIC_total_bats, "AIC_total_bats.csv")
tab_model(fruitbats_glmm_hurdle_bn2_zi) ### tab best model: # negative binomial 2
parameters(fruitbats_glmm_hurdle_bn2_zi)
r2(fruitbats_glmm_hurdle_bn2_zi)

glmm_emmeans <-emmeans(fruitbats_glmm_hurdle_bn2_zi, ~ treatment, type="response")
glmm_emmeans
glmm_emmeans <- as.data.frame(glmm_emmeans)

dwtest(fruitbats_glmm_hurdle_bn2_zi) # Durbin-Watson test
# First, store the residuals
residuals <- resid(fruitbats_glmm_hurdle_bn2_zi)
# Calculate ACF
acf_data <- acf(residuals, plot = FALSE)
# Convert ACF data to dataframe
acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
# Calculate confidence intervals
ci_upper <- qnorm(0.975) / sqrt(length(residuals))
ci_lower <- -qnorm(0.975) / sqrt(length(residuals))

ACF_fruitbats <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", width = 0.1) +
  geom_hline(yintercept = c(ci_upper, ci_lower, 0), linetype = c("dashed", "dashed", "solid"), color = c("darkgrey", "darkgrey", "black")) +
  labs(title = "ACF Fruit bats",
       x = "Lag", y = "Autocorrelation")+ 
  theme_test(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  geom_text(x = 12, y = 0.8, label = "  DW = 1.857,\nP-value = 0.291", size = 2.5)

  
ACF_fruitbats

# Check for overdispesion and zero inflation 
summary(fruitbats_glmm_hurdle_bn2_zi)
check_overdispersion(fruitbats_glmm_hurdle_bn2_zi) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_hurdle_bn2_zi) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(fruitbats_glmm_hurdle_bn2_zi)) # looks like no interaction 

total_bats$prediction <- predict(fruitbats_glmm_hurdle_bn2_zi, total_bats, type="response")
fruitbatstotal <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_test(base_size = 15) +  
  geom_point(data = total_bats, aes(x = treatment, y = total_bats, color = site), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical\nlures")) + 
  ylab ("Total count of fruit bats") +
  xlab (" ") +
  theme(legend.position = "none") 
fruitbatstotal

ggsave(file="Figure4.jpg", 
       plot= fruitbatstotal,
       width=7,height=10,units="cm",dpi=500)

### Pivot longer 
data_long <- data %>%
  pivot_longer(cols = -c(date, site, treatment, hours, nets),
               names_to = "bat_species", values_to = "abundance") 
data_long <- data_long  %>% filter(!bat_species %in% c("bats", "fruit_bats", "cperspicillata", "csowelli", "ccastanea")) # eliminate redundant data
head(data_long)
histogram(data_long$abundance) # lots of zeros! 

### Subset the dataset by bat species 

# carollia
data_long
carollia_data <- data_long %>% filter(bat_species %in% c("carollia_spp"))

# models 
carollia_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = carollia_data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = poisson)
carollia_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = truncated_poisson)
carollia_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = carollia_data, family=truncated_nbinom2)
carollia_glmm_null <- glmmTMB(abundance ~ 1, data = carollia_data, family = poisson)

# compare all models
check_distribution(carollia_glmm) # neg. binomial (zero-infl.) 56%

AIC_carollia <- AICtab(carollia_glmm_null,
                       carollia_glmm,
                       carollia_glmm_nb1,
                       #carollia_glmm_nb2,
                       carollia_glmm_zi,
                       #carollia_glmm_nb1_zi,
                       carollia_glmm_nb2_zi,
                       carollia_glmm_hurdle_poisson_zi,
                       carollia_glmm_hurdle_bn2_zi,
                       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_carollia, "AIC_carollia.csv")

tab_model(carollia_glmm) ### tab best model
parameters(carollia_glmm)
r2(carollia_glmm)
check_model(carollia_glmm)
# Check for overdispesion and zero inflation 
summary(carollia_glmm)
check_overdispersion(carollia_glmm) # No overdispersion detected.
check_zeroinflation(carollia_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(carollia_glmm))

# calculate effect size
glmm_emmeans <- emmeans(carollia_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# control 0.29, lures 0.60
0.60 - 0.29 # 0.31
(0.31*100)/0.29 # 106.89
0.60/0.29

dwtest(carollia_glmm)
# store the residuals
residuals <- resid(carollia_glmm)
# Calculate ACF
acf_data <- acf(residuals, plot = FALSE)
# Convert ACF data to dataframe
acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
# Calculate confidence intervals
ci_upper <- qnorm(0.975) / sqrt(length(residuals))
ci_lower <- -qnorm(0.975) / sqrt(length(residuals))

ACF_carollia <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", width = 0.1) +
  geom_hline(yintercept = c(ci_upper, ci_lower, 0), linetype = c("dashed", "dashed", "solid"), color = c("darkgrey", "darkgrey", "black")) +
  labs(title = expression("ACF" * italic(" Carollia")~"spp."),
       x = "Lag", y = " ") + 
  theme_test(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  geom_text(x = 12, y = 0.8, label = "  DW = 2.066,\nP-value = 0.605", size = 2.5)
ACF_carollia

# compute predictions for graph using the function predict
carollia_text <- parse(text = "italic('Carollia') ~ 'spp.'")
glmm_graph_carollia <-  ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_test(base_size = 15) +  
  geom_point(data = carollia_data, aes(x = treatment, y = abundance, color = site), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical\nlures")) +
  ylab (expression("Total count of " * italic("Carollia") ~ 'spp.')) +
  xlab ("") + 
  geom_text(x = 1.5, y = 5, label = "*", size = 10)
glmm_graph_carollia
ggsave(file="carollia.jpg", 
       plot= glmm_graph_carollia,
       width=15,height=15,units="cm",dpi=300)

# Uroderma
uroderma_data <- data_long %>% filter(bat_species %in% c("uroderma_spp"))

# models 
uroderma_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = uroderma_data, family = poisson)
uroderma_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = poisson)
uroderma_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = truncated_poisson)
uroderma_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = uroderma_data, family=truncated_nbinom2)
uroderma_glmm_null <- glmmTMB(abundance ~ 1, data = uroderma_data, family = poisson)
# compare all models
check_distribution(uroderma_glmm) # neg. binomial (zero-infl.) 56%

AIC_uroderma <- AICtab(uroderma_glmm_null,
                       uroderma_glmm,
                       uroderma_glmm_nb1,
                       uroderma_glmm_nb2,
                       uroderma_glmm_zi,
                       #uroderma_glmm_nb1_zi,
                       #uroderma_glmm_nb2_zi,
                       uroderma_glmm_hurdle_poisson_zi,
                       #uroderma_glmm_hurdle_bn2_zi,
                       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_uroderma, "AIC_uroderma.csv")

tab_model(uroderma_glmm_hurdle_poisson_zi) ### tab best model
parameters(uroderma_glmm_hurdle_poisson_zi)
r2(uroderma_glmm_hurdle_poisson_zi)
# Check for overdispesion and zero inflation 
summary(uroderma_glmm_hurdle_poisson_zi)
check_overdispersion(uroderma_glmm_hurdle_poisson_zi) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_hurdle_poisson_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(uroderma_glmm))
check_model(uroderma_glmm)

# calculate effect size
glmm_emmeans <-emmeans(uroderma_glmm_hurdle_poisson_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

dwtest(uroderma_glmm)
# store the residuals
residuals <- resid(uroderma_glmm)
# Calculate ACF
acf_data <- acf(residuals, plot = FALSE)
# Convert ACF data to dataframe
acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
# Calculate confidence intervals
ci_upper <- qnorm(0.975) / sqrt(length(residuals))
ci_lower <- -qnorm(0.975) / sqrt(length(residuals))

ACF_uroderma <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", width = 0.1) +
  geom_hline(yintercept = c(ci_upper, ci_lower, 0), linetype = c("dashed", "dashed", "solid"), color = c("darkgrey", "darkgrey", "black")) +
  labs(title = expression("ACF" * italic(" Uroderma")~"spp."),
       x = "Lag", y = " ") + 
  theme_test(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  geom_text(x = 12, y = 0.8, label = "  DW = 2.205,\nP-value = 0.784", size = 2.5)
ACF_uroderma


uroderma_data$prediction_uroderma <- predict(uroderma_glmm_hurdle_poisson_zi, uroderma_data, type="response")
glmm_graph_uroderma <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_test(base_size = 15) +  
  geom_point(data = uroderma_data, aes(x = treatment, y = abundance, color = site), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  #geom_jitter(data = uroderma_data, aes(x = treatment, y = abundance, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical\nlures")) + 
  ylab (expression("Total count of " * italic("Uroderma") ~ 'spp.')) +
  xlab ("") 
glmm_graph_uroderma
ggsave(file="uroderma.jpg", 
       plot= glmm_graph_uroderma,
       width=15,height=15,units="cm",dpi=300)

# Ectophylla alba
ectophylla.alba_data <- data_long %>% filter(bat_species %in% c("ectophylla.alba"))
# GLMMs
ectophylla.alba_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = ectophylla.alba_data, family = poisson)
ectophylla.alba_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = ectophylla.alba_data, family = nbinom1(link = "log"))
ectophylla.alba_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = ectophylla.alba_data, family = nbinom2(link = "log"))
ectophylla.alba_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = poisson)
ectophylla.alba_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = nbinom1(link = "log"))
ectophylla.alba_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = nbinom2(link = "log"))
ectophylla.alba_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = truncated_poisson)
ectophylla.alba_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = ectophylla.alba_data, family=truncated_nbinom2)
ectophylla.alba_glmm_null <- glmmTMB(abundance ~ 1, data = ectophylla.alba_data, family = poisson)
# Model selection
AIC_ectophylla <- AICtab(ectophylla.alba_glmm_null,
                         ectophylla.alba_glmm,
                         ectophylla.alba_glmm_nb1,
                         ectophylla.alba_glmm_nb2,
                         #ectophylla.alba_glmm_zi,
                         #ectophylla.alba_glmm_nb1_zi,
                         #ectophylla.alba_glmm_nb2_zi,
                         ectophylla.alba_glmm_hurdle_poisson_zi,
                         ectophylla.alba_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_ectophylla, "AIC_ectophylla.csv")
tab_model(ectophylla.alba_glmm_hurdle_bn2_zi) ### tab best model
parameters(ectophylla.alba_glmm_hurdle_bn2_zi)
r2(ectophylla.alba_glmm_hurdle_bn2_zi)
# Check for overdispesion and zero inflation 
summary(ectophylla.alba_glmm_hurdle_bn2_zi)
check_overdispersion(ectophylla.alba_glmm_hurdle_bn2_zi) # No overdispersion detected.
check_zeroinflation(ectophylla.alba_glmm_hurdle_bn2_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(ectophylla.alba_glmm_hurdle_bn2_zi))
check_model(ectophylla.alba_glmm_hurdle_bn2_zi)
# calculate effect size
glmm_emmeans <-emmeans(ectophylla.alba_glmm_hurdle_bn2_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)

dwtest(ectophylla.alba_glmm_hurdle_bn2_zi)
# store the residuals
residuals <- resid(ectophylla.alba_glmm_hurdle_bn2_zi)
# Calculate ACF
acf_data <- acf(residuals, plot = FALSE)
# Convert ACF data to dataframe
acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
# Calculate confidence intervals
ci_upper <- qnorm(0.975) / sqrt(length(residuals))
ci_lower <- -qnorm(0.975) / sqrt(length(residuals))

ACF_ectophylla.alba <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", width = 0.1) +
  geom_hline(yintercept = c(ci_upper, ci_lower, 0), linetype = c("dashed", "dashed", "solid"), color = c("darkgrey", "darkgrey", "black")) +
  labs(title = expression("ACF" * italic(" E. alba")),
       x = "Lag", y = " ") + 
  theme_test(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  geom_text(x = 12, y = 0.8, label = "  DW = 2.110,\nP-value = 0.667", size = 2.5)

ACF_ectophylla.alba

# compute predictions for graph using the function predict
ectophylla.alba_data$prediction_ectophylla.alba <- predict(ectophylla.alba_glmm_hurdle_bn2_zi, ectophylla.alba_data, type="response")
glmm_graph_ectophylla.alba <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_test(base_size = 15) +  
  geom_point(data = ectophylla.alba_data, aes(x = treatment, y = abundance, color = site), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical\nlures")) + 
  ylab (expression("Total count of " * italic("Ectophylla alba"))) +
  xlab ("") 
glmm_graph_ectophylla.alba
ggsave(file="ectophylla.alba.jpg", 
       plot= glmm_graph_ectophylla.alba,
       width=15,height=15,units="cm",dpi=300)

sp_bats <- ggarrange(fruitbatstotal,
                     glmm_graph_carollia,
                     glmm_graph_ectophylla.alba,
                     glmm_graph_uroderma,
                     ncol = 2,
                     nrow = 2,
                     align = "hv",
                     common.legend = TRUE,
                     legend="bottom")
sp_bats

ggsave(file="Figure5.jpg", 
       plot= sp_bats,
       width=15,height=18,units="cm",dpi=600)

#### Temporal variation ###
data$date <- as.Date(data$date, format = "%B/%d")
data$days_since_reference <- as.numeric(data$date - min(data$date))
head(data)

temporal <- glmmTMB(fruit_bats ~ days_since_reference + (1|treatment), data = data, family = poisson)
summary(temporal)

temporal_variation <- ggplot(data, aes(x = days_since_reference, y = fruit_bats, color = site, shape = treatment)) +
  theme_test(base_size = 10) +
  geom_point(size = 1.5, color = "darkgrey") +  
  geom_point(size = 1) + 
  scale_shape_manual(values = c("control" = "circle", "lures" = "triangle"),
                     name = "Treatment", labels = c("Control", "Chemical\nlures")) + 
  scale_color_viridis_d(option = "D", name = "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  xlab("Capture night") +
  ylab ("Total fruit bats") +
  scale_y_continuous(breaks = pretty_breaks(n = 2))
             
temporal_variation 

temporal_variation_carollia <- ggplot(data, aes(x = days_since_reference, y = carollia_spp, color = site, shape = treatment)) +
  theme_test(base_size = 10) +
  geom_point(size = 1.5, color = "darkgrey") +  
  geom_point(size = 1) + 
  scale_shape_manual(values = c("control" = "circle", "lures" = "triangle"),
                     name = "Treatment", labels = c("Control", "Chemical\nlures")) + 
  scale_color_viridis_d(option = "D", name = "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  xlab("Capture night") +
  ylab (expression(italic("Carollia") ~ 'spp.')) +
  scale_y_continuous(breaks = pretty_breaks(n = 2))

temporal_variation_carollia

temporal_variation_ectophylla <- ggplot(data, aes(x = days_since_reference, y = ectophylla.alba, color = site, shape = treatment)) +
  theme_test(base_size = 10) +
  geom_point(size = 1.5, color = "darkgrey") +  
  geom_point(size = 1) + 
  scale_shape_manual(values = c("control" = "circle", "lures" = "triangle"),
                     name = "Treatment", labels = c("Control", "Chemical\nlures")) + 
  scale_color_viridis_d(option = "D", name = "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  xlab("Capture night") +
  ylab (expression(italic("E. alba"))) + 
  scale_y_continuous(breaks = pretty_breaks(n = 2))

temporal_variation_ectophylla

temporal_variation_uroderma <- ggplot(data, aes(x = days_since_reference, y = uroderma_spp, color = site, shape = treatment)) +
  theme_test(base_size = 10) +
  geom_point(size = 1.5, color = "darkgrey") +  
  geom_point(size = 1) + 
  scale_shape_manual(values = c("control" = "circle", "lures" = "triangle"),
                     name = "Treatment", labels = c("Control", "Chemical\nlures")) + 
  scale_color_viridis_d(option = "D", name = "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  xlab("Capture night") +
  ylab (expression(italic("Uroderma") ~ 'spp.')) + 
  scale_y_continuous(breaks = pretty_breaks(n = 1))

temporal_variation_uroderma

temporal_all <- ggarrange(temporal_variation,
                          temporal_variation_carollia,
                          temporal_variation_ectophylla,
                          temporal_variation_uroderma,
                          ncol = 4,
                          nrow = 1,
                          align = "hv",
                          common.legend = TRUE)


temporal_ACF <- ggarrange(ACF_fruitbats,
                          ACF_carollia,
                          ACF_ectophylla.alba,
                          ACF_uroderma,
                          ncol = 4,
                          nrow = 1,
                          align = "hv")

final_temporal_ACF <- ggarrange(temporal_all,
                                 temporal_ACF,
                                 ncol = 1,
                                 nrow = 2,
                                 align = "v")
                                 
final_temporal_ACF


ggsave(file="final_temporal_ACF.jpg", 
       plot= final_temporal_ACF,
       width=7,height=5,units="in",dpi=600)


#######################################
#######################################
#######################################
#######################################