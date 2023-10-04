######################################################
### OBJECTIVE 2: SEED RAIN WITH AND WITHOUT LURES ###
######################################################

### script 3: SEEDS COMMUNITY AND SEED GLMMs ###

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
library(scales)
library(sjPlot)
library(ggvenn)

### read and clean data ###
setwd("~/Desktop/lures/objective2")  # set directory to objective 2 folder
data <- read.csv("seed_data.csv")
data$week <- as.factor(data$week)
data$site_letter <- as.factor(data$site_letter)
data$treatment <- as.factor(data$treatment)
data_cleaned <- data[, c(1, 3, 5, 6:27)]
head(data_cleaned)

### Describe data, get N 
table(data$treatment) # basal = 12; control = 8 ; lures = 8
table(interaction(data$treatment, data$site_name))

### Objective 2.1: characterize the seed rain in the study area

# venn diagram 
species_list <- data_cleaned %>%
  pivot_longer(cols = -c(week, site_letter, treatment), names_to = "species", values_to = "value") %>%
  group_by(treatment) %>%
  filter(value > 0) %>%
  summarise(species_present = list(unique(species))) %>%
  pull(species_present)
names(species_list) <- c("Baseline", "Control", "Chemical lures")

viridis(3)

colors <- c("#21908CFF", "#FDE725FF", "#440154FF")
venn_diagram <- ggvenn(
  species_list,
  fill_color = colors,
  stroke_color = "darkgrey",
  stroke_size = 0.7, set_name_size = 5
)
venn_diagram
ggsave(file="Figure8.jpg",
       plot= venn_diagram,
       width=16,height=16,units="cm",dpi=300)

# pivot longer

data_longer <- data_cleaned %>%
  pivot_longer(cols = -c(week, site_letter, treatment),
               names_to = "plant_species",
               values_to = "seed_count")
head(data_longer)
filtered_df <- data_longer %>% filter(seed_count > 0)
head(filtered_df)
supp.labs <- c("Control", "Chemical lures", "Baseline")
names(supp.labs) <- c("control", "treatment", "basal")

allseeds_total <- ggplot(filtered_df, aes(x = site_letter , y = seed_count, fill = plant_species)) +
  theme_classic(base_size = 25) +
  theme(legend.text = element_text(size = 15)) + 
  theme(panel.border = element_blank(), strip.background = element_rect(colour = "NA", fill = "NA")) +
  scale_fill_viridis_d(direction = -1, name = "Seed ID",
                       labels = c("unknown_1" = "Unknown 1",
                                  "unknown_2" = "Unknown 2",
                                  "cecropia_insignis" = parse(text = "italic('Cecropia insignis')"),
                                  "unknown_3" = "Unknown 3",
                                  "unknown_4" = "Unknown 4",
                                  "vismia_sp" = parse(text = "italic('Vismia') ~'sp.'"),
                                  "piper_multiplinervium" = parse(text = "italic('Piper multiplinervium')"),
                                  "cecropia_obtusifolia" = parse(text = "italic('Cecropia obtusifolia')"),
                                  "piper_santifelicis" = parse(text = "italic('Piper santifelicis')"),
                                  "solanum_sp" = parse(text = "italic('Solanum') ~ 'sp.'"),
                                  "anturium_sp" = parse(text = "italic('Anturium') ~'sp.'"),
                                  "ficus_columbrinae" = parse(text = "italic('Ficus columbrinae')"),
                                  "poaceae" = "Poaceae",
                                  "paspalum_conjugatum" = parse(text = "italic('Paspalum conjugatum')"),
                                  "piper_umbricola" = parse(text = "italic('Piper umbricola')"),
                                  "unknown_6" = "Unknown 6",
                                  "achenio_sp" = parse(text = "italic('Achenio') ~'sp.'"),
                                  "solanaceae" = "Solanaceae",
                                  "unknown_5" = "Unknown 5",
                                  "unknown_7" = "Unknown 7",
                                  "unknown_8" = "Unknown 8",
                                  "ficus_insipida" = parse(text = "italic('Ficus insipida')"),
                                  "unknown_9" = "Unknown 9")) +
  geom_bar(stat = "identity") +
  ylab("Total seed number") +
  xlab("") +
  theme(legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

allseeds_total
ggsave(file="allseeds_total.jpg", 
       plot= allseeds_total,
       width=18,height=16,units="cm",dpi=300)

#######################################
### ORGANIZE DATA IN PLANT FAMILIES ###
#######################################
data_family_filtered <- data_cleaned %>%
  rowwise() %>%
  mutate(Piperaceae = piper_multiplinervium + piper_umbricola + piper_santifelicis,
         Urticaceae = cecropia_insignis + cecropia_obtusifolia,
         Hypericaceae = vismia_sp,
         Poaceae = paspalum_conjugatum + poaceae,
         Moraceae = ficus_columbrinae + ficus_insipida,
         Araceae = anturium_sp,
         Solanaceae = solanaceae + solanum_sp,
         Unknown = sum(c_across(starts_with("unknown")))) %>%
  filter(!treatment %in% "basal") %>%
  select_if(function(col) any(col != 0))
head(data_family_filtered)
data_family <- data_family_filtered[, c(1:3, 18:23)]
data_family$total <- data_family$Piperaceae + data_family$Urticaceae + data_family$Poaceae + data_family$Moraceae + data_family$Solanaceae + data_family$Unknown

### SUMMARY OF TOTAL SEEDS COLLETED IN CONTROL AND TREATMENT ###
head(data_family)
data_family$species_richness <- rowSums(data_family[,4:9] > 0)

total_seed_count_per_family <- data_family %>%
  pivot_longer(cols = c(Piperaceae, Urticaceae, Poaceae, Moraceae, Solanaceae, Unknown, total),
               names_to = "Plant_family",
               values_to = "value") %>%
  group_by(treatment, Plant_family) %>%
  summarize(total_seed_count = sum(value)) %>%
  spread(treatment, total_seed_count, sep = "_")
total_seed_count_per_family <- as.data.frame(total_seed_count_per_family)
total_seed_count_per_family
write.csv(total_seed_count_per_family, file = "total_seed_count_per_family.csv")

### TOTAL NUMBER OF SEEDS ###
# GLMMs
total_glmm <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
total_glmm_nb1 <- glmmTMB(total ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
total_glmm_nb2 <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
total_glmm_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
total_glmm_nb1_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
total_glmm_nb2_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
total_glmm_hurdle_poisson_zi  <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
total_glmm_hurdle_bn2_zi <- glmmTMB(total ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
total_glmm_null <- glmmTMB(total ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_total <- AICtab(total_glmm_null,
                    total_glmm, 
                    total_glmm_nb1,
                    total_glmm_nb2,
                    total_glmm_zi,
                    #total_glmm_nb1_zi,
                    total_glmm_nb2_zi,
                    total_glmm_hurdle_poisson_zi,
                    total_glmm_hurdle_bn2_zi,
                    base=TRUE, weights=TRUE, logLik=TRUE)
AIC_total
write.csv(AIC_total, "AIC_total.csv")
tab_model(total_glmm_zi) ### tab best model
parameters(total_glmm_zi)
r2(total_glmm_zi)

# Check for overdispesion and zero inflation 
summary(total_glmm_zi)
check_overdispersion(total_glmm_zi) # No overdispersion detected.
check_zeroinflation(total_glmm_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(total_glmm_zi))
check_model(total_glmm_zi)

# calculate effect size
glmm_emmeans <-emmeans(total_glmm_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# control 31.38, lures 72.64
72.64 - 31.38 # 41.26
(41.26*100)/31.38 # 131.485
72.64/31.38

head(data_family)
color_palette <- viridis(20, option = "G")
selected_colors <- color_palette[c(5, 10, 15, 18)]

glmm_graph_total <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = total, color = site_letter, size = species_richness), position = position_jitterdodge(dodge.width = 0.4), alpha = 0.7) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Total seed count") +
  xlab ("") +
  theme(legend.position = "right") +
  geom_text(x = 1.5, y = 700, label = "***", size = 5) + 
  scale_size(name = "Family\nrichness")
glmm_graph_total
ggsave(file="totalseed_glmm.jpg", 
       plot= glmm_graph_total,
       width=12,height=12,units="cm",dpi=500)

### PIPERACEAE ###
# GLMMs
piperaceae_glmm <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
piperaceae_glmm_nb1 <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
piperaceae_glmm_nb2 <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
piperaceae_glmm_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
piperaceae_glmm_nb1_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
piperaceae_glmm_nb2_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
piperaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
piperaceae_glmm_hurdle_bn2_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
piperaceae_glmm_null <- glmmTMB(Piperaceae ~ 1, data = data_family, family = poisson)
# Model comparison 
AIC_piperaceae <- AICtab(piperaceae_glmm_null,
                         piperaceae_glmm, 
                         piperaceae_glmm_nb1,
                         #piperaceae_glmm_nb2,
                         piperaceae_glmm_zi,
                         #piperaceae_glmm_nb1_zi,
                         #piperaceae_glmm_nb2_zi,
                         piperaceae_glmm_hurdle_poisson_zi,
                         piperaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_piperaceae
write.csv(AIC_piperaceae, "AIC_piperaceae.csv")
tab_model(piperaceae_glmm) ### tab best model
parameters(piperaceae_glmm)
r2(piperaceae_glmm)

# Check for overdispesion and zero inflation 
summary(piperaceae_glmm)
check_overdispersion(piperaceae_glmm) # No overdispersion detected.
check_zeroinflation(piperaceae_glmm) # Model seems ok 
plot(allEffects(piperaceae_glmm))
check_model(piperaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(piperaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)

glmm_graph_piperaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = Piperaceae, color = site_letter), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Piperaceae") +
  xlab ("")

glmm_graph_piperaceae
ggsave(file="piperaceae.jpg", 
       plot= glmm_graph_piperaceae,
       width=15,height=15,units="cm",dpi=300)

### URTUCACEAE ###
# GLMMs
Urticaceae_glmm <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Urticaceae_glmm_nb1 <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Urticaceae_glmm_nb2 <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Urticaceae_glmm_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Urticaceae_glmm_nb1_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Urticaceae_glmm_nb2_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Urticaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Urticaceae_glmm_hurdle_bn2_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Urticaceae_glmm_null <- glmmTMB(Urticaceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Urticaceae <- AICtab(Urticaceae_glmm_null,
                         Urticaceae_glmm, 
                         Urticaceae_glmm_nb1,
                         Urticaceae_glmm_nb2,
                         #Urticaceae_glmm_zi,
                         #Urticaceae_glmm_nb1_zi,
                         #Urticaceae_glmm_nb2_zi,
                         #Urticaceae_glmm_hurdle_poisson_zi,
                         #Urticaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Urticaceae
write.csv(AIC_Urticaceae, "AIC_Urticaceae.csv")
tab_model(Urticaceae_glmm_nb2) ### tab best model
parameters(Urticaceae_glmm_nb2)
r2(Urticaceae_glmm_nb2)

# Check for overdispesion and zero inflation 
summary(Urticaceae_glmm_nb2)
check_overdispersion(Urticaceae_glmm_nb2) # No overdispersion detected.
check_zeroinflation(Urticaceae_glmm_nb2) # Model seems ok 
plot(allEffects(Urticaceae_glmm_nb2))
check_model(Urticaceae_glmm_nb2)

# calculate effect size
glmm_emmeans <-emmeans(Urticaceae_glmm_nb2,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans$treatment
glmm_emmeans$response

glmm_graph_Urticaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = Urticaceae, color = site_letter), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Urticaceae") +
  xlab ("")

glmm_graph_Urticaceae

ggsave(file="Urticaceae.jpg", 
       plot= glmm_graph_Urticaceae,
       width=15,height=15,units="cm",dpi=300)

### POACEAE ###
# GLMMs
Poaceae_glmm <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Poaceae_glmm_nb1 <- glmmTMB(Poaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Poaceae_glmm_nb2 <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Poaceae_glmm_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Poaceae_glmm_nb1_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Poaceae_glmm_nb2_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Poaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Poaceae_glmm_hurdle_bn2_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Poaceae_glmm_null <- glmmTMB(Poaceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Poaceae <- AICtab(Poaceae_glmm_null,
                         Poaceae_glmm, 
                         Poaceae_glmm_nb1,
                         Poaceae_glmm_nb2,
                         Poaceae_glmm_zi,
                         #Poaceae_glmm_nb1_zi,
                         #Poaceae_glmm_nb2_zi,
                         Poaceae_glmm_hurdle_poisson_zi,
                         Poaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Poaceae
write.csv(AIC_Poaceae, "AIC_Poaceae.csv")
tab_model(Poaceae_glmm) ### tab best model
parameters(Poaceae_glmm)
r2(Poaceae_glmm)

# Check for overdispesion and zero inflation 
summary(Poaceae_glmm)
check_overdispersion(Poaceae_glmm) # No overdispersion detected.
check_zeroinflation(Poaceae_glmm) # Model seems ok 
plot(allEffects(Poaceae_glmm))
check_model(Poaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Poaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans$treatment
glmm_emmeans$rate

glmm_graph_Poaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = Poaceae, color = site_letter), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Poaceae") +
  xlab ("")

glmm_graph_Poaceae
ggsave(file="Poaceae.jpg", 
       plot= glmm_graph_Poaceae,
       width=15,height=15,units="cm",dpi=300)

### MORACEAE ###
# GLMMs
Moraceae_glmm <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Moraceae_glmm_nb1 <- glmmTMB(Moraceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Moraceae_glmm_nb2 <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Moraceae_glmm_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Moraceae_glmm_nb1_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Moraceae_glmm_nb2_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Moraceae_glmm_hurdle_poisson_zi  <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Moraceae_glmm_hurdle_bn2_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Moraceae_glmm_null <- glmmTMB(Moraceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Moraceae <- AICtab(Moraceae_glmm_null,
                      Moraceae_glmm, 
                      Moraceae_glmm_nb1,
                      Moraceae_glmm_nb2,
                      #Moraceae_glmm_zi,
                      #Moraceae_glmm_nb1_zi,
                      #Moraceae_glmm_nb2_zi,
                      Moraceae_glmm_hurdle_poisson_zi,
                      #Moraceae_glmm_hurdle_bn2_zi,
                      base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Moraceae
write.csv(AIC_Moraceae, "AIC_Moraceae.csv")
tab_model(Moraceae_glmm) ### tab best model
parameters(Moraceae_glmm)
r2(Moraceae_glmm)

# Check for overdispesion and zero inflation 
summary(Moraceae_glmm)
check_overdispersion(Moraceae_glmm) # No overdispersion detected.
check_zeroinflation(Moraceae_glmm) # Model seems ok 
plot(allEffects(Moraceae_glmm))
check_model(Moraceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Moraceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

glmm_graph_Moraceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = Moraceae, color = site_letter), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Moraceae") +
  xlab ("")

glmm_graph_Moraceae
ggsave(file="Moraceae.jpg", 
       plot= glmm_graph_Moraceae,
       width=15,height=15,units="cm",dpi=300)

### SOLANACEAE ###
# GLMMs
Solanaceae_glmm <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Solanaceae_glmm_nb1 <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Solanaceae_glmm_nb2 <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Solanaceae_glmm_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Solanaceae_glmm_nb1_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Solanaceae_glmm_nb2_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Solanaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Solanaceae_glmm_hurdle_bn2_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Solanaceae_glmm_null <- glmmTMB(Solanaceae ~ 1, data = data_family, family = poisson)
# Model comparison 
AIC_Solanaceae <- AICtab(Solanaceae_glmm_null,
                         Solanaceae_glmm, 
                         Solanaceae_glmm_nb1,
                         Solanaceae_glmm_nb2,
                         #Solanaceae_glmm_zi,
                         #Solanaceae_glmm_nb1_zi,
                         #Solanaceae_glmm_nb2_zi,
                         #Solanaceae_glmm_hurdle_poisson_zi,
                         Solanaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Solanaceae
write.csv(AIC_Solanaceae, "AIC_Solanaceae.csv")
tab_model(Solanaceae_glmm) ### tab best model
parameters(Solanaceae_glmm)
r2(Solanaceae_glmm)

# Check for overdispesion and zero inflation 
summary(Solanaceae_glmm)
check_overdispersion(Solanaceae_glmm) # No overdispersion detected.
check_zeroinflation(Solanaceae_glmm) # Model seems ok 
plot(allEffects(Solanaceae_glmm))
check_model(Solanaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Solanaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

glmm_graph_Solanaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_point(data = data_family, aes(x = treatment, y = Solanaceae, color = site_letter), position = position_jitterdodge(dodge.width = 0.2), size = 2, alpha = 0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Solanaceae") +
  xlab ("")

glmm_graph_Solanaceae
ggsave(file="Solanaceae.jpg", 
       plot= glmm_graph_Solanaceae,
       width=15,height=15,units="cm",dpi=300)

family_seeds <- ggarrange(glmm_graph_total, 
                          glmm_graph_Moraceae,
                          glmm_graph_Poaceae,
                          glmm_graph_piperaceae,
                          glmm_graph_Solanaceae,
                          glmm_graph_Urticaceae,
                          ncol = 3,
                          nrow = 2,
                          align = "hv",
                          common.legend = TRUE,
                          legend="right")
family_seeds
ggsave(file="Figure10.jpg", 
       plot= family_seeds,
       width=25,height=16,units="cm",dpi=300)

#######################
#######################
#######################

### NMDS ###
head(data_family)
mmatrix <- data_family[, c(4:9)]
mmatrix <- mmatrix + 1   #add a small constant to deal with the excess of zeros. Without the constant the NMDS does not run.
matrix <- as.matrix(mmatrix) # turn data frame into matrix
head(matrix)
nmds_results <- metaMDS(comm = matrix,
                        autotransform = FALSE,
                        distance = "bray",
                        trymax = 100)
nmds_results$stress
plot(nmds_results, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$site <- data_family$site_letter
data.scores$treatment <- data_family$treatment
data.scores

# PERMANOVA
adonis_site <- adonis2(matrix ~ data.scores$site, distance = "bray", perm=999)
adonis_site
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data_family$site_letter
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

color_palette <- viridis(20, option = "G")
selected_colors <- color_palette[c(5, 10, 15, 18)]

nmdsgraph_site <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2, color = site)) +
  geom_point(size = 3) + 
  theme_classic(base_size = 25) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right") + 
  annotate("text", x = -2, y = 2.5, size = 6, label = paste("PERMDISP2, P = 0.654\nPERMANOVA, P = 0.001 ***")) +
  stat_ellipse(level = 0.95, aes(color = site)) 
nmdsgraph_site 

# PERMANOVA
adonis_t <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=999)
adonis_t
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data_family$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)

nmdsgraph_treatment <- ggplot(data.scores, aes(x = MDS1, y = MDS2, color = treatment)) +
  geom_point(size = 3) + 
  theme_classic(base_size = 25) +
  scale_color_manual(values = selected_colors, name = "Treatment", labels = c("control" = "Control", "treatment" = "Lures")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right") + 
  annotate("text", x = -2.5, y = 2, size = 6, label = paste("PERMDISP2, P = 0.170\nPERMANOVA, P = 0.847")) +
  stat_ellipse(level = 0.95)
nmdsgraph_treatment

### Group all community figures together
nmds_seeds_hor <- ggarrange(nmdsgraph_site,
                            nmdsgraph_treatment,
                            ncol = 1,
                            nrow = 2,
                            align = "hv")

seed_community <- ggarrange(nmds_seeds_hor,
                            allseeds_total,
                            ncol = 2,
                            nrow = 1)
seed_community
ggsave(file="seed_community.jpg", 
       plot= seed_community,
       width=55,height=30,units="cm",dpi=300)

##################
### chi_square ###
##################

head(data_family)
colnames(data_family)
filter_data <- data_family[, c("treatment", "Piperaceae", "Urticaceae", "Poaceae", "Moraceae", "Solanaceae", "Unknown")]
filter_data <- filter_data %>%
  group_by(treatment) %>%
  summarise(across(everything(), sum))

# Create a matrix with the counts for 'control' and 'treatment' for each family
counts_matrix <- matrix(c(19, 188, 926, 742, 9, 16, 5, 97, 15, 6, 21, 12), nrow = 2, byrow = TRUE)

# List of family names
family_names <- c("Piperaceae", "Urticaceae", "Poaceae", "Moraceae", "Solanaceae", "Unknown")
family_names
# List to store individual chi-square test results
chi_square_results <- list()

# Perform chi-square tests for each family
for (i in 1:length(family_names)) {
  family_name <- family_names[i]
  family_counts <- counts_matrix[, i]
  
  # Create a contingency table
  contingency_table <- matrix(family_counts, nrow = 2, byrow = TRUE)
  rownames(contingency_table) <- c("control", "treatment")
  
  # Perform chi-square test
  chi_square_test_result <- chisq.test(contingency_table, simulate.p.value = TRUE)
  
  # Store the result in the list
  chi_square_results[[family_name]] <- chi_square_test_result
}

# Print the results for each family
for (family_name in family_names) {
  cat("Family:", family_name, "\n")
  print(chi_square_results[[family_name]])
  cat("\n")
}

#########################
### diversity indices ###
#########################

head(data_cleaned)
# Select columns 4 to 25 and calculate species richness
sp_matrix <- data_cleaned[, 4:25]
data_cleaned$morpho_number <- specnumber(sp_matrix)
# Calculate Shannon and Simpson diversity indices
data_cleaned$shannon <- diversity(sp_matrix)
data_cleaned$simpson <- diversity(sp_matrix, "simpson")

# Distribution? create histograms for species richness, Shannon, and Simpson indices
hist(data_cleaned$morpho_number)
hist(data_cleaned$shannon)
hist(data_cleaned$simpson)

# Filter out rows where 'treatment' is not equal to basal
data_cleaned <- data_cleaned %>%
  filter(!treatment %in% "basal")

# Analyze species richness using a gamma GLMM model
data_cleaned$morpho_number_gamma <- data_cleaned$morpho_number + 1
morpho_number_glmm <- glmmTMB(morpho_number_gamma ~ treatment + (1|site_letter) + (1|week), data = data_cleaned, family=Gamma(link = "log"))
summary(morpho_number_glmm)
tab_model(morpho_number_glmm)
r2(morpho_number_glmm)
plot(allEffects(morpho_number_glmm))

# Create a violin plot for morpho species richness by treatment
morpho_number <- ggplot(data_cleaned, aes(x = treatment, y = morpho_number)) +
  theme_classic(base_size = 25) + 
  geom_violin() +
  geom_point() 
morpho_number

# Analyze Shannon diversity using a gamma GLMM model
data_cleaned$shannon_gamma <- data_cleaned$shannon + 1
shannon_glmm <- glmmTMB(shannon_gamma ~ treatment + (1|site_letter) + (1|week), data = data_cleaned, family=Gamma(link = "log"))
summary(shannon_glmm)
tab_model(shannon_glmm)
r2(shannon_glmm)
plot(allEffects(shannon_glmm))

# Create a violin plot for Shannon diversity by treatment
shannon <- ggplot(data_cleaned, aes(x = treatment, y = shannon)) +
  theme_classic(base_size = 25) + 
  geom_violin() +
  geom_point() 
shannon

# Analyze Simpson diversity using a gamma GLMM model
data_cleaned$simpson_gamma <- data_cleaned$simpson + 1
simpson_glmm <- glmmTMB(simpson ~ treatment + (1|site_letter) + (1|week), data = data_cleaned)
summary(simpson_glmm)
tab_model(simpson_glmm)
r2(simpson_glmm)
plot(allEffects(simpson_glmm))

# Create a violin plot for Simpson diversity by treatment
simpson <- ggplot(data_cleaned, aes(x = treatment, y = simpson)) +
  theme_classic(base_size = 25) + 
  geom_violin() +
  geom_point() 
simpson

### No significant differences in the diversity indices calculated here 

#######################################
#######################################
#######################################
#######################################



