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
library(AICcmodavg)
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
supp.labs <- c("Control", "Chemical\nlures", "Baseline")
names(supp.labs) <- c("control", "treatment", "basal")

allseeds_total <- ggplot(filtered_df, aes(x = site_letter , y = seed_count, fill = plant_species)) +
  theme_test(base_size = 15) +
  ggtitle("c")+ 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  theme(strip.background = element_rect(colour = "NA", fill = "NA")) +
  scale_fill_viridis_d(option = "inferno", direction = -1, name = "Seed ID",
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
  ylab("Total count of seeds") +
  xlab("Sites") +
  guides(fill = guide_legend(ncol = 1)) + 
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

allseeds_total

### merge unknowns 

# Combine all unknowns into a single category in filtered_df
filtered_df <- filtered_df %>%
  mutate(plant_species = ifelse(grepl("unknown", plant_species), "Unknown", plant_species))

# Update the labels for the merged unknown category
updated_labels <- c(
  "Unknown" = "Unknown",
  "cecropia_insignis" = parse(text = "italic('Cecropia insignis')"),
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
  "achenio_sp" = parse(text = "italic('Achenio') ~'sp.'"),
  "solanaceae" = "Solanaceae",
  "ficus_insipida" = parse(text = "italic('Ficus insipida')")
)

# Create the graph with merged unknowns
allseeds_total <- ggplot(filtered_df, aes(x = site_letter, y = seed_count, fill = plant_species)) +
  theme_test(base_size = 15) +
  ggtitle("c") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10)) +
  theme(strip.background = element_rect(colour = "NA", fill = "NA")) +
  scale_fill_viridis_d(option = "inferno", direction = -1, name = "Seed ID", labels = updated_labels) +
  geom_bar(stat = "identity") +
  ylab("Total count of seeds") +
  xlab("Sites") +
  guides(fill = guide_legend(ncol = 1)) + 
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

allseeds_total

allseeds_total_log <- ggplot(filtered_df, aes(x = site_letter, y = seed_count, fill = plant_species)) +
  theme_test(base_size = 15) +
  ggtitle("c") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10)) +
  theme(strip.background = element_rect(colour = "NA", fill = "NA")) +
  scale_fill_viridis_d(option = "inferno", direction = -1, name = "Seed ID", labels = updated_labels) +
  geom_bar(stat = "identity") +
  ylab("Total seed number (log scale)") +
  xlab("Sites") +
  scale_y_log10() +  # Log-transform the y-axis
  guides(fill = guide_legend(ncol = 1)) + 
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

allseeds_total_log

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
write.csv(data_family, file = "data_family.csv")

summary_data <- data_family %>%
  group_by(treatment, week) %>%
  summarize(
    mean_Piperaceae = round(mean(Piperaceae), 2),
    sd_Piperaceae = round(sd(Piperaceae), 2),
    mean_Urticaceae = round(mean(Urticaceae), 2),
    sd_Urticaceae = round(sd(Urticaceae), 2),
    mean_Poaceae = round(mean(Poaceae), 2),
    sd_Poaceae = round(sd(Poaceae), 2),
    mean_Moraceae = round(mean(Moraceae), 2),
    sd_Moraceae = round(sd(Moraceae), 2),
    mean_Solanaceae = round(mean(Solanaceae), 2),
    sd_Solanaceae = round(sd(Solanaceae), 2),
    mean_Unknown = round(mean(Unknown), 2),
    sd_Unknown = round(sd(Unknown), 2),
    mean_total = round(mean(total), 2),
    sd_total = round(sd(total), 2)
  ) %>%
  mutate(
    Piperaceae = paste0(mean_Piperaceae, " ± ", sd_Piperaceae),
    Urticaceae = paste0(mean_Urticaceae, " ± ", sd_Urticaceae),
    Poaceae = paste0(mean_Poaceae, " ± ", sd_Poaceae),
    Moraceae = paste0(mean_Moraceae, " ± ", sd_Moraceae),
    Solanaceae = paste0(mean_Solanaceae, " ± ", sd_Solanaceae),
    Unknown = paste0(mean_Unknown, " ± ", sd_Unknown),
    total = paste0(mean_total, " ± ", sd_total)
  )

summary_data_2 <- summary_data[, c(1,16,17,18,19,20,21,22)]

summary_data_long <- summary_data_2 %>%
  pivot_longer(
    cols = -treatment,
    names_to = "family",
    values_to = "mean_sd"
  )
summary_data_long

summary_data_wide <- summary_data_long %>%
  pivot_wider(
    names_from = treatment,
    values_from = mean_sd
  )

print(summary_data_wide)
write.csv(summary_data_wide, file = "summary_data_wide.csv")

### just total 

summary_data <- data_family %>%
  group_by(treatment, week) %>%
  summarize(
    mean_total = round(mean(total), 2),
    sd_total = round(sd(total), 2)
  ) %>%
  mutate(total = paste0(mean_total, " ± ", sd_total)
  )

summary_data_2 <- summary_data[, c(1, 2, 5)]

summary_data_wide_total <- summary_data_2%>%
  pivot_wider(
    names_from = treatment,
    values_from = total
  )
write.csv(summary_data_wide_total, file = "summary_data_wide_total.csv")

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
summary(total_glmm)
check_overdispersion(total_glmm) # Overdispersion detected.
check_zeroinflation(total_glmm) # Model is underfitting zeros (probable zero-inflation).
tab_model(total_glmm)

total_glmm_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
summary(total_glmm_zi)
check_overdispersion(total_glmm_zi) # Overdispersion detected.
check_zeroinflation(total_glmm_zi) # Model is underfitting zeros (probable zero-inflation).
tab_model(total_glmm_zi)


glmm_emmeans <-emmeans(total_glmm_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# control 31.38, lures 72.64
72.64 - 31.38 # 41.26
(41.26*100)/31.38 # 131.485
72.64/31.38

dwtest(total_glmm_zi) # Durbin-Watson test
# First, store the residuals
residuals <- resid(total_glmm_zi)
# Calculate ACF
acf_data <- acf(residuals, plot = FALSE)
# Convert ACF data to dataframe
acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
# Calculate confidence intervals
ci_upper <- qnorm(0.975) / sqrt(length(residuals))
ci_lower <- -qnorm(0.975) / sqrt(length(residuals))

ACF_seeds <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_bar(stat = "identity", width = 0.1) +
  geom_hline(yintercept = c(ci_upper, ci_lower, 0), linetype = c("dashed", "dashed", "solid"), color = c("darkgrey", "darkgrey", "black")) +
  labs(title = "ACF Total seed count",
       x = "Lag", y = "Autocorrelation")+ 
  theme_test(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  geom_text(x = 9, y = 0.8, label = "DW = 2.477,\np-value = 0.798", size = 3)

ACF_seeds

temporal_variation_seeds <- ggplot(data_family, aes(x = week, y = total, color = site_letter, shape = treatment)) +
  theme_test(base_size = 10) +
  geom_point(size = 2.5, color = "darkgrey") +  
  geom_point(size = 2) + 
  scale_shape_manual(values = c("control" = "circle", "treatment" = "triangle"),
                     name = "Treatment", labels = c("Control", "Chemical\nlures")) + 
  scale_color_viridis_d(option = "D", name = "Sites") +
  xlab("Collection week") +
  ylab ("Total seed count")

temporal_variation_seeds

final_temporal_ACF_seeds <- ggarrange(temporal_variation_seeds,
                                      ACF_seeds,
                                      ncol = 2,
                                      nrow = 1,
                                      align = "hv")

final_temporal_ACF_seeds

ggsave(file="final_temporal_ACF_seeds.jpg", 
       plot= final_temporal_ACF_seeds,
       width=6,height=3,units="in",dpi=500)

head(data_family)
color_palette <- viridis(20, option = "A")
selected_colors <- color_palette[c(7, 13, 16, 19)]

data_family$totaladd2 <- data_family$total + 2

glmm_graph_total <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_test(base_size = 15) + 
  geom_violin(data = data_family, color = "white", fill = "light grey", alpha = 0.3, aes(x = treatment, y = totaladd2)) +
  geom_point(data = data_family, aes(x = treatment, y = totaladd2, color = site_letter, size = species_richness), position = position_jitterdodge(dodge.width = 0.1)) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 5) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Total count of seeds") +
  xlab ("") +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  geom_text(x = 1.5, y = 2.8, label = "***", size = 8) + 
  scale_size(name = "Family\nrichness") + 
  scale_y_log10(breaks = c(0, 10, 50, 100, 200, 400, 800), limits = c(2, 1000))
glmm_graph_total
ggsave(file="totalseed_glmm.jpg", 
       plot= glmm_graph_total,
       width=10,height=10,units="cm",dpi=600)

#######################
#######################
#######################

### NMDS ###
head(data_family)
data_family <- data_family %>% filter(total != 0)
data_family

mmatrix <- data_family[, c(4:9)]
mmatrix
mmatrix_filtered <- log(mmatrix + 1)
dist_matrix <- vegdist(mmatrix_filtered, method = "bray")
nmds_results <- metaMDS(comm = dist_matrix,
                        autotransform = FALSE,
                        distance = "bray",
                        trymax = 300)
nmds_results$stress
plot(nmds_results, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$site <- data_family$site_letter
data.scores$treatment <- data_family$treatment
data.scores

# PERMANOVA
adonis_site <- adonis2(dist_matrix ~ data.scores$site, distance = "bray", perm=999)
adonis_site
# BETADISP
dist_matrix <- vegdist(dist_matrix, method = "bray")
groups <- data_family$site_letter
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

color_palette <- viridis(20, option = "G")
selected_colors <- color_palette[c(5, 10, 15, 18)]

nmdsgraph_site <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2, color = site)) +
  theme_test(base_size = 15) +
  ggtitle("a")+ 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  geom_point(size = 2) + 
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  annotate("text", x = -1.5, y = 1.8, size = 3, label = paste("PERMDISP2, P = 0.029*\nPERMANOVA, P = 0.001***"), hjust = 0) +
  ylim(-2, 2) +
  xlim(-3, 3)
nmdsgraph_site 

?adonis2
# PERMANOVA
adonis_t <- adonis2(dist_matrix ~ data.scores$treatment, strata = data_family$site_letter, distance = "bray", perm=999)
adonis_t
# BETADISP
groups <- data_family$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)

nmdsgraph_treatment <- ggplot(data.scores, aes(x = MDS1, y = MDS2, color = treatment)) +
  theme_test(base_size = 15) +
  ggtitle("b")+ 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(size = 2) + 
  scale_color_viridis_d(option="D", name= "Treatment", labels = c("control" = "Control", "treatment" = "Chemical\nlures")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10)) + 
  annotate("text", x = -1, y = 1.8, size = 3, label = paste("PERMDISP2, P = 0.553\nPERMANOVA, P = 0.547"),  hjust = 0) +
  stat_ellipse(level = 0.95) + 
  ylim(-2, 2) +
  xlim(-3, 3)
nmdsgraph_treatment

### Group all community figures together
nmds_seeds_hor <- ggarrange(nmdsgraph_site,
                            nmdsgraph_treatment,
                            ncol = 1,
                            nrow = 2,
                            align = "hv")

seed_community <- ggarrange(nmds_seeds_hor,
                            allseeds_total,
                            widths = c(2, 3),
                            ncol = 2,
                            nrow = 1)
seed_community
ggsave(file="seed_community.jpg", 
       plot= seed_community,
       width=9.5,height=5.5,units="in",dpi=600)

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


### For reviwers, removing the outliers
data_family <- data_family %>%
  filter(total <= 400)

total_glmm <- glmmTMB(total ~ treatment, data = data_family, family = poisson)
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
glmm_emmeans <-emmeans(total_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans


head(data_family)
color_palette <- viridis(20, option = "A")
selected_colors <- color_palette[c(7, 13, 16, 19)]

glmm_graph_total <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_test(base_size = 15) + 
  geom_violin(data = data_family, color = "white", fill = "light grey", alpha = 0.3, aes(x = treatment, y = total)) +
  geom_point(data = data_family, aes(x = treatment, y = total, color = site_letter), size = 5,position = position_jitterdodge(dodge.width = 0.2)) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 5) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical\nlures")) + 
  ylab ("Total seed count") +
  xlab ("") +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))+
  geom_text(x = 1.5, y = 2.5, label = "***", size = 8) + 
  scale_size(name = "Family\nrichness") +
  scale_y_log10(breaks = c(0, 10, 50, 100, 200, 400), limits = c(5, 500))
glmm_graph_total
ggsave(file="totalseed_glmm_NOOUTLIERS.jpg", 
       plot= glmm_graph_total,
       width=10,height=10,units="cm",dpi=600)

#######################################
#######################################
#######################################
#######################################

