#  Synthetic chemical lures increase *Carollia* spp. (Chiroptera: Phyllostomidae) activity and seed rain in a Neotropical forest 

This repository contains scripts, information, and figures related to the statistical analyses conducted for the research paper titled "Synthetic chemical lures increase Carollia spp. (Chiroptera: Phyllostomidae) activity and seed rain in a Neotropical forest" The analyses were performed in R v. 4.2.1. 

## Scripts

### 1. Objective 1:  Assess the effectiveness of the chemical lure to increase bat activity in open and semi-open areas (`script1.R` and `script2.R`)

These scripts details the process of analyzing the impact of chemical lures on bat activity. `script1.R` presents the process of comparing bat communities across different sites and treatments. Non-metric multidimensional scaling (NMDS) was employed for visualization. Homogeneity of variances was checked using the 'betadisper()' function, followed by permutational multivariate analysis of variance (PERMANOVA) using the 'adonis2()' function with 999 permutations.

In `script2.R` GLMMs were employed using the glmmTMB package. The models included the interaction between bat abundance and bat species as fixed effects, with site and capture date as random effects. The analysis was performed using various count data distributions from the glmmTMB package. Model comparison was done using the Akaike Information Criterion (AIC). The 'Anova()' function from the car package was used to evaluate the effects and interactions of each factor. Overdispersion and zero inflation were assessed using the 'check_overdispersion()' and 'check_zeroinflation()' functions from the performance package. Effect sizes were computed based on estimated marginal means using the 'emmeans()' function from the emmeans package.

### 1. Objective 2: Assess the effectiveness of the chemical lure to increase seed rain of open and semi-open areas (`script3.R`)

This script outlines the analysis of the impact of chemical lures on seed rain. GLMMs were constructed using the glmmTMB package, considering the interaction between seed count and plant family as fixed effects, with site and collection date as random effects. Similar to Objective 1, count data distributions were used, and AIC was employed for model comparison. The 'Anova()' function from the car package assessed the effects and interactions of factors. Overdispersion and zero inflation were checked using the 'check_overdispersion()' and 'check_zeroinflation()' functions. Effect sizes were calculated using the 'emmeans()' function. Then, presents the comparison of seed communities across sites and treatments. NMDS was used for visualization, and homogeneity of variances was checked with 'betadisper()'. PERMANOVA was conducted using the 'adonis2()' function with 999 permutations to test for statistical significance.

## Data Files

### Folder Objective 1: `data.csv` and `data_nodates.csv`

Contains the data used for analyzing bat activity. The columns in the dataset are as follows:

- `date`: The date of bat capture.
- `site`: The site where the capture took place.
- `bats`: Total number of bats captured.
- `fruit_bats`: Total number of captured fruit bats.
- `cperspicillata`: Number of captured *Carollia perspicillata* bats.
- `csowelli`: Number of captured *Carollia sowelli* bats.
- `ccastanea`: Number of captured *Carollia castanea* bats.
- `carollia_spp`: Total number of captured bats from the *Carollia* genus.
- `uroderma_spp`: Number of captured bats from the *Uroderma* genus.
- `sturnira_spp`: Number of captured bats from the *Sturnira* genus.
- `ectophylla_alba`: Number of captured *Ectophylla alba* bats.
- `artibeus_spp`: Number of captured bats from the *Artibeus* genus.
- `desmodus_rotundus`: Number of captured *Desmodus rotundus* bats.
- `nectarivorous_bats`: Number of captured nectarivorous bats.
- `insectivorous_bats`: Number of captured insectivorous bats.
- `treatment`: The treatment applied ("control" or "lures").
- `hours`: Total hours of mist nesting. 
- `nets`: Total number of nets used for bat capture. 

### Folder Objective 2: `seed_data.csv`

Contains the data used for analyzing seed rain. The columns in the dataset are as follows:

- `week`: The week of the seed collection. Seeds were collected every 15 days.
- `site_name`: The name of the site where the observation took place at La Selva ("Zompopa", "Lab", "STR - Sendero Tres Rios", "PS - Parcelas de Sucesion").
- `site_letter`: The site's letter designation ("A", "B", "C", "D")
- `collection_date`: The date of seed collection.
- `treatment`: The treatment applied ("baseline", "control", "treatment").
- Columns for various plant species/plant families, indicating the count of seeds for each species.
- `total`: The total count of seeds for all plant species per collection week.
- `comments`: Additional comments or notes about the observation.

## Figures folder

The 'Figures' folder contains various output files generated from the analyses conducted in the main scripts. These figures visually represent the results and insights obtained from the data exploration and statistical modeling.