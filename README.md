#  Testing the effectiveness of synthetic chemical lures to increase fruit bat-mediated seed dispersal in a tropical forest

This repository contains scripts, information, and figures related to the statistical analyses conducted for the research paper titled "Testing the effectiveness of synthetic chemical lures to increase fruit bat-mediated seed dispersal in a tropical forest."
Our research explores the effectiveness of synthetic chemical lures as a novel strategy to attract fruit bats and enhance seed rain in a lowland rainforest in northeastern Costa Rica (La Selva Biological Station). We investigated the impact of chemical lures on increasing bat activity and seed rain. All the analyses were performed in R v. 4.2.1.

## Scripts

### 1. Objective 1:  Assess the effectiveness of the chemical lure to increase bat activity in open and semi-open areas (`script1.R` and `script2.R`)

These scripts details the process of analyzing the impact of chemical lures on bat activity. `script1.R` compares bat communities across different sites and treatments. Non-metric multidimensional scaling (NMDS) was employed for visualization. Homogeneity of variances was checked using the 'betadisper()' function, followed by permutational multivariate analysis of variance (PERMANOVA) using the 'adonis2()' function with 999 permutations.

In `script2.R` GLMMs were employed using the glmmTMB package. The models included the bat abundance as fixed effect, and site and capture date as random effects. The analysis was performed using various count data distributions from the glmmTMB package Overdispersion and zero inflation were assessed using the 'check_overdispersion()' and 'check_zeroinflation()' functions from the performance package. Effect sizes were computed based on estimated marginal means using the 'emmeans()' function from the emmeans package. We performed an autocorrelation analysis on the residuals of each model fitted. First, we performed a Durbin-Watson (DW) using the function 'dwtest()' from the lmtest package to test to assess temporal autocorrelation in these residuals. Then, we generated a visual representation of the autocorrelation function (ACF). Our results indicate that there is no temporal autocorrelation present in our bat data.

### 1. Objective 2: Assess the effectiveness of the chemical lure to increase seed rain of open and semi-open areas (`script3.R`)

This script outlines the analysis of the impact of chemical lures on seed rain NMDS was used for visualization, and homogeneity of variances was checked with 'betadisper()'. PERMANOVA was conducted using the 'adonis2()' function with 999 permutations to test for statistical significance.

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