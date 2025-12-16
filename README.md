# SRS plant community assembly

Contact: Katherine Hulting, hultingk@msu.edu

### Scripts
#### 00_functions.R
Script containing custom functions used throughout analysis.

#### 01_wrangling.R
Script for data cleaning and wrangling.
* Requires: access to "LO_original" folder containing data (not committed). 
* Products: srs_plant_all.csv (cleaned plant community data)

#### 02_pcoa_permanova.R
Script for visualizing PCoA of plant community trajectories 
* Requires: 01_wrangling.R 
* Products: figureS1.pdf 

#### 03_convergence.R
Script for calculating convergence/divergence between patch type communities across time, repeat within dispersal mode groups
* Requires: 00_functions.R, 02_pcoa_permanova.R
* Products: figure2.pdf, tableS1.html, tableS2.html, tableS3.html

#### 03b_colonization_extinction.R
Script for partitioning changes in spatial beta diversity over time into changes due to colonization and extinction
* Requires: 00_functions.R, 01_wrangling.R
* Products: figureS5.pdf

#### 04_cta_segments.R
Script for calculating interannual trajectory distances, repeat within dispersal mode groups
* Requires: 00_functions.R, 01_wrangling.R
* Products: figure3.pdf, tableS4.html, tableS5.html, table S6.html

#### 05_cta_directionality.R
Script for calculating trajectory directionality, repeat within dispersal mode groups
* Requires: 04_cta_segments.R
* Products: tableS7.html, tableS8.html, figure4.pdf

#### 06_supplement.R
Script for supplemental exploratory plots: species richness over time, proportion of dispersal mode, turnover
* Requires: 00_functions.R, 01_wrangling.R
* Products: figureS2.pdf, figureS3.pdf, figureS6.pdf



