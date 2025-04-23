# SRS community stability

Contact: Katherine Hulting, hultingk@msu.edu

### Scripts
#### 01_wrangling.R
Script for data cleaning and wrangling. Output is cleaned plant occurence per year. 

#### 02_local-stability.R
Script for calculating temporal beta diversity for each local patch. Requires 01_wrangling.R to be run. Output is summarized at one value per patch per year. 

#### 03_analysis.R
Script for analysis. Requires 02_local-stability.R to be run.
