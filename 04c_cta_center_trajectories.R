### community trajectory analysis - segement lengths
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")
