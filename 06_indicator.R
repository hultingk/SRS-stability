librarian::shelf(tidyverse, vegan, indicspecies, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

source(here::here("02_pcoa_permanova.R"))

patch_info <- patch_info %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"))


indicator <- multipatt(sp_info, cluster = patch_info$patch_type,
                       func = "IndVal.g",
                       control = how(nperm = 999))

summary(indicator)
coverage(sp_info, indicator, At = 0.8, alpha = 0.05)
