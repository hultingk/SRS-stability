# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, glmmTMB, performance, DHARMa, emmeans) # Install missing packages and load needed libraries


# loading data
jaccard_patch <- read_csv(file = file.path("data", "L2_summarized", "patch_jaccard.csv"))



#### Q1: How is connectivity related to community stability across succession? ####
## Once spp lists are finalized: set up glmms, patch type*time, exclude center patch for this analysis, 
### block and patch nested random intercepts, time between survey as random intercept?
### interactive model, drop interaction if not significant (emtrends)











