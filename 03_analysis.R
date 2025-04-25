# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, glmmTMB, performance, DHARMa, emmeans, ggeffects) # Install missing packages and load needed libraries


# loading data
jaccard_patch <- read_csv(file = file.path("data", "L2_summarized", "patch_jaccard.csv"))

jaccard_patch <- jaccard_patch %>%
  filter(patch_type != "center") %>% # removing center patch from analysis
  mutate(patch_type = as.factor(patch_type))

#### Q1: How is connectivity related to community stability across succession? ####
# need to add alpha diversity into this to check

m1 <- glmmTMB(jaccard_dissimilarity ~ patch_type * time2 + (1|block/patch),
              data = jaccard_patch)
summary(m1)
plot(simulateResiduals(m1))

# posthoc for interaction
m1_posthoc <- emtrends(m1, specs = pairwise ~ patch_type, var = "time2")
m1_posthoc

# posthoc for patch type
m1_posthoc <- emmeans(m1, specs = pairwise ~ patch_type * time2)
pairs(m1_posthoc, simple = "patch_type")


m1_predict <- ggpredict(m1, terms = c("time2", "patch_type"))
m1_predict %>%
  ggplot() +
  geom_line(aes(x, predicted, color = group), linewidth = 2) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) + 
  geom_point(aes(time2, jaccard_dissimilarity, color = patch_type), data = jaccard_patch) +
  theme_classic() + 
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")
  



jaccard_plot <- jaccard_patch %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time2, jaccard_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_ID), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)
jaccard_plot




