# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, glmmTMB, performance, DHARMa, emmeans, ggeffects, car) # Install missing packages and load needed libraries


# loading data
jaccard_patch <- read_csv(file = file.path("data", "L2_summarized", "patch_jaccard.csv"))

jaccard_patch <- jaccard_patch %>%
  filter(patch_type != "center") %>% # removing center patch from analysis
  mutate(patch_type = as.factor(patch_type))

#### Q1: How is connectivity related to community stability across succession? ####
m1 <- glmmTMB(COMMON_jaccard ~ patch_type * time2 + (1|block/patch),
              data = jaccard_patch)
summary(m1)
plot(simulateResiduals(m1))
Anova(m1)

# posthoc for interaction
m1_posthoc <- emtrends(m1, specs = pairwise ~ patch_type, var = "time2")
m1_posthoc

# posthoc for patch type
m1_posthoc <- emmeans(m1, specs = pairwise ~ patch_type * time2)
pairs(m1_posthoc, simple = "patch_type")

# plotting
m1_predict <- ggpredict(m1, terms = c("time2", "patch_type"))
m1_plot <- m1_predict %>%
  rename(patch_type = group) %>%
  ggplot() +
  geom_line(aes(x, predicted, color = patch_type), linewidth = 2) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = patch_type), alpha = 0.3) + 
  geom_point(aes(time2, jaccard_dissimilarity, color = patch_type), data = jaccard_patch, size = 2.5, alpha = 0.3) +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  xlab("Year of succession") +
  ylab("Jaccard dissimilarity") +
  labs(color = "Patch type", fill = "Patch type") +
  theme(panel.spacing = unit(1.2, "lines"))

pdf(file = "jaccard.pdf", width = 10, height = 6)
m1_plot
dev.off()
#jaccard_plot <- jaccard_patch %>% # plotting jaccard dissimilarity across year 
#  ggplot(aes(time2, jaccard_dissimilarity, color = patch_type)) + 
#  facet_wrap(~patch_type) + 
#  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_ID), linewidth = 1, alpha = 0.3) +
#  theme_bw() + 
#  scale_color_brewer(palette = "Set2") +
 # stat_smooth(method = "lm", se = F, linewidth = 3)
#jaccard_plot


# How does temporal gamma diversity play into this? 
# Would expect that a higher temporal beta diversity would increase gamma diversity



# initial turnover did not differ between patch types 
initial_jaccard <- jaccard_patch %>%
  filter(time2 == 2)

m2 <- glmmTMB(jaccard_dissimilarity ~ patch_type + (1|block/patch),
              data = initial_jaccard)
summary(m2)



#### dispersal mode ####
# gravity dispersed 
m.gravity <- glmmTMB(GRAVITY_jaccard ~ patch_type * time2 + (1|block/patch),
                     data = jaccard_patch)
summary(m.gravity)
# posthoc for interaction
m.gravity_posthoc <- emtrends(m.gravity, specs = pairwise ~ patch_type, var = "time2")
m.gravity_posthoc

# posthoc for patch type
m.gravity_posthoc <- emmeans(m.gravity, specs = pairwise ~ patch_type * time2)
pairs(m.gravity_posthoc, simple = "patch_type")


# wind dispersed 
m.wind <- glmmTMB(WIND_jaccard ~ patch_type * time2 + (1|block/patch),
                     data = jaccard_patch)
summary(m.wind)
# posthoc for interaction
m.wind_posthoc <- emtrends(m.wind, specs = pairwise ~ patch_type, var = "time2")
m.wind_posthoc

# posthoc for patch type
m.wind_posthoc <- emmeans(m.wind, specs = pairwise ~ patch_type * time2)
pairs(m.wind_posthoc, simple = "patch_type")




# animal dispersed 
m.animal <- glmmTMB(ANIMAL_jaccard ~ patch_type * time2 + (1|block/patch),
                  data = jaccard_patch)
summary(m.animal)
# posthoc for interaction
m.animal_posthoc <- emtrends(m.animal, specs = pairwise ~ patch_type, var = "time2")
m.animal_posthoc

# posthoc for patch type
m.animal_posthoc <- emmeans(m.animal, specs = pairwise ~ patch_type * time2)
pairs(m.animal_posthoc, simple = "patch_type")





#### regressing beta diversity against mean alpha diversity to control for alpha diversity ####
m_beta_corrected <- glmmTMB(COMMON_jaccard ~ COMMON_alpha,
                            data = jaccard_patch)
summary(m_beta_corrected)
resid <- residuals(m_beta_corrected, type = "response")
jaccard_patch$resid_beta_corrected <- resid


m_beta <- glmmTMB(resid_beta_corrected ~ patch_type + time2 + (1|block/patch),
                  data = jaccard_patch)
summary(m_beta)
Anova(m_beta)

# posthoc for interaction
m_beta_posthoc <- emtrends(m_beta, specs = pairwise ~ patch_type, var = "time2")
m_beta_posthoc

# posthoc for patch type
m_beta_posthoc <- emmeans(m_beta, specs = pairwise ~ patch_type)
pairs(m_beta_posthoc, simple = "patch_type")

# plotting
# plotting
m_resid_predict <- ggpredict(m_beta, terms = c("time2", "patch_type"))
m_resid_plot <- m_resid_predict %>%
  rename(patch_type = group) %>%
  ggplot() +
  geom_line(aes(x, predicted, color = patch_type), linewidth = 2) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = patch_type), alpha = 0.3) + 
  geom_point(aes(time2, resid_beta_corrected, color = patch_type), data = jaccard_patch, size = 2.5, alpha = 0.3) +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  xlab("Year of succession") +
  ylab("Jaccard dissimilarity (residuals)") +
  labs(color = "Patch type", fill = "Patch type") +
  theme(panel.spacing = unit(1.2, "lines"))
m_resid_plot

pdf(file = "jaccard_resid.pdf", width = 10, height = 6)
m_resid_plot
dev.off()

jaccard_results %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time2, resid_beta_corrected, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_id), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)

