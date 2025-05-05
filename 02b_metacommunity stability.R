# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, vegan) # Install missing packages and load needed libraries

source("00_functions.R")
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)


srs_data_sub <- srs_data %>%
  filter(patch %in% c("A", "B", "C", "D")) %>%
  filter(!year %in% c("2004", "2013", "2024")) # years that center patch was not sampled in any block


srs_metacommunity <- srs_data_sub %>%
  mutate(metacommunity = if_else(str_detect(unique_id, "A-center|B-connected"), 
                                 "connected-metacommunity", "unconnected-metacommunity")) %>%
  mutate(unique_id = paste(block, metacommunity, sep = "-"))


srs_metacommunity_split <- srs_metacommunity %>%
  dplyr::count(unique_id, time, sppcode) %>%
  mutate(n = 1) %>%
  group_by(unique_id) %>%
  group_split()

metacommunity_jaccard <- srs_metacommunity_split %>%
  lapply(compute_jaccard) %>%
  bind_rows()

metacommunity_jaccard <- metacommunity_jaccard %>%
  separate(unique_id, into = c("block", "patch_type", "metacommunity"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) %>% # converting to numeric
  mutate(years_bw_surveys = time2 - time1)


metacommunity_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time2, jaccard_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_point(size = 2, alpha = 0.7) + 
 #geom_line(aes(time, jaccard_dissimilarity, group = unique_id), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
 stat_smooth(method = "lm", se = F, linewidth = 3)

#### regressing against alpha diversity ####
m_resid <- glmmTMB(jaccard_dissimilarity ~ mean_alpha_diversity,
            data = metacommunity_jaccard)
summary(m_resid)
resid <- residuals(m_resid, type = "response")
metacommunity_jaccard$resid_beta_corrected <- resid


m_metacommunity_beta <- glmmTMB(resid_beta_corrected ~ patch_type*time2 + (1|block),
                  data = metacommunity_jaccard)
summary(m_metacommunity_beta)
plot(simulateResiduals(m_metacommunity_beta))
