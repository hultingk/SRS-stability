### community trajectory analysis - segement lengths
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg, performance) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode, soil_moisture, year_since_fire) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.numeric(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year, soil_moisture, year_since_fire)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))


# Jaccard distance matrix
jaccard_dist <- vegdist(sp_info, method = "jaccard")

# defining trajectories
srs_trajectory <- defineTrajectories(jaccard_dist, sites = patch_info$unique_id, surveys = patch_info$time)


# segment lengths of trajectories between consectutive years
segment_lengths <- trajectoryLengths(srs_trajectory)
segment_lengths <- segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  ))

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys <- patch_info %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths <- cbind(segment_lengths, time_surveys)
segment_lengths$s.time <- as.numeric(scale(segment_lengths$time)) # scaling time


  
# segment lengths model 
# linear
m_length <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                    data = segment_lengths)
summary(m_length)
# quadratic
m_length_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                    data = segment_lengths)
# null
m_length_null <- glmmTMB(distance ~ 1 + (1|block/patch), # null model
                         data = segment_lengths)
# AIC comparison
a <- list(m_length, m_length_quad, m_length_null)
aictab(a) # quadratic much better fit

summary(m_length_quad)
plot(simulateResiduals(m_length_quad))
check_model(m_length_quad)
performance::r2(m_length_quad)



# percent change in segment length from year 2-22 (20 years)
# time 2 = -1.52128914
# time 22 = 1.56397837
# intercept + time_estimate (time) + time^2_estimate (time)
# at time 2
0.242909 + -0.016521*(-1.52128914) + 0.013890*((-1.52128914)^2) # 0.3001881
# at time 22
0.242909 + -0.016521*(1.56397837) + 0.013890*((1.56397837)^2) # 0.2510458

# percent change = time 22 - time 2 / time 2 * 100
(0.2510458 - 0.3001881)/0.3001881 * 100 #-16.3705 % decrease

# 95% CI, percent change
confint(m_length_quad)
# lower 95% CI
# at time 1
0.216317421 + -0.024855072*(-1.52128914) + 0.005869060*((-1.52128914)^2) # 0.2677121
# at time 21
0.216317421 + -0.024855072*(1.56397837) + 0.005869060*((1.56397837)^2) # 0.1918005
# percent change = time 21 - time 1 / time 1 * 100
(0.1918005 - 0.2677121)/0.2677121 * 100 #-28.35569 % 

# upper 95% CI
# at time 1
0.269499648 + -0.008186654*(-1.52128914) + 0.021910606*((-1.52128914)^2) # 0.3326621
# at time 21
0.269499648 + -0.008186654*(1.56397837) + 0.021910606*((1.56397837)^2) # 0.3102899
# percent change = time 21 - time 1 / time 1 * 100
(0.3102899 - 0.3326621)/0.3326621 * 100 #-6.725203 % 




# posthoc
m_length_posthoc <- emmeans(m_length_quad, ~ patch_type*s.time + patch_type * I(s.time^2))
pairs(m_length_posthoc, simple = "patch_type")
m_length_posthoc
m_length_posthoc <- emtrends(m_length_quad, "patch_type", var = "s.time")
pairs(m_length_posthoc)

# prediction plot
# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

m_length_predict <- ggpredict(m_length_quad, terms = c("s.time [all]", "patch_type"))
m_length_predict$group <- factor(m_length_predict$group, levels = c("Connected", "Rectangular", "Winged"))

segment_lengths_plot <- m_length_predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 5.5, alpha = 0.19, data = segment_lengths) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme_minimal(base_size = 24) +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") +
  annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)
segment_lengths_plot


# pdf(file = file.path("plots", "segment_lengths.pdf"), width = 11, height = 8)
# segment_lengths_plot
# dev.off()


# # # plotting segment lengths 
# segment_lengths_plot <- segment_lengths %>%
#   filter(!block %in% c("75W", "75E")) %>%
#   ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
#   geom_point(size = 4.5, alpha = 0.3) +
#   theme_minimal(base_size = 28) +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
# 
#   #scale_color_brewer(palette = "Set2", name = "Patch Type") +
#  # scale_fill_brewer(palette = "Set2", name = "Patch Type") +
#   ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
#   #ylab("Trajectory distance between consecutive surveys") +
#   xlab("Time since site creation (years)") +
#   annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)
# 
# segment_lengths_plot

# pdf(file = file.path("plots", "segment_lengths.pdf"), width = 12, height = 8)
# segment_lengths_plot
# dev.off()







