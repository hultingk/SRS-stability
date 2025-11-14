########
## SCRIPT NAME: 04_cta_segments.R
## AUTHOR: Katherine Hulting
## PURPOSE: Calculate interannual trajectory distances, repeat within dispersal mode groups
## PRODUCTS: 
#########

### community trajectory analysis - segment lengths
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot, kableExtra) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% 
  filter(transplant != TRUE) %>% # removing experimentally planted species 
  filter(patch_type != "Center") # removing center patch from analysis



########################
#### ALL SPECIES ####
########################
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
  dplyr::select(!time)

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys <- patch_info %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths <- cbind(segment_lengths, time_surveys)
segment_lengths$dispersal_mode <- "All Species"
segment_lengths$s.time <- as.numeric(scale(segment_lengths$time)) # scaling time


###### MODELS ######
# segment lengths model 
# linear
m_length <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                    data = segment_lengths)
# quadratic
m_length_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                    data = segment_lengths)
# null
m_length_null <- glmmTMB(distance ~ 1 + (1|block/patch), # null model
                         data = segment_lengths)

# AIC comparison
a <- list(m_length, m_length_quad, m_length_null)
aictab(a) # quadratic much better fit

# model fit
summary(m_length_quad)
plot(simulateResiduals(m_length_quad))
check_model(m_length_quad)
performance::r2(m_length_quad)

# posthoc
m_length_posthoc <- emmeans(m_length_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0, 1.8)))
m_length_pairs <- pairs(m_length_posthoc, simple = "patch_type")
m_length_pairs

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





posthoc_table <- as.data.frame(m_length_pairs)
posthoc_table %>%
  dplyr::select(-s.time) %>%
  kbl(digits = 3, caption = "Post-hoc Pairwise Comparisons (emmeans)") %>%
  kable_classic(full_width = FALSE) %>%
  kable_styling(html_font = "Times New Roman")




# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

m_length_predict <- ggpredict(m_length_quad, terms = c("s.time [all]", "patch_type"))
m_length_predict$group <- factor(m_length_predict$group, levels = c("Connected", "Rectangular", "Winged"))
m_length_predict <- m_length_predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

segment_lengths_plot <- m_length_predict %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 5.5, alpha = 0.19, data = segment_lengths) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme_minimal(base_size = 24) +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") 
  #annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)
segment_lengths_plot

# pdf(file = file.path("plots", "segment_lengths.pdf"), width = 12, height = 8)
# segment_lengths_plot
# dev.off()



#####################
#### ANIMAL ####
#####################

##### animal CTA segments #####
## animal dispersed 
animal_data <- srs_data %>%
  filter(dispersal_mode == "Animal") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
animal_data$time <- as.factor(animal_data$time)
animal_data$unique_id <- as.factor(animal_data$unique_id)
animal_data$year <- as.factor(animal_data$year)

# patch data
animal_patch_info <- animal_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_sp_info <- animal_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist <- vegdist(animal_sp_info, method = "jaccard")

animal_patch_info$time <- as.numeric(animal_patch_info$time)
# defining trajectories
animal_srs_trajectory <- defineTrajectories(animal_jaccard_dist, sites = animal_patch_info$unique_id, surveys = animal_patch_info$time)

# segment lengths of trajectories between consectutive years
animal_segment_lengths <- trajectoryLengths(animal_srs_trajectory)
animal_segment_lengths <- animal_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
animal_time_surveys <- animal_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  #filter(!(block == "08" & time == 25)) %>%
  # filter(!(block == "10" & time == 25)) %>%
  filter(!(block == "52" & time == 18)) %>%
  #filter(!(block == "54N" & time == 18)) %>%
  filter(!(block == "57" & time == 18)) %>%
  filter(!(block == "75E" & time == 18)) %>%
  filter(!(block == "75W" & time == 18)) %>%
  # filter(!(unique_id == "53N-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "53S-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "54N-E-rectangle" & time == 17)) %>%
  # filter(!(unique_id == "54S-C-rectangle" & time == 25)) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))

# joining with segment lengths
animal_segment_lengths <- cbind(animal_segment_lengths, animal_time_surveys)
animal_segment_lengths$dispersal_mode <- "Animal"
animal_segment_lengths$s.time <- as.numeric(scale(animal_segment_lengths$time)) # scaling time


###### MODELS ######
# segment lengths model 
# quadratic
m_length_animal_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                         data = animal_segment_lengths)

# model fit
summary(m_length_animal_quad)
plot(simulateResiduals(m_length_animal_quad))
check_model(m_length_animal_quad)
performance::r2(m_length_animal_quad)

# posthoc
m_length_animal_posthoc <- emmeans(m_length_animal_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0, 1.8)))
m_length_animal_pairs <- pairs(m_length_animal_posthoc, simple = "patch_type")










#####################
#### GRAVITY ####
#####################

##### gravity CTA segments #####
gravity_data <- srs_data %>%
  filter(dispersal_mode == "Gravity") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
gravity_data$time <- as.factor(gravity_data$time)
gravity_data$unique_id <- as.factor(gravity_data$unique_id)
gravity_data$year <- as.factor(gravity_data$year)

# patch data
gravity_patch_info <- gravity_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_sp_info <- gravity_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))
# Jaccard distance matrix
gravity_jaccard_dist <- vegdist(gravity_sp_info, method = "jaccard")

gravity_patch_info$time <- as.numeric(gravity_patch_info$time)
# defining trajectories
gravity_srs_trajectory <- defineTrajectories(gravity_jaccard_dist, sites = gravity_patch_info$unique_id, surveys = gravity_patch_info$time)

# segment lengths of trajectories between consectutive years
gravity_segment_lengths <- trajectoryLengths(gravity_srs_trajectory)
gravity_segment_lengths <- gravity_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
gravity_time_surveys <- gravity_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  #filter(!(block == "08" & time == 25)) %>%
  # filter(!(block == "10" & time == 25)) %>%
  filter(!(block == "52" & time == 18)) %>%
  #filter(!(block == "54N" & time == 18)) %>%
  filter(!(block == "57" & time == 18)) %>%
  filter(!(block == "75E" & time == 18)) %>%
  filter(!(block == "75W" & time == 18)) %>%
  # filter(!(unique_id == "53N-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "53S-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "54N-E-rectangle" & time == 17)) %>%
  # filter(!(unique_id == "54S-C-rectangle" & time == 25)) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))


# joining with segment lengths
gravity_segment_lengths <- cbind(gravity_segment_lengths, gravity_time_surveys)
gravity_segment_lengths$dispersal_mode <- "Gravity"
gravity_segment_lengths$s.time <- as.numeric(scale(gravity_segment_lengths$time)) # scaling time


###### MODELS ######
# segment lengths model 
# quadratic
m_length_gravity_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                                data = gravity_segment_lengths)

# model fit
summary(m_length_gravity_quad)
plot(simulateResiduals(m_length_gravity_quad))
check_model(m_length_gravity_quad)
performance::r2(m_length_gravity_quad)

# posthoc
m_length_gravity_posthoc <- emmeans(m_length_gravity_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0, 1.8)))
m_length_gravity_pairs <- pairs(m_length_gravity_posthoc, simple = "patch_type")









#####################
#### WIND ####
#####################

##### wind CTA segments #####
## wind dispersed 
wind_data <- srs_data %>%
  filter(dispersal_mode == "Wind") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
wind_data$time <- as.factor(wind_data$time)
wind_data$unique_id <- as.factor(wind_data$unique_id)
wind_data$year <- as.factor(wind_data$year)

# patch data
wind_patch_info <- wind_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_sp_info <- wind_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))
# Jaccard distance matrix
wind_jaccard_dist <- vegdist(wind_sp_info, method = "jaccard")

wind_patch_info$time <- as.numeric(wind_patch_info$time)
# defining trajectories
wind_srs_trajectory <- defineTrajectories(wind_jaccard_dist, sites = wind_patch_info$unique_id, surveys = wind_patch_info$time)

# segment lengths of trajectories between consectutive years
wind_segment_lengths <- trajectoryLengths(wind_srs_trajectory)
wind_segment_lengths <- wind_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
wind_time_surveys <- wind_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  #filter(!(block == "08" & time == 25)) %>%
  # filter(!(block == "10" & time == 25)) %>%
  filter(!(block == "52" & time == 18)) %>%
  #filter(!(block == "54N" & time == 18)) %>%
  filter(!(block == "57" & time == 18)) %>%
  filter(!(block == "75E" & time == 18)) %>%
  filter(!(block == "75W" & time == 18)) %>%
  # filter(!(unique_id == "53N-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "53S-E-wing" & time == 25)) %>%
  #filter(!(unique_id == "54N-E-rectangle" & time == 17)) %>%
  # filter(!(unique_id == "54S-C-rectangle" & time == 25)) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))


# joining with segment lengths
wind_segment_lengths <- cbind(wind_segment_lengths, wind_time_surveys)
wind_segment_lengths$dispersal_mode <- "Wind"
wind_segment_lengths$s.time <- as.numeric(scale(wind_segment_lengths$time)) # scaling time


###### MODELS ######
# segment lengths model 
# quadratic
m_length_wind_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                                 data = wind_segment_lengths)

# model fit
summary(m_length_wind_quad)
plot(simulateResiduals(m_length_wind_quad))
check_model(m_length_wind_quad)
performance::r2(m_length_wind_quad)

# posthoc
m_length_wind_posthoc <- emmeans(m_length_wind_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0, 1.8)))
m_length_wind_pairs <- pairs(m_length_wind_posthoc, simple = "patch_type")





########################
#### TABLES ####
########################
# emmeans posthoc tables

# creating dataframes of results
# all species
m_length_pairs_df <- as.data.frame(m_length_pairs)
m_length_pairs_df <- m_length_pairs_df %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(Time = if_else(s.time == 0, "Midpoint of time series (12 years)", "End of time series (24 years)")) 
# animal dispersed
m_length_animal_pairs_df <- as.data.frame(m_length_animal_pairs)
m_length_animal_pairs_df <- m_length_animal_pairs_df %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(Time = if_else(s.time == 0, "Midpoint of time series (12 years)", "End of time series (24 years)")) 
# gravity dispersed
m_length_gravity_pairs_df <- as.data.frame(m_length_gravity_pairs)
m_length_gravity_pairs_df <- m_length_gravity_pairs_df %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") %>%
  mutate(Time = if_else(s.time == 0, "Midpoint of time series (12 years)", "End of time series (24 years)")) 
# wind dispersed
m_length_wind_pairs_df <- as.data.frame(m_length_wind_pairs)
m_length_wind_pairs_df <- m_length_wind_pairs_df %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(Time = if_else(s.time == 0, "Midpoint of time series (12 years)", "End of time series (24 years)")) 

m_length_table_all <- rbind(
  m_length_pairs_df, m_length_animal_pairs_df, m_length_gravity_pairs_df, m_length_wind_pairs_df
)

tableS2 <- m_length_table_all %>% 
  dplyr::select(`Dispersal mode`, Time, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  row_spec(seq(3, nrow(m_length_table_all), 3), extra_css = "border-bottom: 2px solid;") %>%
  row_spec(seq(6, nrow(m_length_table_all), 6), extra_css = "border-bottom: 5px double;") %>%
  row_spec(1, extra_css = "border-top: 5px double;") %>%
  row_spec(0:nrow(m_length_table_all), extra_css = "padding-bottom: 10px;") 
tableS2

# exporting
# save_kable(tableS2, file = file.path("plots", "tableS2.html"))




####################
#### PLOTS ####
###################

# model predictions
# All species predictions
m_length_predict <- ggpredict(m_length_quad, terms = c("s.time [all]", "patch_type"))
m_length_predict <- as.data.frame(m_length_predict)
m_length_predict$dispersal_mode <- "All Species"
scaled_time_key <- segment_lengths %>% # creating key of scaled times to join to predictions for easy visualization
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m_length_predict <- m_length_predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

# Animal dispersed predictions
m.animal_segments.predict <- ggpredict(m_length_animal_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.animal_segments.predict <- as.data.frame(m.animal_segments.predict)
m.animal_segments.predict$dispersal_mode <- "Animal"
animal_time_key <- animal_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.animal_segments.predict <- m.animal_segments.predict %>%
  left_join(animal_time_key, by = c("x" = "s.time"))

m.gravity_segments.predict <- ggpredict(m_length_gravity_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.gravity_segments.predict <- as.data.frame(m.gravity_segments.predict)
m.gravity_segments.predict$dispersal_mode <- "Gravity"
gravity_time_key <- gravity_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.gravity_segments.predict <- m.gravity_segments.predict %>%
  left_join(gravity_time_key, by = c("x" = "s.time"))

m.wind_segments.predict <- ggpredict(m_length_wind_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.wind_segments.predict <- as.data.frame(m.wind_segments.predict)
m.wind_segments.predict$dispersal_mode <- "Wind"
wind_time_key <- wind_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.wind_segments.predict <- m.wind_segments.predict %>%
  left_join(wind_time_key, by = c("x" = "s.time"))


# FACET BY ROWS - total and animal together and gravity and wind together
# joining together predictions
predict_segments_1 <- rbind(
  m_length_predict, m.animal_segments.predict
)
predict_segments_2 <- rbind(
  m.gravity_segments.predict, m.wind_segments.predict
)

predict_segments_1$dispersal_mode <- factor(predict_segments_1$dispersal_mode, levels = c("All Species", "Animal"))
predict_segments_2$dispersal_mode <- factor(predict_segments_2$dispersal_mode, levels = c("Gravity", "Wind"))

# joining together data points
segment_lengths <- segment_lengths %>%
  dplyr::select(-soil_moisture, -year_since_fire)

dispersal_mode_segments_1 <- rbind(
  segment_lengths, animal_segment_lengths
)
dispersal_mode_segments_1$dispersal_mode <- factor(dispersal_mode_segments_1$dispersal_mode, levels = c("All Species", "Animal"))

dispersal_mode_segments_2 <- rbind(
  gravity_segment_lengths, wind_segment_lengths
)
dispersal_mode_segments_2$dispersal_mode <- factor(dispersal_mode_segments_2$dispersal_mode, levels = c("Gravity", "Wind"))

# two faceted plots
# first set of plots
segments_plot_1 <- predict_segments_1 %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 3, alpha = 0.15, data = dispersal_mode_segments_1) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.4) +
  facet_wrap(~dispersal_mode, scales = "free") +
  ylim(c(0, 0.5)) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab(NULL) +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "none") 
segments_plot_1

# second set of plots
segments_plot_2 <- predict_segments_2 %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 3, alpha = 0.15, data = dispersal_mode_segments_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.4) +
  facet_wrap(~dispersal_mode, scales = "free") +
  ylim(c(0, 0.6)) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Years since site creation") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "none") 
segments_plot_2

# get legend
pL <- predict_segments_2 %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 4, alpha = 0.2, data = dispersal_mode_segments_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
l <- get_legend(pL)

# put together
figure3 <- cowplot::plot_grid(segments_plot_1, l, segments_plot_2, 
                                      ncol = 2, nrow = 2, rel_widths = c(1, 0.3), rel_heights = c(1, 1.1),
                                      label_size = 20, label_x = 0.2, label_y = 0.95)
figure3
# exporting
pdf(file = file.path("plots", "figure3.pdf"), width = 10.5, height = 9)
figure3
dev.off()

# 
# ### plotting model predictions
# # creating key of scaled times to join to predictions for easy visualization
# scaled_time_key <- dispersal_mode_segments %>%
#   count(time, s.time) %>%
#   dplyr::select(-n) %>%
#   mutate(s.time = round(s.time, 2))
# 
# # model predictions
# m.dispersal_segments.predict <- ggpredict(m.dispersal_segments_quad, terms=c("s.time [all]", "patch_type [all]", "dispersal_mode [all]"), back_transform = T)
# m.dispersal_segments.predict <- as.data.frame(m.dispersal_segments.predict)
# 
# # plotting
# dispersal_segments_plot <- m.dispersal_segments.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   rename(dispersal_mode = facet) %>%
#   ggplot() +
#   geom_point(aes(time, distance, color = patch_type), size = 3.5, alpha = 0.15, data = dispersal_mode_segments) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 3) +
#   theme_minimal(base_size = 24) +
#   facet_wrap(~dispersal_mode) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) 
# dispersal_segments_plot
# 
# 
# pdf(file = file.path("plots", "dispersal_segments.pdf"), width = 14, height = 6)
# dispersal_segments_plot
# dev.off()
# 
# 
# 
# 
# 
# 
