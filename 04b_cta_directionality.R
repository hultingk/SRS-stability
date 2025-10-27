### community trajectory analysis - directionality
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, 
                 ggeffects, AICcmodavg, multcomp, multcompView) # Install missing packages and load needed libraries

source(here::here("04a_cta_segments.R"))


##### directionality broken into time periods ####
###### directionality in first 12 years #####
sp_info_1_12 <- srs_data_wider %>%
  filter(time <= 12)

patch_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_1_12 <- vegdist(sp_info_1_12, method = "jaccard")

# defining trajectories
srs_trajectory_1_12 <- defineTrajectories(jaccard_dist_1_12, sites = patch_info_1_12$unique_id, surveys = patch_info_1_12$time)

# directionality
segment_direction_1_12 <- trajectoryDirectionality(srs_trajectory_1_12)
segment_direction_1_12 <- data.frame(segment_direction_1_12)
segment_direction_1_12 <- segment_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = segment_direction_1_12)


###### directionality in second 12 years #####
sp_info_13_24 <- srs_data_wider %>%
  filter(time >= 13)

patch_info_13_24 <- sp_info_13_24 %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_13_24 <- sp_info_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_13_24 <- vegdist(sp_info_13_24, method = "jaccard")

# defining trajectories
srs_trajectory_13_24 <- defineTrajectories(jaccard_dist_13_24, sites = patch_info_13_24$unique_id, surveys = patch_info_13_24$time)

# directionality
segment_direction_13_24 <- trajectoryDirectionality(srs_trajectory_13_24)
segment_direction_13_24 <- data.frame(segment_direction_13_24)
segment_direction_13_24 <- segment_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = segment_direction_13_24)


#### putting all together
segment_direction_all <- rbind(
  segment_direction_1_12,
  segment_direction_13_24
)


# model
m.direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = segment_direction_all)
summary(m.direction)
plot(simulateResiduals(m.direction))
performance::r2(m.direction)

# percent change from decade 1 to decade 2
(-0.028890)/0.369338 * 100 # -7.822103% decrease in directionality
confint(m.direction)
(-0.042316593)/(0.359607329) *100 # -11.76744%
(-0.0154639504)/(0.3790690622) *100 # -4.079455%


m.direction.posthoc <- emmeans(m.direction, ~ patch_type*time)
pairs(m.direction.posthoc, simple = "patch_type")
m.direction_cld <- cld(object = m.direction.posthoc,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
m.direction_cld

### plotting 
# segment_direction_plot <- segment_direction_all %>%
#   #filter(!block %in% c("75E", "75W")) %>%
#   mutate(patch_type = dplyr::case_when(
#     patch_type %in% c("connected") ~ "Connected",
#     patch_type %in% c("rectangle") ~ "Rectangular",
#     patch_type %in% c("wing") ~ "Winged"
#   )) %>%
#   ggplot(aes(time, directionality, fill = patch_type)) +
#   #geom_point(aes(color = patch_type), size=2, alpha=0.9, position = position_jitterdodge(dodge.width = 1)) +
#   geom_boxplot(aes(fill = patch_type)) +
#   theme_minimal(base_size = 28) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   #scale_fill_brewer(palette = "Set2", name = "Patch Type") +
#   ylab("Trajectory directionality") +
#   xlab("Time period") +
#   annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)
# segment_direction_plot

# prediction plot
m.direction.predict <- ggpredict(m.direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
segment_direction_plot <- m.direction.predict %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type), 
              data = segment_direction_all, alpha = 0.2, size = 7, 
              position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = m.direction.predict, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
  theme_minimal(base_size = 28) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 10, 
             data = m.direction.predict,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  labs(title = NULL,
       x = NULL,
       y = "Trajectory directionality") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") #+
 # annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.455')), size=7)
segment_direction_plot



# pdf(file = file.path("plots", "segment_direction.pdf"), width = 12, height = 8)
# segment_direction_plot
# dev.off()

#### dispersal mode ####
###### animal directionality #####
animal_data$time <- as.numeric(animal_data$time)
# first 12 years
animal_1_12 <- animal_data %>%
  filter(time <= 12)

# patch info
animal_patch_info_1_12 <- animal_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_1_12 <- animal_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist_1_12 <- vegdist(animal_1_12, method = "jaccard")
# defining trajectories
animal_trajectory_1_12 <- defineTrajectories(animal_jaccard_dist_1_12, sites = animal_patch_info_1_12$unique_id, surveys = animal_patch_info_1_12$time)

# directionality
animal_direction_1_12 <- trajectoryDirectionality(animal_trajectory_1_12)
animal_direction_1_12 <- data.frame(animal_direction_1_12)
animal_direction_1_12 <- animal_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = animal_direction_1_12)


# second half of succession
animal_13_24 <- animal_data %>%
  filter(time >= 13)

# patch info
animal_patch_info_13_24 <- animal_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_13_24 <- animal_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist_13_24 <- vegdist(animal_13_24, method = "jaccard")
# defining trajectories
animal_trajectory_13_24 <- defineTrajectories(animal_jaccard_dist_13_24, sites = animal_patch_info_13_24$unique_id, surveys = animal_patch_info_13_24$time)

# directionality
animal_direction_13_24 <- trajectoryDirectionality(animal_trajectory_13_24)
animal_direction_13_24 <- data.frame(animal_direction_13_24)
animal_direction_13_24 <- animal_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = animal_direction_13_24)

#### putting all together
animal_direction_all <- rbind(
  animal_direction_1_12,
  animal_direction_13_24
)
animal_direction_all$dispersal_mode <- "Animal"


###### gravity directionality #####
gravity_data$time <- as.numeric(gravity_data$time)
# first 12 years
gravity_1_12 <- gravity_data %>%
  filter(time <= 12)

# patch info
gravity_patch_info_1_12 <- gravity_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_1_12 <- gravity_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
gravity_jaccard_dist_1_12 <- vegdist(gravity_1_12, method = "jaccard")
# defining trajectories
gravity_trajectory_1_12 <- defineTrajectories(gravity_jaccard_dist_1_12, sites = gravity_patch_info_1_12$unique_id, surveys = gravity_patch_info_1_12$time)

# directionality
gravity_direction_1_12 <- trajectoryDirectionality(gravity_trajectory_1_12)
gravity_direction_1_12 <- data.frame(gravity_direction_1_12)
gravity_direction_1_12 <- gravity_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = gravity_direction_1_12)


# second half of succession
gravity_13_24 <- gravity_data %>%
  filter(time >= 13)

# patch info
gravity_patch_info_13_24 <- gravity_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_13_24 <- gravity_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
gravity_jaccard_dist_13_24 <- vegdist(gravity_13_24, method = "jaccard")
# defining trajectories
gravity_trajectory_13_24 <- defineTrajectories(gravity_jaccard_dist_13_24, sites = gravity_patch_info_13_24$unique_id, surveys = gravity_patch_info_13_24$time)

# directionality
gravity_direction_13_24 <- trajectoryDirectionality(gravity_trajectory_13_24)
gravity_direction_13_24 <- data.frame(gravity_direction_13_24)
gravity_direction_13_24 <- gravity_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = gravity_direction_13_24)

#### putting all together
gravity_direction_all <- rbind(
  gravity_direction_1_12,
  gravity_direction_13_24
)
gravity_direction_all$dispersal_mode <- "Gravity"






###### wind directionality #####
wind_data$time <- as.numeric(wind_data$time)
# first 12 years
wind_1_12 <- wind_data %>%
  filter(time <= 12)

# patch info
wind_patch_info_1_12 <- wind_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_1_12 <- wind_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
wind_jaccard_dist_1_12 <- vegdist(wind_1_12, method = "jaccard")
# defining trajectories
wind_trajectory_1_12 <- defineTrajectories(wind_jaccard_dist_1_12, sites = wind_patch_info_1_12$unique_id, surveys = wind_patch_info_1_12$time)

# directionality
wind_direction_1_12 <- trajectoryDirectionality(wind_trajectory_1_12)
wind_direction_1_12 <- data.frame(wind_direction_1_12)
wind_direction_1_12 <- wind_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = wind_direction_1_12)


# second half of succession
wind_13_24 <- wind_data %>%
  filter(time >= 13)

# patch info
wind_patch_info_13_24 <- wind_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_13_24 <- wind_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
wind_jaccard_dist_13_24 <- vegdist(wind_13_24, method = "jaccard")
# defining trajectories
wind_trajectory_13_24 <- defineTrajectories(wind_jaccard_dist_13_24, sites = wind_patch_info_13_24$unique_id, surveys = wind_patch_info_13_24$time)

# directionality
wind_direction_13_24 <- trajectoryDirectionality(wind_trajectory_13_24)
wind_direction_13_24 <- data.frame(wind_direction_13_24)
wind_direction_13_24 <- wind_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = wind_direction_13_24)

#### putting all together
wind_direction_all <- rbind(
  wind_direction_1_12,
  wind_direction_13_24
)
wind_direction_all$dispersal_mode <- "Wind"


###### Directionality models #####
segment_direction_all$dispersal_mode <- "Total"
directionality_all <- rbind(
  segment_direction_all,
  animal_direction_all,
  gravity_direction_all,
  wind_direction_all
)


# all species
m.direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = segment_direction_all)
summary(m.direction)
# animal
m.direction.animal <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = animal_direction_all)
summary(m.direction.animal)
# gravity
m.direction.gravity <- glmmTMB(directionality ~ patch_type * time + (1|block),
                              data = gravity_direction_all)
summary(m.direction.gravity)
# wind
m.direction.wind <- glmmTMB(directionality ~ patch_type * time + (1|block),
                               data = wind_direction_all)
summary(m.direction.wind)



# predictions 
m.direction.predict <- ggpredict(m.direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.predict$dispersal_mode <- "Total"
m.direction.animal.predict <- ggpredict(m.direction.animal, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.animal.predict$dispersal_mode <- "Animal"
m.direction.gravity.predict <- ggpredict(m.direction.gravity, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.gravity.predict$dispersal_mode <- "Gravity"
m.direction.wind.predict <- ggpredict(m.direction.wind, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.wind.predict$dispersal_mode <- "Wind"

direction_predict_all <- rbind(
  m.direction.predict, 
  m.direction.animal.predict,
  m.direction.gravity.predict,
  m.direction.wind.predict
)

# reordering factors
direction_predict_all$dispersal_mode <- fct_relevel(direction_predict_all$dispersal_mode, "Total", "Animal", "Gravity", "Wind")
directionality_all$dispersal_mode <- fct_relevel(directionality_all$dispersal_mode, "Total", "Animal", "Gravity", "Wind")

# plotting
direction_predict_plot <- direction_predict_all %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type), 
              data = directionality_all, alpha = 0.2, size = 7, 
              position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = direction_predict_all, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_mode, scales = "free") +
  theme_minimal(base_size = 26) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 10, 
             data = direction_predict_all,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  labs(title = NULL,
       x = NULL,
       y = "Trajectory directionality") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.justification = "top")

pdf(file = file.path("plots", "direction_predict_plot.pdf"), width = 12.5, height = 9)
direction_predict_plot
dev.off()





















# 
# 
# ###### trajectory angles instead? ######
# srs_angles <- trajectoryAngles(srs_trajectory)
# srs_angles <- srs_angles %>%
#   rownames_to_column("uniqueID")
# srs_angles <- srs_angles %>%
#   pivot_longer(cols = 2:25, names_to = "time_period", values_to = "angle")
# 
# # assinging last year of sequence as time
# srs_angles <- srs_angles %>%
#   filter(!time_period %in% c("mean", "sd", "rho")) %>%
#   mutate(time_period = dplyr::case_when(
#     time_period %in% c("S1-S2") ~ 2,
#     time_period %in% c("S2-S3") ~ 3,
#     time_period %in% c("S3-S4") ~ 4,
#     time_period %in% c("S4-S5") ~ 5,
#     time_period %in% c("S5-S6") ~ 6,
#     time_period %in% c("S6-S7") ~ 7,
#     time_period %in% c("S7-S8") ~ 8,
#     time_period %in% c("S8-S9") ~ 9,
#     time_period %in% c("S9-S10") ~ 10,
#     time_period %in% c("S10-S11") ~ 11,
#     time_period %in% c("S11-S12") ~ 12,
#     time_period %in% c("S12-S13") ~ 13,
#     time_period %in% c("S13-S14") ~ 14,
#     time_period %in% c("S14-S15") ~ 15,
#     time_period %in% c("S15-S16") ~ 16,
#     time_period %in% c("S16-S17") ~ 17,
#     time_period %in% c("S17-S18") ~ 18,
#     time_period %in% c("S18-S19") ~ 19,
#     time_period %in% c("S19-S20") ~ 20,
#     time_period %in% c("S20-S21") ~ 21,
#     time_period %in% c("S21-S22") ~ 22,
#     .default = 0
#   ))
#   
# srs_angles <- srs_angles %>%
#   separate(uniqueID, into = c("block", "patch_rep", "patch_type"))
# 
# # plotting
# srs_angle_plot <- srs_angles %>%
#   mutate(patch_type = dplyr::case_when(
#     patch_type %in% c("connected") ~ "Connected",
#     patch_type %in% c("rectangle") ~ "Rectangular",
#     patch_type %in% c("wing") ~ "Winged"
#   )) %>%
#   ggplot(aes(time_period, angle, color = patch_type, fill = patch_type)) +
#   geom_point(size = 4.5, alpha = 0.3) +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 28) +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   theme(plot.margin = margin(1, 1, 1, 2, "cm"))+
#   xlab("Time since site creation (years)") +
#   #ylab(expression(paste("Trajectory directionality \nbetween consecutive years", " (angle ", theta, ")")))
#   ylab("Trajectory directionality \nbetween consecutive years")
# srs_angle_plot
# 
# pdf(file = file.path("plots", "srs_angle_plot.pdf"), width = 12, height = 8)
# srs_angle_plot
# dev.off()
# 
# 
# # model comparison
# m_angle_linear <- glmmTMB(angle ~ patch_type + time_period + (1|block),
#                         data = srs_angles)
# summary(m_angle_linear)
# 
# m_angle_quad <- glmmTMB(angle ~ patch_type * time_period + patch_type * I(time_period^2) + (1|block),
#               data = srs_angles)
# 
# m_angle_asym <- nls(angle ~ SSasymp(time_period, Asym, R0, lrc),
#               data = srs_angles)
# summary(m_angle_asym)
# 
# 
# # AIC comparison
# a <- AIC(m_angle_quad, m_angle_asym, m_angle_linear)
# 
# 
# 
# 
