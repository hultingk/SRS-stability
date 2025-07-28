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
m.direction <- glmmTMB(directionality ~ patch_type * time + (1|block/patch),
                       data = segment_direction_all)
summary(m.direction)
plot(simulateResiduals(m.direction))
performance::r2(m_length_quad)

m.direction.posthoc <- emmeans(m.direction, ~ patch_type*time)
pairs(m.direction.posthoc)
m.direction_cld <- cld(object = m.direction.posthoc,
                       adjust = "Tukey",
                       Letters = letters,
                       alpha = 0.05)
m.direction_cld

### plotting 
segment_direction_plot <- segment_direction_all %>%
  #filter(!block %in% c("75E", "75W")) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, directionality, fill = patch_type)) +
  #geom_point(aes(color = patch_type), size=2, alpha=0.9, position = position_jitterdodge(dodge.width = 1)) +
  geom_boxplot(aes(fill = patch_type)) +
  theme_minimal(base_size = 28) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  #scale_fill_brewer(palette = "Set2", name = "Patch Type") +
  ylab("Trajectory directionality") +
  xlab("Time period") +
  annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)
segment_direction_plot

pdf(file = file.path("plots", "segment_direction.pdf"), width = 12, height = 8)
segment_direction_plot
dev.off()









###### trajectory angles instead? ######
srs_angles <- trajectoryAngles(srs_trajectory)
srs_angles <- srs_angles %>%
  rownames_to_column("uniqueID")
srs_angles <- srs_angles %>%
  pivot_longer(cols = 2:25, names_to = "time_period", values_to = "angle")

# assinging last year of sequence as time
srs_angles <- srs_angles %>%
  filter(!time_period %in% c("mean", "sd", "rho")) %>%
  mutate(time_period = dplyr::case_when(
    time_period %in% c("S1-S2") ~ 2,
    time_period %in% c("S2-S3") ~ 3,
    time_period %in% c("S3-S4") ~ 4,
    time_period %in% c("S4-S5") ~ 5,
    time_period %in% c("S5-S6") ~ 6,
    time_period %in% c("S6-S7") ~ 7,
    time_period %in% c("S7-S8") ~ 8,
    time_period %in% c("S8-S9") ~ 9,
    time_period %in% c("S9-S10") ~ 10,
    time_period %in% c("S10-S11") ~ 11,
    time_period %in% c("S11-S12") ~ 12,
    time_period %in% c("S12-S13") ~ 13,
    time_period %in% c("S13-S14") ~ 14,
    time_period %in% c("S14-S15") ~ 15,
    time_period %in% c("S15-S16") ~ 16,
    time_period %in% c("S16-S17") ~ 17,
    time_period %in% c("S17-S18") ~ 18,
    time_period %in% c("S18-S19") ~ 19,
    time_period %in% c("S19-S20") ~ 20,
    time_period %in% c("S20-S21") ~ 21,
    time_period %in% c("S21-S22") ~ 22,
    .default = 0
  ))
  
srs_angles <- srs_angles %>%
  separate(uniqueID, into = c("block", "patch_rep", "patch_type"))

# plotting
srs_angle_plot <- srs_angles %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time_period, angle, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 28) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Angle")
srs_angle_plot

pdf(file = file.path("plots", "srs_angle_plot.pdf"), width = 12, height = 8)
srs_angle_plot
dev.off()


# model comparison
m_angle_linear <- glmmTMB(angle ~ patch_type + time_period + (1|block),
                        data = srs_angles)
summary(m_angle_linear)

m_angle_quad <- glmmTMB(angle ~ patch_type * time_period + patch_type * I(time_period^2) + (1|block),
              data = srs_angles)

m_angle_asym <- nls(angle ~ SSasymp(time_period, Asym, R0, lrc),
              data = srs_angles)
summary(m_angle_asym)


# AIC comparison
a <- AIC(m_angle_quad, m_angle_asym, m_angle_linear)




