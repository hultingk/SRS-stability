### looking at dispersal mode info
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg, ape, BiodiversityR) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))
patch_pair_ID <- read.csv(file = file.path("data", "L2_summarized", "patch_pair_ID.csv")) # reading in key to pairs of patches


srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center") %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  ))



# srs_dispersal_prop <- srs_data %>%
#   count(block, patch, patch_type, unique_id, year, time, dispersal_mode) %>%
#   pivot_wider(names_from = dispersal_mode, values_from = n) %>%
#   dplyr::select(!c("NA")) %>%
#   mutate(sp_richness = Animal + Gravity + Wind) %>%
#   mutate(prop_animal = Animal/sp_richness,
#          prop_gravity = Gravity/sp_richness,
#          prop_wind = Wind/sp_richness)
# 
# srs_dispersal_prop <- srs_dispersal_prop %>%
#   pivot_longer(cols = 11:13, names_to = "dispersal_mode", values_to = "proportion") 
# 
# m1 <- glmmTMB(proportion ~ patch_type*time*dispersal_mode + (1|block),
#               data = srs_dispersal_prop,
#               family = "binomial",
#               weights = sp_richness)
# summary(m1)
# plot(simulateResiduals(m1))
# 
# 
# srs_dispersal_prop %>%
#   ggplot(aes(time, proportion, color = patch_type, fill = patch_type)) +
#   facet_wrap(~dispersal_mode) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm") +
#   theme_minimal(base_size = 12) +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") 
# 
# srs_dispersal_prop %>%
#   ggplot(aes(time, sp_richness, color = patch_type, fill = patch_type)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm") +
#   theme_minimal(base_size = 12) +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") 


###### PCOA by dispersal mode #######
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
  select(unique_id, time, year)

# species matrix
animal_sp_info <- animal_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist <- vegdist(animal_sp_info, method = "jaccard")
animal_jaccard_dist_all <- as.matrix(animal_jaccard_dist)
animal_jaccard_dist_all_df <- as.data.frame(animal_jaccard_dist_all)

animal_jaccard_dist_all_df <- cbind(animal_patch_info, animal_jaccard_dist_all_df) # merge with patch info

# pcoa
animal_pcoa <- pcoa(animal_jaccard_dist)

# add axes to patch and time
animal_pcoa_axes <- animal_pcoa$vectors[,c(1,2)]
animal_pcoa_axes <- cbind(animal_patch_info, animal_pcoa_axes)


## gravity dispersed 
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
  select(unique_id, time, year)

# species matrix
gravity_sp_info <- gravity_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
gravity_jaccard_dist <- vegdist(gravity_sp_info, method = "jaccard")
gravity_jaccard_dist_all <- as.matrix(gravity_jaccard_dist)
gravity_jaccard_dist_all_df <- as.data.frame(gravity_jaccard_dist_all)

gravity_jaccard_dist_all_df <- cbind(gravity_patch_info, gravity_jaccard_dist_all_df) # merge with patch info

# pcoa
gravity_pcoa <- pcoa(gravity_jaccard_dist)

# add axes to patch and time
gravity_pcoa_axes <- gravity_pcoa$vectors[,c(1,2)]
gravity_pcoa_axes <- cbind(gravity_patch_info, gravity_pcoa_axes)



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
  select(unique_id, time, year)

# species matrix
wind_sp_info <- wind_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
wind_jaccard_dist <- vegdist(wind_sp_info, method = "jaccard")
wind_jaccard_dist_all <- as.matrix(wind_jaccard_dist)
wind_jaccard_dist_all_df <- as.data.frame(wind_jaccard_dist_all)

wind_jaccard_dist_all_df <- cbind(wind_patch_info, wind_jaccard_dist_all_df) # merge with patch info

# pcoa
wind_pcoa <- pcoa(wind_jaccard_dist)

# add axes to patch and time
wind_pcoa_axes <- wind_pcoa$vectors[,c(1,2)]
wind_pcoa_axes <- cbind(wind_patch_info, wind_pcoa_axes)


#### Convergence/divergence between patch types ####
##### animal dispersed convergence/divergence #####
# use PCoA axes
animal_pcoa_dist_bw_patch <- animal_pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  mutate(block_time = paste(block, time, sep = "-"))

# separating by patch replicate
animal_dist_bw_b <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
  rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
animal_dist_bw_c <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
  rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
animal_dist_bw_d <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
  rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
animal_dist_bw_e <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
  rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)


# joining all together
animal_dist_bw_all <- animal_dist_bw_b %>%
  left_join(animal_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
  left_join(animal_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
  left_join(animal_dist_bw_e, by = c("block_time", "block", "time", "year"))

# calculate distance using pythagorean theorem
animal_dist_bw_all <- animal_dist_bw_all %>%
  mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
  mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
  mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
  mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
  mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
  mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))

# putting together
animal_dist_bw_all <- animal_dist_bw_all %>%
  select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
animal_dist_bw_all <- animal_dist_bw_all %>%
  pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
               names_to = "patch_pair", values_to = "distance")
animal_dist_bw_all$time <- as.numeric(animal_dist_bw_all$time)


# joining to data 
animal_dist_bw_all <- animal_dist_bw_all %>%
  left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
  filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  mutate(dispersal_mode = "Animal")

animal_dist_bw_all %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Distance between patch type communities") #+
 # annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)




##### gravity dispersed convergence/divergence #####
# use PCoA axes
gravity_pcoa_dist_bw_patch <- gravity_pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  mutate(block_time = paste(block, time, sep = "-"))

# separating by patch replicate
gravity_dist_bw_b <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
  rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
gravity_dist_bw_c <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
  rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
gravity_dist_bw_d <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
  rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
gravity_dist_bw_e <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
  rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)


# joining all together
gravity_dist_bw_all <- gravity_dist_bw_b %>%
  left_join(gravity_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
  left_join(gravity_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
  left_join(gravity_dist_bw_e, by = c("block_time", "block", "time", "year"))

# calculate distance using pythagorean theorem
gravity_dist_bw_all <- gravity_dist_bw_all %>%
  mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
  mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
  mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
  mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
  mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
  mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))

# putting together
gravity_dist_bw_all <- gravity_dist_bw_all %>%
  select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
gravity_dist_bw_all <- gravity_dist_bw_all %>%
  pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
               names_to = "patch_pair", values_to = "distance")
gravity_dist_bw_all$time <- as.numeric(gravity_dist_bw_all$time)


# joining to data 
gravity_dist_bw_all <- gravity_dist_bw_all %>%
  left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
  filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  mutate(dispersal_mode = "Gravity")

gravity_dist_bw_all %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Distance between patch type communities") #+
# annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)





##### wind dispersed convergence/divergence #####
# use PCoA axes
wind_pcoa_dist_bw_patch <- wind_pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  mutate(block_time = paste(block, time, sep = "-"))

# separating by patch replicate
wind_dist_bw_b <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
  rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
wind_dist_bw_c <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
  rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
wind_dist_bw_d <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
  rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
wind_dist_bw_e <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
  rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)


# joining all together
wind_dist_bw_all <- wind_dist_bw_b %>%
  left_join(wind_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
  left_join(wind_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
  left_join(wind_dist_bw_e, by = c("block_time", "block", "time", "year"))

# calculate distance using pythagorean theorem
wind_dist_bw_all <- wind_dist_bw_all %>%
  mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
  mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
  mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
  mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
  mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
  mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))

# putting together
wind_dist_bw_all <- wind_dist_bw_all %>%
  select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
wind_dist_bw_all <- wind_dist_bw_all %>%
  pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
               names_to = "patch_pair", values_to = "distance")
wind_dist_bw_all$time <- as.numeric(wind_dist_bw_all$time)


# joining to data 
wind_dist_bw_all <- wind_dist_bw_all %>%
  left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
  filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  mutate(dispersal_mode = "Wind")

wind_dist_bw_all %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Distance between patch type communities") #+
# annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)



##### convergence all dispersal modes together #####
dispersal_mode_convergence <- rbind(
  wind_dist_bw_all, gravity_dist_bw_all, animal_dist_bw_all
)


dispersal_mode_convergence_plot <- dispersal_mode_convergence %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  facet_wrap(~dispersal_mode) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab(expression(atop("Distance between ", paste("patch type communities")))) 
dispersal_mode_convergence_plot

# pdf(file = file.path("plots", "dispersal_convergence.pdf"), width = 14, height = 6)
# dispersal_mode_convergence_plot
# dev.off()


##### CTA segments ####
##### animal CTA segments #####
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



animal_segment_lengths %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  theme_minimal(base_size = 28) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") #+
  #annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)





##### gravity CTA segments #####
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



gravity_segment_lengths %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  theme_minimal(base_size = 28) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") #+
#annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)







##### wind CTA segments #####
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


wind_segment_lengths %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  theme_minimal(base_size = 28) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") #+
#annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)


##### CTA segments all dispersal modes together #####
dispersal_mode_segments <- rbind(
  wind_segment_lengths, gravity_segment_lengths, animal_segment_lengths
)

dispersal_mode_segments_plot <- dispersal_mode_segments %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  theme_minimal(base_size = 28) +
  facet_wrap(~dispersal_mode) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") #+
#annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)
dispersal_mode_segments_plot

# pdf(file = file.path("plots", "dispersal_segments.pdf"), width = 14, height = 6)
# dispersal_mode_segments_plot
# dev.off()



##### Directionality ####
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


animal_direction_all %>%
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
  xlab("Time period") #+
 # annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)






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


gravity_direction_all %>%
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
  xlab("Time period") #+
# annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)





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

wind_direction_all %>%
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
  xlab("Time period") #+
# annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)


###### Directionality all dispersal modes together #####
directionality_all <- rbind(
  wind_direction_all,
  gravity_direction_all,
  animal_direction_all
)



directionality_all %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, directionality, fill = patch_type)) +
  geom_boxplot(aes(fill = patch_type)) +
  theme_minimal(base_size = 28) +
  facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  ylab("Trajectory directionality") +
  xlab("Time period") 


#model
m.direction_dispersal <- glmmTMB(directionality ~ patch_type * dispersal_mode * time + (1|block/patch),
                       data = directionality_all)
summary(m.direction_dispersal)
plot(simulateResiduals(m.direction_dispersal))
performance::r2(m_length_quad)

m.direction_dispersal.posthoc <- emmeans(m.direction_dispersal, ~ patch_type*time*dispersal_mode)
pairs(m.direction_dispersal.posthoc, simple = "time")
# trends are primarily driven by wind dispersed species


# prediction plot
m.direction_dispersal.predict <- ggpredict(m.direction_dispersal, terms=c("time [all]", "patch_type [all]", "dispersal_mode [all]"), back_transform = T)
m.direction_dispersal.predict$dispersal_mode <- m.direction_dispersal.predict$facet

dispersal_direction_plot <- m.direction_dispersal.predict %>%
  ggplot() +
  geom_point(aes(x = x, y = predicted, color = group), size = 10, data = m.direction_dispersal.predict,  position = position_dodge(width = 0.7))+ 
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group),
                data = m.direction_dispersal.predict, width = 0, linewidth = 4.5,  position = position_dodge(width = 0.7)) +
  theme_minimal(base_size = 28) +
  facet_wrap(~dispersal_mode) +
  geom_jitter(aes(x = time, y = directionality, color = patch_type), 
              data = directionality_all, alpha = 0.2, size = 7, 
              position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  #geom_text(data = m1.stat.test, aes(x = ptype, y = height, label = significance), size = 6) +
  labs(title = NULL,
       x = NULL,
       y = "Trajectory directionality") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") #+
  #annotate("text", x = 2.25, y=0.395, label = expression(paste('R'^2*' = 0.431')), size=7)
dispersal_direction_plot

# pdf(file = file.path("plots", "dispersal_direction.pdf"), width = 15, height = 7)
# dispersal_direction_plot
# dev.off()











