


### looking at dispersal mode info
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, ape, BiodiversityR, performance) # Install missing packages and load needed libraries

source(here::here("00_functions.R"))
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


srs_data_count <- srs_data %>%
  count(block, patch, patch_type, unique_id, time, dispersal_mode)

srs_data_count %>%
  ggplot() +
  geom_point(aes(time, n, color = dispersal_mode), alpha = 0.3, size = 2) +
  facet_wrap(~patch_type) +
  geom_smooth(aes(time, n, color = dispersal_mode, fill = dispersal_mode), size = 2) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")
  #scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  #scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison")


srs_dispersal_prop <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, dispersal_mode) %>%
  pivot_wider(names_from = dispersal_mode, values_from = n) %>%
  dplyr::select(!c("NA")) %>%
  mutate(sp_richness = Animal + Gravity + Wind) %>%
  mutate(prop_animal = Animal/sp_richness,
         prop_gravity = Gravity/sp_richness,
         prop_wind = Wind/sp_richness)

srs_dispersal_prop <- srs_dispersal_prop %>%
  pivot_longer(cols = 11:13, names_to = "dispersal_mode", values_to = "proportion")


srs_dispersal_prop %>%
  ggplot(aes(time, proportion, color = patch_type, fill = patch_type)) +
  facet_wrap(~dispersal_mode) +
  geom_point(alpha = 0.2, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")



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
  dplyr::select(unique_id, time, year)

# species matrix
animal_sp_info <- animal_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

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
  dplyr::select(unique_id, time, year)

# species matrix
gravity_sp_info <- gravity_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

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
  dplyr::select(unique_id, time, year)

# species matrix
wind_sp_info <- wind_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

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
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
animal_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Animal") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
animal_convergence_jaccard <- animal_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Animal")



# use PCoA axes
# animal_pcoa_dist_bw_patch <- animal_pcoa_axes %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
#   #filter(!block %in% c("75W", "75E")) %>%
#   mutate(block_time = paste(block, time, sep = "-"))
# 
# # separating by patch replicate
# animal_dist_bw_b <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
#   rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
# animal_dist_bw_c <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
#   rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
# animal_dist_bw_d <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
#   rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
# animal_dist_bw_e <- animal_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
#   rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)
# 
# 
# # joining all together
# animal_dist_bw_all <- animal_dist_bw_b %>%
#   left_join(animal_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
#   left_join(animal_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
#   left_join(animal_dist_bw_e, by = c("block_time", "block", "time", "year"))
# 
# # calculate distance using pythagorean theorem
# animal_dist_bw_all <- animal_dist_bw_all %>%
#   mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
#   mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
#   mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
#   mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
#   mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
#   mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))
# 
# # putting together
# animal_dist_bw_all <- animal_dist_bw_all %>%
#   select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
# animal_dist_bw_all <- animal_dist_bw_all %>%
#   pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
#                names_to = "patch_pair", values_to = "distance")
# animal_dist_bw_all$time <- as.numeric(animal_dist_bw_all$time)
# 
# 
# # joining to data 
# animal_dist_bw_all <- animal_dist_bw_all %>%
#   left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
#   filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
#   mutate(dispersal_mode = "Animal")
# 
# animal_dist_bw_all %>%
#   ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 24) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab("Distance between patch type communities") #+
#  # annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)
# 
# 
# 
# 
# ##### gravity dispersed convergence/divergence #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
gravity_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Gravity") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
gravity_convergence_jaccard <- gravity_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Gravity")

# # use PCoA axes
# gravity_pcoa_dist_bw_patch <- gravity_pcoa_axes %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
#   #filter(!block %in% c("75W", "75E")) %>%
#   mutate(block_time = paste(block, time, sep = "-"))
# 
# # separating by patch replicate
# gravity_dist_bw_b <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
#   rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
# gravity_dist_bw_c <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
#   rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
# gravity_dist_bw_d <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
#   rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
# gravity_dist_bw_e <- gravity_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
#   rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)
# 
# 
# # joining all together
# gravity_dist_bw_all <- gravity_dist_bw_b %>%
#   left_join(gravity_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
#   left_join(gravity_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
#   left_join(gravity_dist_bw_e, by = c("block_time", "block", "time", "year"))
# 
# # calculate distance using pythagorean theorem
# gravity_dist_bw_all <- gravity_dist_bw_all %>%
#   mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
#   mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
#   mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
#   mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
#   mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
#   mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))
# 
# # putting together
# gravity_dist_bw_all <- gravity_dist_bw_all %>%
#   select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
# gravity_dist_bw_all <- gravity_dist_bw_all %>%
#   pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
#                names_to = "patch_pair", values_to = "distance")
# gravity_dist_bw_all$time <- as.numeric(gravity_dist_bw_all$time)
# 
# 
# # joining to data 
# gravity_dist_bw_all <- gravity_dist_bw_all %>%
#   left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
#   filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
#   mutate(dispersal_mode = "Gravity")
# 
# gravity_dist_bw_all %>%
#   ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 24) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab("Distance between patch type communities") #+
# # annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)
# 
# 
# 
# 
# 
# ##### wind dispersed convergence/divergence #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
wind_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Wind") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
wind_convergence_jaccard <- wind_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Wind")
# # use PCoA axes
# wind_pcoa_dist_bw_patch <- wind_pcoa_axes %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
#   #filter(!block %in% c("75W", "75E")) %>%
#   mutate(block_time = paste(block, time, sep = "-"))
# 
# # separating by patch replicate
# wind_dist_bw_b <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
#   rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
# wind_dist_bw_c <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
#   rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
# wind_dist_bw_d <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
#   rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
# wind_dist_bw_e <- wind_pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
#   rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)
# 
# 
# # joining all together
# wind_dist_bw_all <- wind_dist_bw_b %>%
#   left_join(wind_dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
#   left_join(wind_dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
#   left_join(wind_dist_bw_e, by = c("block_time", "block", "time", "year"))
# 
# # calculate distance using pythagorean theorem
# wind_dist_bw_all <- wind_dist_bw_all %>%
#   mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
#   mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
#   mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
#   mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
#   mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
#   mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))
# 
# # putting together
# wind_dist_bw_all <- wind_dist_bw_all %>%
#   select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
# wind_dist_bw_all <- wind_dist_bw_all %>%
#   pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
#                names_to = "patch_pair", values_to = "distance")
# wind_dist_bw_all$time <- as.numeric(wind_dist_bw_all$time)
# 
# 
# # joining to data 
# wind_dist_bw_all <- wind_dist_bw_all %>%
#   left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
#   filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
#   mutate(dispersal_mode = "Wind")
# 
# wind_dist_bw_all %>%
#   ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 24) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab("Distance between patch type communities") #+
# # annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)
# 
# 
# 
# ##### convergence all dispersal modes together #####
dispersal_mode_convergence <- rbind(
    animal_convergence_jaccard, gravity_convergence_jaccard, wind_convergence_jaccard
  )
dispersal_mode_convergence$s.time <- as.numeric(scale(dispersal_mode_convergence$time))
animal_convergence_jaccard$s.time <- as.numeric(scale(animal_convergence_jaccard$time))
gravity_convergence_jaccard$s.time <- as.numeric(scale(gravity_convergence_jaccard$time))
wind_convergence_jaccard$s.time <- as.numeric(scale(wind_convergence_jaccard$time))

# models
# linear model
m.dispersal_convergence <- glmmTMB(jaccard ~ patch_pair*s.time*dispersal_mode + (1|block),
                                   data = dispersal_mode_convergence)
# quadratic model
m.dispersal_convergence_quad <- glmmTMB(jaccard ~ patch_pair*dispersal_mode*(s.time )+ patch_pair*I(s.time^2)*dispersal_mode + (1|block),
                                   data = dispersal_mode_convergence)
# not converging, excluding from AIC model list
# null model
m.dispersal_convergence_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                        data = dispersal_mode_convergence)
# AIC comparison
a <- list(m.dispersal_convergence, m.dispersal_convergence_quad, m.dispersal_convergence_null)
aictab(a) # quadratic much better fit


# model checking
summary(m.dispersal_convergence)
plot(simulateResiduals(m.dispersal_convergence)) # some issues
check_model(m.dispersal_convergence) # some issues
performance::r2(m.dispersal_convergence)


### plotting model predictions
# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- dispersal_mode_convergence %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# model predictions
m.dispersal_converge.predict <- ggpredict(m.dispersal_convergence, terms=c("s.time [all]", "patch_pair [all]", "dispersal_mode [all]"), back_transform = T)
m.dispersal_converge.predict <- as.data.frame(m.dispersal_converge.predict)

# plotting
dispersal_convergence_plot <- m.dispersal_converge.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  rename(dispersal_mode = facet) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3.5, alpha = 0.09, data = dispersal_mode_convergence) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  theme_minimal(base_size = 24) +
  facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab("Time since site creation (years)") +
  ylab(expression(atop("Dissimilarity between ", paste("patch type communities"))))
dispersal_convergence_plot

# pdf(file = file.path("plots", "dispersal_convergence.pdf"), width = 14, height = 6)
# dispersal_convergence_plot
# dev.off()


### individual dispersal mode models
m.converge_animal <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                             data = animal_convergence_jaccard)
summary(m.converge_animal)
m.converge_gravity <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                              data = gravity_convergence_jaccard)
summary(m.converge_gravity)
m.converge_wind <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                           data = wind_convergence_jaccard)
summary(m.converge_wind)

## animal dispersed plot
# model predictions
m.converge_animal.predict <- ggpredict(m.converge_animal, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_animal.predict <- as.data.frame(m.converge_animal.predict)

# plotting
animal_convergence_plot <- m.converge_animal.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  #rename(dispersal_mode = facet) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3.5, alpha = 0.09, data = animal_convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  theme_minimal(base_size = 24) +
  #facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab(NULL) +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  theme(legend.position = "none") +
  ylim(c(0, 0.8))
animal_convergence_plot


## gravity dispersed plot
# model predictions
m.converge_gravity.predict <- ggpredict(m.converge_gravity, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_gravity.predict <- as.data.frame(m.converge_gravity.predict)

# plotting
gravity_convergence_plot <- m.converge_gravity.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  #rename(dispersal_mode = facet) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3.5, alpha = 0.09, data = gravity_convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  theme_minimal(base_size = 24) +
  #facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab("Time since site creation (years)") +
  ylab(NULL) +
  theme(legend.position = "none") +
  ylim(c(0, 0.8))
gravity_convergence_plot


## wind dispersed plot
# model predictions
m.converge_wind.predict <- ggpredict(m.converge_wind, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_wind.predict <- as.data.frame(m.converge_wind.predict)

# plotting
wind_convergence_plot <- m.converge_wind.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  #rename(dispersal_mode = facet) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3.5, alpha = 0.09, data = wind_convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  theme_minimal(base_size = 24) +
  #facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab(NULL) +
  ylab(NULL) +
  ylim(c(0, 0.8))
wind_convergence_plot
# 
# dispersal_mode_convergence %>%
#     ggplot(aes(time, jaccard, color = patch_pair, fill = patch_pair)) +
#     geom_point(size = 3, alpha = 0.3) +
#     geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
#     facet_wrap(~dispersal_mode) +
#     theme_minimal(base_size = 24) +
#     scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#     scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#     xlab("Time since site creation (years)") +
#     ylab(expression(atop("Distance between ", paste("patch type communities"))))


# dispersal_mode_convergence <- rbind(
#   wind_dist_bw_all, gravity_dist_bw_all, animal_dist_bw_all
# )
# 
# 
# dispersal_mode_convergence_plot <- dispersal_mode_convergence %>%
#   ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
#   facet_wrap(~dispersal_mode) +
#   theme_minimal(base_size = 24) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab(expression(atop("Distance between ", paste("patch type communities")))) 
# dispersal_mode_convergence_plot

### need to join all the model predictions together than facet!!!!!!!!!
all_dispersal_convergence_plot <- cowplot::plot_grid(animal_convergence_plot, gravity_convergence_plot, wind_convergence_plot, 
                                                     nrow = 1, rel_widths = c(1.1, 1, 2), axis = "t")
all_dispersal_convergence_plot
# 
# # pdf(file = file.path("plots", "dispersal_convergence.pdf"), width = 14, height = 6)
# # dispersal_mode_convergence_plot
# # dev.off()


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
dispersal_mode_segments <- dispersal_mode_segments %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
    ))
dispersal_mode_segments$s.time <- as.numeric(scale(dispersal_mode_segments$time))

# dispersal_mode_segments_plot <- dispersal_mode_segments %>%
#   mutate(patch_type = dplyr::case_when(
#     patch_type %in% c("connected") ~ "Connected",
#     patch_type %in% c("rectangle") ~ "Rectangular",
#     patch_type %in% c("wing") ~ "Winged"
#   )) %>%
#   ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
#   geom_point(size = 4.5, alpha = 0.3) +
#   theme_minimal(base_size = 28) +
#   facet_wrap(~dispersal_mode) +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
#   xlab("Time since site creation (years)") #+
# #annotate("text", x = 20, y=0.47, label = expression(paste('R'^2*' = 0.431')), size=7)
# dispersal_mode_segments_plot


# models
# linear model
m.dispersal_segments <- glmmTMB(distance ~ patch_type*s.time*dispersal_mode + (1|block),
                                   data = dispersal_mode_segments)
# quadratic model
m.dispersal_segments_quad <- glmmTMB(distance ~ patch_type*dispersal_mode*s.time + patch_type*I(s.time^2)*dispersal_mode + (1|block),
                                        data = dispersal_mode_segments)
# not converging, excluding from AIC model list
# null model
m.dispersal_segments_null <- glmmTMB(distance ~ 1 + (1|block),
                                        data = dispersal_mode_segments)
# AIC comparison
a <- list(m.dispersal_segments, m.dispersal_segments_quad, m.dispersal_segments_null)
aictab(a) # quadratic much better fit


# model checking
summary(m.dispersal_segments_quad)
plot(simulateResiduals(m.dispersal_segments_quad)) # some issues
check_model(m.dispersal_segments_quad) # some issues
check_collinearity(m.dispersal_segments_quad)
performance::r2(m.dispersal_segments_quad)


### individual models
# animal
animal_mode_segments <- dispersal_mode_segments %>%
  filter(dispersal_mode == "Animal")
m.animal_segments_quad <- glmmTMB(distance ~ patch_type*s.time + patch_type*I(s.time^2) + (1|block),
                                     data = animal_mode_segments)
summary(m.animal_segments_quad)
animal_segments.posthoc <- emmeans(m.animal_segments_quad, ~patch_type*s.time, at = list(s.time = c(0, 1.11640239)))
pairs(animal_segments.posthoc, simple = "patch_type")

# gravity
gravity_mode_segments <- dispersal_mode_segments %>%
  filter(dispersal_mode == "Gravity")
m.gravity_segments_quad <- glmmTMB(distance ~ patch_type*s.time + patch_type*I(s.time^2) + (1|block),
                                  data = gravity_mode_segments)
summary(m.gravity_segments_quad)
gravity_segments.posthoc <- emmeans(m.gravity_segments_quad, ~patch_type*s.time, at = list(s.time = c(0, 1.11640239)))
pairs(gravity_segments.posthoc, simple = "patch_type")


# wind
wind_mode_segments <- dispersal_mode_segments %>%
  filter(dispersal_mode == "Wind")
m.wind_segments_quad <- glmmTMB(distance ~ patch_type*s.time + patch_type*I(s.time^2) + (1|block),
                                   data = wind_mode_segments)
summary(m.wind_segments_quad)
wind_segments.posthoc <- emmeans(m.wind_segments_quad, ~patch_type*s.time, at = list(s.time = c(0, 1.11640239)))
pairs(wind_segments.posthoc, simple = "patch_type")





### plotting model predictions
# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- dispersal_mode_segments %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# model predictions
m.dispersal_segments.predict <- ggpredict(m.dispersal_segments_quad, terms=c("s.time [all]", "patch_type [all]", "dispersal_mode [all]"), back_transform = T)
m.dispersal_segments.predict <- as.data.frame(m.dispersal_segments.predict)

# plotting
dispersal_segments_plot <- m.dispersal_segments.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  rename(dispersal_mode = facet) %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 3.5, alpha = 0.15, data = dispersal_mode_segments) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  theme_minimal(base_size = 24) +
  facet_wrap(~dispersal_mode) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) 
dispersal_segments_plot


pdf(file = file.path("plots", "dispersal_segments.pdf"), width = 14, height = 6)
dispersal_segments_plot
dev.off()



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
performance::r2(m.direction_dispersal)
Anova(m.direction_dispersal)

m.direction_dispersal.posthoc <- emmeans(m.direction_dispersal, ~ patch_type*time*dispersal_mode)
pairs(m.direction_dispersal.posthoc, simple = c("patch_type", "time"))
pairs(m.direction_dispersal.posthoc)
# trends are primarily driven by wind dispersed species


# prediction plot
m.direction_dispersal.predict <- ggpredict(m.direction_dispersal, terms=c("time [all]", "patch_type [all]", "dispersal_mode [all]"), back_transform = T)
m.direction_dispersal.predict$dispersal_mode <- m.direction_dispersal.predict$facet

dispersal_direction_plot <- m.direction_dispersal.predict %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type), 
              data = directionality_all, alpha = 0.2, size = 7, 
              position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = m.direction_dispersal.predict, width = 0, linewidth = 4.5,  position = position_dodge(width = 0.7)) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 10, data = m.direction_dispersal.predict, 
             position = position_dodge(width = 0.7), colour="black", pch=21, stroke = 2)+ 
  theme_minimal(base_size = 28) +
  facet_wrap(~dispersal_mode) +
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

pdf(file = file.path("plots", "dispersal_direction.pdf"), width = 15.5, height = 7)
dispersal_direction_plot
dev.off()










######## trajectory angles?? ########
# animal
animal_angles <- trajectoryAngles(animal_srs_trajectory)
animal_angles <- animal_angles %>%
  rownames_to_column("uniqueID")
animal_angles <- animal_angles %>%
  pivot_longer(cols = 2:25, names_to = "time_period", values_to = "angle")

# assinging last year of sequence as time
animal_angles <- animal_angles %>%
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

animal_angles <- animal_angles %>%
  separate(uniqueID, into = c("block", "patch_rep", "patch_type")) %>%
  mutate(dispersal_mode = "Animal")

# plotting
animal_angles %>%
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
  theme(plot.margin = margin(1, 1, 1, 2, "cm"))+
  xlab("Time since site creation (years)") +
  #ylab(expression(paste("Trajectory directionality \nbetween consecutive years", " (angle ", theta, ")")))
  ylab("Trajectory directionality \nbetween consecutive years")





# gravity
gravity_angles <- trajectoryAngles(gravity_srs_trajectory)
gravity_angles <- gravity_angles %>%
  rownames_to_column("uniqueID")
gravity_angles <- gravity_angles %>%
  pivot_longer(cols = 2:25, names_to = "time_period", values_to = "angle")

# assinging last year of sequence as time
gravity_angles <- gravity_angles %>%
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

gravity_angles <- gravity_angles %>%
  separate(uniqueID, into = c("block", "patch_rep", "patch_type")) %>%
  mutate(dispersal_mode = "Gravity")

# plotting
gravity_angles %>%
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
  theme(plot.margin = margin(1, 1, 1, 2, "cm"))+
  xlab("Time since site creation (years)") +
  #ylab(expression(paste("Trajectory directionality \nbetween consecutive years", " (angle ", theta, ")")))
  ylab("Trajectory directionality \nbetween consecutive years")





# wind
wind_angles <- trajectoryAngles(wind_srs_trajectory)
wind_angles <- wind_angles %>%
  rownames_to_column("uniqueID")
wind_angles <- wind_angles %>%
  pivot_longer(cols = 2:25, names_to = "time_period", values_to = "angle")

# assinging last year of sequence as time
wind_angles <- wind_angles %>%
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

wind_angles <- wind_angles %>%
  separate(uniqueID, into = c("block", "patch_rep", "patch_type")) %>%
  mutate(dispersal_mode = "Wind")

# plotting
wind_angles %>%
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
  theme(plot.margin = margin(1, 1, 1, 2, "cm"))+
  xlab("Time since site creation (years)") +
  #ylab(expression(paste("Trajectory directionality \nbetween consecutive years", " (angle ", theta, ")")))
  ylab("Trajectory directionality \nbetween consecutive years")



### all together
all_angles <- rbind(
  animal_angles, gravity_angles, wind_angles
)

all_angles %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time_period, angle, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  facet_wrap(~dispersal_mode) +
  theme_minimal(base_size = 28) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme(plot.margin = margin(1, 1, 1, 2, "cm"))+
  xlab("Time since site creation (years)") +
  #ylab(expression(paste("Trajectory directionality \nbetween consecutive years", " (angle ", theta, ")")))
  ylab("Trajectory directionality \nbetween consecutive years")




m1 <- glmmTMB(angle ~ patch_type * time_period + (1|block),
                                 data = gravity_angles)
summary(m1)
plot(simulateResiduals(m1))
performance::r2(m1)
Anova(m1)

m1.posthoc <- emtrends(m1, "patch_type", var = "time_period")
pairs(m1.posthoc)

