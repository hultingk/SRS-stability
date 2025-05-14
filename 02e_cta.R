### community trajectory analysis
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.numeric(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))


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

# plotting segment lengths 
segment_lengths_plot <- segment_lengths %>%
  filter(!block %in% c("75E", "75W")) %>%
  ggplot(aes(time, distance, color = patch_type)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal(base_size = 20) +
  geom_smooth() +
  scale_color_brewer(palette = "Set2", name = "Patch Type") +
  ylab("Trajectory distance between consecutive surveys") +
  xlab("Time since experiment")
segment_lengths_plot

pdf(file = "segment_lengths.pdf", width = 11, height = 8)
segment_lengths_plot
dev.off()


# trajectory directionality
segment_direction <- trajectoryDirectionality(srs_trajectory, nperm = 999)
segment_direction <- data.frame(segment_direction)
segment_direction <- segment_direction %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) 

###### directionality in first decade since patch creation #####
sp_info_1_10 <- srs_data_wider %>%
  filter(time < 11)

patch_info_1_10 <- sp_info_1_10 %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info_1_10 <- sp_info_1_10 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_1_10 <- vegdist(sp_info_1_10, method = "jaccard")

# defining trajectories
srs_trajectory_1_10 <- defineTrajectories(jaccard_dist_1_10, sites = patch_info_1_10$unique_id, surveys = patch_info_1_10$time)

# directionality
segment_direction_1_10 <- trajectoryDirectionality(srs_trajectory_1_10)
segment_direction_1_10 <- data.frame(segment_direction_1_10)
segment_direction_1_10 <- segment_direction_1_10 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-10") %>%
  rename(directionality = segment_direction_1_10)





###### directionality in SECOND decade since patch creation #####
sp_info_11_24 <- srs_data_wider %>%
  filter(time > 11)

patch_info_11_24 <- sp_info_11_24 %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info_11_24 <- sp_info_11_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_11_24 <- vegdist(sp_info_11_24, method = "jaccard")

# defining trajectories
srs_trajectory_11_24 <- defineTrajectories(jaccard_dist_11_24, sites = patch_info_11_24$unique_id, surveys = patch_info_11_24$time)

# directionality
segment_direction_11_24 <- trajectoryDirectionality(srs_trajectory_11_24)
segment_direction_11_24 <- data.frame(segment_direction_11_24)
segment_direction_11_24 <- segment_direction_11_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 11-24") %>%
  rename(directionality = segment_direction_11_24)


#### putting all together
segment_direction_all <- rbind(
  segment_direction_1_10,
  segment_direction_11_24
)

segment_direction_all %>%
  ggplot(aes(time, directionality)) +
  #geom_point() +
  geom_boxplot(aes(fill = patch_type)) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2", name = "Patch Type") +
  ylab("Trajectory directionality") +
  xlab("Decade")


# model
m.direction <- glmmTMB(directionality ~ time*patch_type + (1|block/patch),
                       data = segment_direction_all)
summary(m.direction)
