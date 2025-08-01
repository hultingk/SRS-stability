### community trajectory analysis 
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

#source(here::here("04b_cta_directionality.R"))

# loading data
srs_data_core <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_core.csv"))

srs_data_core <- srs_data_core %>% # removing experimentally planted species 
  filter(transplant != TRUE) #%>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) #%>%
  #filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data_core %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.numeric(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

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

# centering trajectories
center_srs_trajectory <- centerTrajectories(srs_trajectory)
# calculating dissimilarity between trajectories
srs_traj_dist <- trajectoryDistances(center_srs_trajectory)


# converting to dataframe
srs_traj_dist <- as.matrix(srs_traj_dist)
srs_traj_dist <- as.data.frame(srs_traj_dist)
# pivoting
srs_traj_dist <- srs_traj_dist %>%
  rownames_to_column(var = "traj1") %>%
  pivot_longer(2:41, names_to = "traj2", values_to = "resemblance")

# removing trajectories compared to itself
srs_traj_dist <- srs_traj_dist %>%
  filter(resemblance != 0)
# removing duplicate comparisons
srs_traj_dist <- unique(srs_traj_dist)


# keeping only unique pairs of patches
srs_traj_dist$pair_min <- pmin(srs_traj_dist$traj1, srs_traj_dist$traj2)
srs_traj_dist$pair_max <- pmax(srs_traj_dist$traj1, srs_traj_dist$traj2)
# keep only the first occurrence of each unordered pair - should have 496 pairs
srs_traj_dist <- srs_traj_dist[!duplicated(srs_traj_dist[c("pair_min", "pair_max")]), ]

# wrangling patch name columns
srs_traj_dist <- srs_traj_dist %>%
  separate(traj1, into = c("block1", "patch_rep1", "patch_type1"), sep = "-") %>%
  separate(traj2, into = c("block2", "patch_rep2", "patch_type2"), sep = "-") %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-"))

# removing unneeded columns
srs_traj_dist$pair_min <- NULL
srs_traj_dist$pair_max <- NULL
srs_traj_dist %>%
  count(patch_pair)

#srs_traj_dist <- srs_traj_dist %>%
#  filter(block1 == block2)
# combining patch pairs
srs_traj_dist <- srs_traj_dist %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("connected-rectangle", "rectangle-connected") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing", "wing-connected") ~ "Connected-Winged",
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Rectangular-Winged",
    patch_pair %in% c("center-connected", "connected-center") ~ "Center-Connected",
    patch_pair %in% c("center-rectangle", "rectangle-center") ~ "Center-Rectangular",
    patch_pair %in% c("center-wing", "wing-center") ~ "Center-Winged",
    .default = patch_pair
  )) %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle", "center-center", "connected-connected"))
  

srs_traj_dist %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle", "center-center", "connected-connected")) %>%
  #filter(patch_pair %in% c("Center-Winged", "Center-Rectangular", "Center-Connected")) %>%
  mutate(patch_pair = fct_reorder(patch_pair, resemblance, mean, .na_rm = T)) %>%
  ggplot(aes(patch_pair, resemblance, fill = patch_pair)) +
  geom_boxplot() +
 # scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 60,  hjust=1))


center_comparison <- srs_traj_dist %>%
  #filter(patch_pair %in% c("Connected-Winged", "Connected-Rectangular", "Rectangular-Winged"))
  filter(patch_pair %in% c("Center-Winged", "Center-Rectangular", "Center-Connected"))
srs_traj_dist$patch_pair <- as.factor(srs_traj_dist$patch_pair)
m1 <- glmmTMB(resemblance ~ patch_pair + (1|block1), 
              data = srs_traj_dist,
              family = "gaussian")
summary(m1)
pairs(emmeans(m1, ~patch_pair))
plot(simulateResiduals(m1))
check_model(m1)
# trajectories in connected patches resemble center patches more than center patches compared to other patch types
m1.predict <- ggpredict(m1, terms=c("patch_pair [all]"), back_transform = T)
srs_traj_dist %>% 
  mutate(patch_pair = fct_reorder(patch_pair, resemblance, mean, .na_rm = T)) %>%
  ggplot(aes(x = fct_reorder(patch_pair, resemblance, mean, .na_rm = T), y = resemblance, color = patch_pair)) +
  geom_jitter(aes(x = fct_reorder(patch_pair, resemblance, mean, .na_rm = T), y = resemblance, color = patch_pair), data = srs_traj_dist, 
              alpha = 0.05, size = 3.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 0.5)) +
  geom_point(aes(x = fct_reorder(x, predicted, mean, .na_rm = T), y = predicted, color = x), size = 4.5, data = m1.predict,  position = position_dodge(0.5))+ 
  geom_errorbar(aes(x = fct_reorder(x, predicted, mean, .na_rm = T), y = predicted, ymin = conf.low, ymax = conf.high, color = x),data = m1.predict, width = 0, linewidth = 2,  position = position_dodge(0.5)) +
  theme_classic() +
  #geom_text(data = m1.stat.test, aes(x = ptype, y = height, label = significance), size = 6) +
  labs(title = NULL,
       x = NULL,
       y = "Trajectory dissimilarity between patch pairs") +
 # scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Pair") +
 # scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Pair") +
  theme(axis.text = element_text(size = 22)) +
  theme(axis.title = element_text(size = 26)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 60,  hjust=1))




###### broken up by first and second half of succession #####

#### FIRST HALF
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

# centering trajectories
center_srs_trajectory_1_12 <- centerTrajectories(srs_trajectory_1_12)
# calculating dissimilarity between trajectories
srs_traj_dist_1_12 <- trajectoryDistances(center_srs_trajectory_1_12)


# converting to dataframe
srs_traj_dist_1_12 <- as.matrix(srs_traj_dist_1_12)
srs_traj_dist_1_12 <- as.data.frame(srs_traj_dist_1_12)
# pivoting
srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  rownames_to_column(var = "traj1") %>%
  pivot_longer(2:41, names_to = "traj2", values_to = "resemblance")

# removing trajectories compared to itself
srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  filter(resemblance != 0)
# removing duplicate comparisons
srs_traj_dist_1_12 <- unique(srs_traj_dist_1_12)

# keeping only unique pairs of patches
srs_traj_dist_1_12$pair_min <- pmin(srs_traj_dist_1_12$traj1, srs_traj_dist_1_12$traj2)
srs_traj_dist_1_12$pair_max <- pmax(srs_traj_dist_1_12$traj1, srs_traj_dist_1_12$traj2)
# keep only the first occurrence of each unordered pair - should have 496 pairs
srs_traj_dist_1_12 <- srs_traj_dist_1_12[!duplicated(srs_traj_dist_1_12[c("pair_min", "pair_max")]), ]

# wrangling patch name columns
srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  separate(traj1, into = c("block1", "patch_rep1", "patch_type1"), sep = "-") %>%
  separate(traj2, into = c("block2", "patch_rep2", "patch_type2"), sep = "-") %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-"))

# removing unneeded columns
srs_traj_dist_1_12$pair_min <- NULL
srs_traj_dist_1_12$pair_max <- NULL
srs_traj_dist_1_12 %>%
  count(patch_pair)

srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  filter(block1 == block2)
# combining patch pairs
srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("connected-rectangle", "rectangle-connected") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing", "wing-connected") ~ "Connected-Winged",
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Rectangular-Winged",
    patch_pair %in% c("center-connected", "connected-center") ~ "Center-Connected",
    patch_pair %in% c("center-rectangle", "rectangle-center") ~ "Center-Rectangular",
    patch_pair %in% c("center-wing", "wing-center") ~ "Center-Winged",
    .default = patch_pair
  ))

srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  dplyr::select(block1, resemblance, patch_pair)




#### SECOND HALF
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

# centering trajectories
center_srs_trajectory_13_24 <- centerTrajectories(srs_trajectory_13_24)
# calculating dissimilarity between trajectories
srs_traj_dist_13_24 <- trajectoryDistances(center_srs_trajectory_13_24)


# converting to dataframe
srs_traj_dist_13_24 <- as.matrix(srs_traj_dist_13_24)
srs_traj_dist_13_24 <- as.data.frame(srs_traj_dist_13_24)
# pivoting
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  rownames_to_column(var = "traj1") %>%
  pivot_longer(2:33, names_to = "traj2", values_to = "resemblance")

# removing trajectories compared to itself
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  filter(resemblance != 0)
# removing duplicate comparisons
srs_traj_dist_13_24 <- unique(srs_traj_dist_13_24)

# keeping only unique pairs of patches
srs_traj_dist_13_24$pair_min <- pmin(srs_traj_dist_13_24$traj1, srs_traj_dist_13_24$traj2)
srs_traj_dist_13_24$pair_max <- pmax(srs_traj_dist_13_24$traj1, srs_traj_dist_13_24$traj2)
# keep only the first occurrence of each unordered pair - should have 496 pairs
srs_traj_dist_13_24 <- srs_traj_dist_13_24[!duplicated(srs_traj_dist_13_24[c("pair_min", "pair_max")]), ]

# wrangling patch name columns
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  separate(traj1, into = c("block1", "patch_rep1", "patch_type1"), sep = "-") %>%
  separate(traj2, into = c("block2", "patch_rep2", "patch_type2"), sep = "-") %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-"))

# removing unneeded columns
srs_traj_dist_13_24$pair_min <- NULL
srs_traj_dist_13_24$pair_max <- NULL
srs_traj_dist_13_24 %>%
  count(patch_pair)

srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  filter(block1 == block2)
# combining patch pairs
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("connected-rectangle", "rectangle-connected") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing", "wing-connected") ~ "Connected-Winged",
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Rectangular-Winged",
    patch_pair %in% c("center-connected", "connected-center") ~ "Center-Connected",
    patch_pair %in% c("center-rectangle", "rectangle-center") ~ "Center-Rectangular",
    patch_pair %in% c("center-wing", "wing-center") ~ "Center-Winged",
    .default = patch_pair
  ))
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  dplyr::select(block1, resemblance, patch_pair)



### combining 
srs_traj_dist_1_12 <- srs_traj_dist_1_12 %>%
  mutate(decade = "Years 1-12")
srs_traj_dist_13_24 <- srs_traj_dist_13_24 %>%
  mutate(decade = "Years 13-24")


srs_dist_decade <- rbind(
  srs_traj_dist_1_12,
  srs_traj_dist_13_24
)


srs_dist_decade <- srs_dist_decade %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) 
srs_dist_decade %>%
  mutate(patch_pair = fct_reorder(patch_pair, resemblance, mean, .na_rm = T)) %>%
  ggplot() +
  geom_boxplot(aes(patch_pair, resemblance, fill = patch_pair)) +
  facet_wrap(~decade) +
  #scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
 # scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60,  hjust=1))

m2 <- glmmTMB(resemblance ~ patch_pair*decade + (1|block1), 
              data = srs_dist_decade)
summary(m2)
plot(simulateResiduals(m2))
pairs(emmeans(m2, ~patch_pair*decade), simple = "patch_pair")



















