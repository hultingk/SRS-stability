### community trajectory analysis 
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) #%>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  #filter(patch_type != "center")


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
# combining patch pairs
srs_traj_dist <- srs_traj_dist %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("connected-rectangle", "rectangle-connected") ~ "connected-rectangle",
    patch_pair %in% c("connected-wing", "wing-connected") ~ "connected-wing",
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "rectangle-wing",
    patch_pair %in% c("center-connected", "connected-center") ~ "center-connected",
    patch_pair %in% c("center-rectangle", "rectangle-center") ~ "center-rectangle",
    patch_pair %in% c("center-wing", "wing-center") ~ "center-wing",
    .default = patch_pair
  ))

srs_traj_dist %>%
  ggplot(aes(patch_pair, resemblance, fill = patch_pair)) +
  geom_boxplot()

center_comparison <- srs_traj_dist %>%
  filter(patch_pair %in% c("center-connected", "center-rectangle", "center-wing"))

m1 <- glmmTMB(resemblance ~ patch_pair, 
              data = center_comparison)
summary(m1)
pairs(emmeans(m1, ~patch_pair))
# trajectories in connected patches resemble center patches more than center patches compared to other patch types
