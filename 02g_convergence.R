librarian::shelf(tidyverse, vegan, ape, BiodiversityR)

source(here::here("02f_pcoa_permanova.R"))
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  filter(!time %in% c("0", "4")) %>%
  dplyr::count(unique_id, time, year, sppcode, soil_moisture, year_since_fire) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
srs_data_wider$time <- as.factor(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)
srs_data_wider$year_since_fire <- as.numeric(srs_data_wider$year_since_fire)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year, soil_moisture, year_since_fire)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year", "soil_moisture", "year_since_fire"))






patch_info <- patch_info %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  mutate(patch_time = paste(patch_type, time, sep = "-"))


# calculate distance to centroid within a patch type
bd_patch_list <- betadisper(vegdist(sp_info, method = "jaccard"), patch_info$patch_time, type = "centroid")   

# make a dataframe
bd_patch <- as.data.frame(bd_patch_list$group.distances, row.names = rownames(bd_patch_list$group))
bd_patch <- bd_patch %>%
  rownames_to_column("patch_time") %>%
  separate(patch_time, into = c("patch_type", "time"), sep = "-")

# renaming column
bd_patch$distance_centroids <- bd_patch$`bd_patch_list$group.distances`
bd_patch$`bd_patch_list$group.distances` <- NULL
bd_patch$time <- as.numeric(bd_patch$time)

##### trying with individual patch distances to group centroid
# make a dataframe
dist_patch <- as.data.frame(bd_patch_list$distances)
dist_patch <- dist_patch %>%
  rownames_to_column("patch_time") %>%
  separate(patch_time, into = c("block", "patch_rep", "patch_type", "time", "year"), sep = "-")

# renaming column
dist_patch$distance_centroids <- dist_patch$`bd_patch_list$distances`
dist_patch$`bd_patch_list$distances` <- NULL
dist_patch$time <- as.numeric(dist_patch$time)




#############################################################
######## convergence/divergence WITHIN a patch type #########
#### plot distance to patch type centroid (convergence within a patch type)
bd_patch %>%
  ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
  theme_minimal(base_size = 16) +
  xlab("Time since site creation (years)")+
  ylab("Average distance to patch type centroid") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette ="Set2")

# alternative - broken up by block, not an average
dist_patch %>%
  ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.2) +
  theme_minimal(base_size = 16) +
  xlab("Time since site creation (years)")+
  ylab("Distance to patch type centroid") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette ="Set2")


# models
m_centroid <- glmmTMB(distance_centroids ~ patch_type + time + (1|block/patch_rep),
                      data = dist_patch)
summary(m_centroid)






#############################################################
####### calculating distance BETWEEN patch type centroids #######
# pulling out centroids
centroids <- as.data.frame(bd_patch_list[["centroids"]])
centroids <- centroids[,1:2] # only keeping first 2 axes
centroids <- centroids %>%
  rownames_to_column("patch_time") %>%
  separate(patch_time, into = c("patch_type", "time"), sep = "-", remove = F)

# separating by patch type
centroids_c <- centroids %>% filter(patch_type == "connected")
centroids_w <- centroids %>% filter(patch_type == "wing")
centroids_r <- centroids %>% filter(patch_type == "rectangle")

# putting all together
centroids_all <- merge(centroids_c, centroids_w, by = "time", all.x = TRUE)
centroids_all <- merge(centroids_all, centroids_r, by = "time", all.x = TRUE)
centroids_all <- centroids_all %>% # renaming columns for clarity
  rename(c_PCoA1 = PCoA1.x, c_PCoA2 = PCoA2.x,
         w_PCoA1 = PCoA1.y, w_PCoA2 = PCoA2.y,
         r_PCoA1 = PCoA1, r_PCoA2 = PCoA2)

# calculate distance between centroids in pairs of patch types
centroids_all <- centroids_all %>%
  mutate(dist_c_w = sqrt((w_PCoA1 - c_PCoA1)^2 + (w_PCoA2 - c_PCoA2)^2)) %>%
  mutate(dist_c_r = sqrt((r_PCoA1 - c_PCoA1)^2 + (r_PCoA2 - c_PCoA2)^2)) %>%
  mutate(dist_r_w = sqrt((w_PCoA1 - r_PCoA1)^2 + (w_PCoA2 - r_PCoA2)^2))

# putting together
bw_group_dist <- centroids_all %>%
  select(time, dist_c_w, dist_c_r, dist_r_w)
bw_group_dist <- bw_group_dist %>%
  pivot_longer(cols = c("dist_c_w", "dist_c_r", "dist_r_w"), names_to = "patch_pair", values_to = "distance")
bw_group_dist$time <- as.numeric(bw_group_dist$time)

# plotting
bw_group_dist %>%
  ggplot(aes(time, distance, color = patch_pair, fill = patch_pair)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
  theme_minimal(base_size = 16) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette ="Set2") +
  xlab("Time since site creation (years)") +
  ylab("Average distance between patch type centroids")

# models
m_centroid_pt <- glmmTMB(distance ~ time * patch_pair,
                         data = bw_group_dist)
summary(m_centroid_pt)


