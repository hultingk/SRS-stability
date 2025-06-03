librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB, AICcmodavg, DHARMa)

source(here::here("02f_pcoa_permanova.R"))
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(!block %in% c("75W", "75E")) %>%
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
conv_within_patch_plot <- dist_patch %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 28) +
  xlab("Time since site creation (years)")+
  ylab("Distance to patch type centroid") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
conv_within_patch_plot



# pdf(file = file.path("plots", "NO75_conv_within_patch_plot.pdf"), width = 10, height =8)
# conv_within_patch_plot
# dev.off()


# models
m_centroid <- glmmTMB(distance_centroids ~ patch_type + time + (1|block/patch_rep),
                      data = dist_patch)
summary(m_centroid)

m_centroid_quad <- glmmTMB(distance_centroids ~ patch_type + time + I(time^2) + (1|block/patch_rep),
                      data = dist_patch)
summary(m_centroid_quad)

m_centroid_null <- glmmTMB(distance_centroids ~ 1,
                           data = dist_patch)
summary(m_centroid_null)

# AIC comparison
a <- list(m_centroid, m_centroid_quad, m_centroid_null)
aictab(a) # quadratic much better fit




#############################################################
####### calculating distance BETWEEN patch type centroids #######
# pulling out centroids
# centroids <- as.data.frame(bd_patch_list[["centroids"]])
# centroids <- centroids[,1:2] # only keeping first 2 axes
# centroids <- centroids %>%
#   rownames_to_column("patch_time") %>%
#   separate(patch_time, into = c("patch_type", "time"), sep = "-", remove = F)
# 
# # separating by patch type
# centroids_c <- centroids %>% filter(patch_type == "connected")
# centroids_w <- centroids %>% filter(patch_type == "wing")
# centroids_r <- centroids %>% filter(patch_type == "rectangle")
# 
# # putting all together
# centroids_all <- merge(centroids_c, centroids_w, by = "time", all.x = TRUE)
# centroids_all <- merge(centroids_all, centroids_r, by = "time", all.x = TRUE)
# centroids_all <- centroids_all %>% # renaming columns for clarity
#   rename(c_PCoA1 = PCoA1.x, c_PCoA2 = PCoA2.x,
#          w_PCoA1 = PCoA1.y, w_PCoA2 = PCoA2.y,
#          r_PCoA1 = PCoA1, r_PCoA2 = PCoA2)
# 
# # calculate distance between centroids in pairs of patch types
# centroids_all <- centroids_all %>%
#   mutate(dist_c_w = sqrt((w_PCoA1 - c_PCoA1)^2 + (w_PCoA2 - c_PCoA2)^2)) %>%
#   mutate(dist_c_r = sqrt((r_PCoA1 - c_PCoA1)^2 + (r_PCoA2 - c_PCoA2)^2)) %>%
#   mutate(dist_r_w = sqrt((w_PCoA1 - r_PCoA1)^2 + (w_PCoA2 - r_PCoA2)^2))
# 
# # putting together
# bw_group_dist <- centroids_all %>%
#   select(time, dist_c_w, dist_c_r, dist_r_w)
# bw_group_dist <- bw_group_dist %>%
#   pivot_longer(cols = c("dist_c_w", "dist_c_r", "dist_r_w"), names_to = "patch_pair", values_to = "distance")
# bw_group_dist$time <- as.numeric(bw_group_dist$time)
# 
# # plotting
# bw_group_dist %>%
#   ggplot(aes(time, distance, color = patch_pair, fill = patch_pair)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
#   theme_minimal(base_size = 16) +
#   scale_fill_brewer(palette = "Set2") +
#   scale_color_brewer(palette ="Set2") +
#   xlab("Time since site creation (years)") +
#   ylab("Average distance between patch type centroids")
# 
# # models
# m_centroid_pt <- glmmTMB(distance ~ time * patch_pair,
#                          data = bw_group_dist)
# summary(m_centroid_pt)
# 
# 
#### calculate distance between patch types within a block for every year
# use PCoA axes
pcoa_dist_bw_patch <- pcoa_axes %>%
  select(!c("soil_moisture", "year_since_fire")) %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
 # filter(!block %in% c("75W", "75E")) %>%
  mutate(block_time = paste(block, time, sep = "-"))


# separating by patch replicate
dist_bw_b <- pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
  rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
dist_bw_c <- pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
  rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
dist_bw_d <- pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
  rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
dist_bw_e <- pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
  rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)

# joining all together
dist_bw_all <- dist_bw_b %>%
  left_join(dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
  left_join(dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
  left_join(dist_bw_e, by = c("block_time", "block", "time", "year"))


dist_bw_all <- dist_bw_all %>%
  mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
  mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
  mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
  mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
  mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
  mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))


# putting together
dist_bw_all <- dist_bw_all %>%
  select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
dist_bw_all <- dist_bw_all %>%
  pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
               names_to = "patch_pair", values_to = "distance")
dist_bw_all$time <- as.numeric(dist_bw_all$time)

# reading in key to pairs of patches
patch_pair_ID <- read.csv(file = file.path("data", "L2_summarized", "patch_pair_ID.csv"))
# joining to data 
dist_bw_all <- dist_bw_all %>%
  left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
  filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular"))
  

dist_bw_all %>%
  count(block)
conv_bw_patch_plot <- dist_bw_all %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Distance between patch type communities") 
conv_bw_patch_plot

pdf(file = file.path("plots", "conv_bw_patch_plot.pdf"), width = 10, height =8)
conv_bw_patch_plot
dev.off()

# models
m_centroid_pt <- glmmTMB(distance ~ patch_pair_ID * time + (1|block),
                         data = dist_bw_all)
summary(m_centroid_pt)
plot(simulateResiduals(m_centroid_pt))

m_centroid_pt_quad <- glmmTMB(distance ~ patch_pair_ID + time + I(time^2) + (1|block),
                         data = dist_bw_all)
summary(m_centroid_pt_quad)

m_centroid_pt_null <- glmmTMB(distance ~ 1,
                              data = dist_bw_all)
summary(m_centroid_pt_null)

# AIC comparison
a <- list(m_centroid_pt, m_centroid_pt_quad, m_centroid_pt_null)
aictab(a) # quadratic much better fit


# posthoc
m_centroid_pt_posthoc <- emmeans(m_centroid_pt, ~ patch_pair_ID)
pairs(m_centroid_pt_posthoc)

