librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB, AICcmodavg, DHARMa, emmeans)

# loading data
srs_data_core <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_core.csv"))

srs_data_core <- srs_data_core %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(!block %in% c("75W", "75E"))

# pivot to wider format
srs_core_wider <- srs_data_core %>%
  filter(!time %in% c("0", "4")) %>%
  dplyr::count(unique_id, time, year, sppcode, soil_moisture, year_since_fire) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
srs_core_wider$time <- as.factor(srs_core_wider$time)
srs_core_wider$unique_id <- as.factor(srs_core_wider$unique_id)
srs_core_wider$year <- as.factor(srs_core_wider$year)
srs_core_wider$year_since_fire <- as.numeric(srs_core_wider$year_since_fire)

# patch data
patch_core_info <- srs_core_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year, soil_moisture, year_since_fire)

patch_core_info <- patch_core_info %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  mutate(patch_time = paste(patch_type, time, sep = "-"))

# species matrix
sp_core_info <- srs_core_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year", "soil_moisture", "year_since_fire"))


# calculate distance to centroid within a patch type
bd_core_list <- betadisper(vegdist(sp_core_info, method = "jaccard"), patch_core_info$patch_time, type = "centroid")   

# make a dataframe
bd_core <- as.data.frame(bd_core_list$group.distances, row.names = rownames(bd_core_list$group))
bd_core <- bd_core %>%
  rownames_to_column("patch_time") %>%
  separate(patch_time, into = c("patch_type", "time"), sep = "-")

# renaming column
bd_core$distance_centroids <- bd_core$`bd_core_list$group.distances`
bd_core$`bd_core_list$group.distances` <- NULL
bd_core$time <- as.numeric(bd_core$time)


##### trying with individual patch distances to group centroid
# make a dataframe
dist_core_patch <- as.data.frame(bd_core_list$distances)
dist_core_patch <- dist_core_patch %>%
  rownames_to_column("patch_time") %>%
  separate(patch_time, into = c("block", "patch_rep", "patch_type", "time", "year"), sep = "-")

# renaming column
dist_core_patch$distance_centroids <- dist_core_patch$`bd_core_list$distances`
dist_core_patch$`bd_core_list$distances` <- NULL
dist_core_patch$time <- as.numeric(dist_core_patch$time)



#############################################################
######## convergence/divergence WITHIN a patch type #########
#### plot distance to patch type centroid (convergence within a patch type)
bd_core %>%
  ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
  theme_minimal(base_size = 16) +
  xlab("Time since site creation (years)")+
  ylab("Average distance to patch type centroid") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette ="Set2")

# alternative - broken up by block, not an average
conv_within_core_plot <- dist_core_patch %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged",
    patch_type %in% c("center") ~ "Center",
  )) %>%
  ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 28) +
  xlab("Time since site creation (years)")+
  ylab("Distance to patch type centroid") +
  scale_fill_manual(values = c("#004D40", "#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#004D40", "#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
conv_within_core_plot


# models
m_centroid_core <- glmmTMB(distance_centroids ~ patch_type + time + (1|block/patch_rep),
                      data = dist_core_patch)
summary(m_centroid_core)

m_centroid_core_quad <- glmmTMB(distance_centroids ~ patch_type + time + I(time^2) + (1|block/patch_rep),
                           data = dist_core_patch)
summary(m_centroid_core_quad)

m_centroid_core_null <- glmmTMB(distance_centroids ~ 1,
                           data = dist_core_patch)
summary(m_centroid_core_null)

# AIC comparison
a <- list(m_centroid_core, m_centroid_core_quad, m_centroid_core_null)
aictab(a) # quadratic much better fit




##### Convergence BETWEEEN patch types ####
# Jaccard distance matrix for pcoa
jaccard_dist_core <- vegdist(sp_core_info, method = "jaccard")
jaccard_dist_core <- as.matrix(jaccard_dist_core)
jaccard_dist_core_df <- as.data.frame(jaccard_dist_core)

jaccard_dist_core_df <- cbind(patch_core_info, jaccard_dist_core_df) # merge with patch info

# pcoa
pcoa_core <- pcoa(jaccard_dist_core)
pcoa_core_cmd <- cmdscale(jaccard_dist_core, eig=TRUE, add=FALSE)  ## Another way of doing the pcoa that gives same results#


# add axes to patch and time
pcoa_axes_core <- pcoa_core$vectors[,c(1,2)]
pcoa_axes_core <- cbind(patch_core_info, pcoa_axes_core)


#### calculate distance between patch types within a block for every year
# use PCoA axes
pcoa_dist_bw_patch_core <- pcoa_axes_core %>%
  select(!c("soil_moisture", "year_since_fire")) %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
  # filter(!block %in% c("75W", "75E")) %>%
  mutate(block_time = paste(block, time, sep = "-"))


# separating by patch replicate
dist_bw_a_core <- pcoa_dist_bw_patch_core %>% filter(patch_rep == "A") %>%
  rename(PCoA1.A = Axis.1, PCoA2.A = Axis.2)
dist_bw_b_core <- pcoa_dist_bw_patch_core %>% filter(patch_rep == "B") %>%
  rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
dist_bw_c_core <- pcoa_dist_bw_patch_core %>% filter(patch_rep == "C") %>%
  rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
dist_bw_d_core <- pcoa_dist_bw_patch_core %>% filter(patch_rep == "D") %>%
  rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
dist_bw_e_core <- pcoa_dist_bw_patch_core %>% filter(patch_rep == "E") %>%
  rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)

# joining all together
dist_bw_core <- dist_bw_b_core %>%
  left_join(dist_bw_a_core, by = c("block_time", "block", "time", "year")) %>%
  left_join(dist_bw_c_core, by = c("block_time", "block", "time", "year")) %>%
  left_join(dist_bw_d_core, by = c("block_time", "block", "time", "year")) %>%
  left_join(dist_bw_e_core, by = c("block_time", "block", "time", "year"))


# calculate distance using pythagorean theorem
dist_bw_core <- dist_bw_core %>%
  mutate(dist_a_b = sqrt((PCoA1.B - PCoA1.A)^2 + (PCoA2.B - PCoA2.A)^2)) %>%
  mutate(dist_a_c = sqrt((PCoA1.C - PCoA1.A)^2 + (PCoA2.C - PCoA2.A)^2)) %>%
  mutate(dist_a_d = sqrt((PCoA1.D - PCoA1.A)^2 + (PCoA2.D - PCoA2.A)^2)) %>%
  mutate(dist_a_e = sqrt((PCoA1.E - PCoA1.A)^2 + (PCoA2.E - PCoA2.A)^2)) %>%
  mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
  mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
  mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
  mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
  mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
  mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))


# putting together
dist_bw_core <- dist_bw_core %>%
  select(block, time, dist_a_b, dist_a_c, dist_a_d, dist_a_e, dist_b_c, dist_b_d, 
         dist_b_e, dist_c_d, dist_c_e, dist_d_e)


dist_bw_core <- dist_bw_core %>%
  pivot_longer(cols = c("dist_a_b", "dist_a_c", "dist_a_d", "dist_a_e", "dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
               names_to = "patch_pair", values_to = "distance")
dist_bw_core$time <- as.numeric(dist_bw_core$time)

# reading in key to pairs of patches
patch_pair_ID <- read.csv(file = file.path("data", "L2_summarized", "patch_pair_ID.csv"))

# joining to data 
dist_bw_core <- dist_bw_core %>%
  left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
  filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  filter(!patch_pair_ID %in% c("Connected-Rectangular", "Connected-Winged", "Winged-Rectangular"))

# plot
conv_bw_patch_core_plot <- dist_bw_core %>%
  ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm",  alpha = 0.5, linewidth = 2) +
  theme_minimal(base_size = 24) +
  #scale_fill_manual(values = c("#004D40", "#5389A4", "#CC6677", "#DCB254", "#E66100", "#332288"), name = "Patch Type") +
  #scale_color_manual(values = c("#004D40", "#5389A4", "#CC6677", "#DCB254", "#E66100", "#332288"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Time since site creation (years)") +
  ylab("Distance between patch type communities") 
conv_bw_patch_core_plot


# models
# linear
m_centroid_core_pt <- glmmTMB(distance ~ patch_pair_ID + time + (1|block),
                         data = dist_bw_core)
summary(m_centroid_core_pt)
plot(simulateResiduals(m_centroid_core_pt))

# quadratic
m_centroid_core_pt_quad <- glmmTMB(distance ~ patch_pair_ID + time + I(time^2) + (1|block),
                              data = dist_bw_core)
summary(m_centroid_core_pt_quad)
plot(simulateResiduals(m_centroid_core_pt_quad))

# null
m_centroid_core_pt_null <- glmmTMB(distance ~ 1,
                              data = dist_bw_core)
summary(m_centroid_core_pt_null)

# AIC comparison
a <- list(m_centroid_core_pt, m_centroid_core_pt_quad, m_centroid_core_pt_null)
aictab(a) # linear is best fit

# posthoc
core_dist_pt_posthoc <- emmeans(m_centroid_core_pt, ~patch_pair_ID)
pairs(core_dist_pt_posthoc)








