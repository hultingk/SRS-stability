### community trajectory analysis - directionality
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

source(here::here("04a_cta_segments.R"))


##### directionality broken into time periods ####
###### directionality in first 12 years #####
sp_info_1_12 <- srs_data_wider %>%
  filter(time <= 12)

patch_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

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
  select(unique_id, time, year)

# species matrix
sp_info_13_24 <- sp_info_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))

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
  xlab("Time period")
segment_direction_plot

# pdf(file = file.path("plots", "segment_direction.pdf"), width = 12, height = 8)
# segment_direction_plot
# dev.off()


# model
m.direction <- glmmTMB(directionality ~ patch_type * time + (1|block/patch),
                       data = segment_direction_all)
summary(m.direction)
plot(simulateResiduals(m.direction))



#### directionality over time ####
# 1. for every three consecutive surveys (e.g., 1-3, 2-4, 3-5)
# ---- calculate jaccard distance
# ---- define trajectories
# ---- calculate directionality
# ---- use last year of three year period in analysis as time
srs_data_wider %>%
  count(time) %>%
  View()

calculate_directionality <- function(df) {
  patch_info <- df %>% 
    arrange(unique_id, time) %>%
    select(unique_id, time, year)
  
  # species matrix
  sp_info <- df %>%
    arrange(unique_id, time) %>%
    mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
    column_to_rownames("unique_id_year") %>%
    select(!c("unique_id", "time", "year"))
  
  jaccard <- vegdist(sp_info, method = "jaccard")
  
  trajectory <- defineTrajectories(jaccard, sites = patch_info$unique_id, surveys = patch_info$time)
  
  segment_direction <- trajectoryDirectionality(trajectory)
  segment_direction <- data.frame(segment_direction)
  segment_direction <- segment_direction %>%
    rownames_to_column("unique_id") %>%
    separate(unique_id, into = c("block", "patch", "patch_type")) %>%
    #mutate(time = "Year 13-24") %>%
    rename(directionality = segment_direction)
  
  return(segment_direction)
}

patch_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))








calculate_trajectory_directionality <- function(comm_matrix, time_vector) {
  # comm_matrix: site × species matrix, with rows as observations across time
  # time_vector: vector of time points corresponding to each row in comm_matrix
  
  directionality_results <- data.frame(
    end_year = numeric(),
    directionality = numeric(),
    stringsAsFactors = FALSE
  )
  
  patch_info <- comm_matrix %>%
    arrange(unique_id, time) %>%
    select(unique_id, time, year)
  
  sp_info <- comm_matrix %>%
    arrange(unique_id, time) %>%
    mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
    column_to_rownames("unique_id_year") %>%
    select(!c("unique_id", "time", "year"))
  
  n_years <- max(patch_info$time)
  
  for (i in 1:(n_years - 2)) {
    # Subset for 3 consecutive years
    sub_matrix <- sp_info[i:(i+2), , drop = FALSE]
    sub_patch <- patch_info[i:(i+2), ]
    time_value <- sub_patch$time
    
    # Skip if any NA in years or community data
    if (any(is.na(time_value)) || any(is.na(sub_matrix))) next
    
    if (length(unique(sub_patch)) == 3 && all(diff(sub_patch$time) > 0)) {
       # Step 1: Calculate Jaccard distances
    jaccard_dist <- vegdist(sub_matrix, method = "jaccard", binary = TRUE)
    
    # Step 2: Define trajectory
    traj <- defineTrajectories(jaccard_dist, sites = sub_patch$unique_id, surveys = sub_patch$time)
    
    # Step 3: Calculate directionality
    dir_val <- trajectoryDirectionality(traj)
    
    # Step 4: Store the result with the last year in the 3-year window
    directionality_results <- rbind(directionality_results, data.frame(
      end_year = max(time_value),
      directionality = dir_val
    ))
    }
  }
  
  return(directionality_results)
}

srs_direction_split <- srs_data_wider %>%
  #separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  #filter(!block %in% c("75E", "75W")) %>%
  #select(!c("block", "patch_rep", "patch_type")) %>%
  group_by(unique_id) %>%
  group_split()

directionality <- srs_direction_split %>%
  lapply(calculate_trajectory_directionality)


directionality <- directionality %>%
  bind_rows() # putting together into a dataframe

directionality <- directionality %>%
  rownames_to_column("unique_ID") %>%
  separate(unique_ID, into = c("block", "patch_rep", "patch_type_time"), sep = "-") %>%
  select(!patch_type_time)


patch_type <- srs_data %>%
  count(block, patch, patch_type) %>%
  select(!n)

directionality <- directionality %>%
  left_join(patch_type, by = c("block" = "block",
                               "patch_rep" = "patch"))


directionality %>%
  ggplot(aes(end_year, directionality, color = patch_type, fill = patch_type)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_minimal() +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")

m1 <- glmmTMB(directionality ~ patch_type + end_year + (1|block),
              data = directionality)
summary(m1)


hist(directionality$directionality)






###### trajectory angles instead? ######

