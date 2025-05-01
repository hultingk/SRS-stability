# loading libraries
librarian::shelf(tidyverse, vegan) # Install missing packages and load needed libraries

setwd("~/Documents/MSU/SRS community stability/SRS-stability")
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)

#### null model trial ####
srs_sub <- srs_data %>%
  filter(unique_id == "53S-B-connected") %>%
  dplyr::count(unique_id, time, sppcode)

srs_sub_wider <- srs_sub %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
  dplyr::select(!c("unique_id")) %>% # remove unique ID column
  arrange(time) %>% # Ensure years are sorted properly
  column_to_rownames("time") #convert years to rownames



#### TRIAL ####
custom_raup_crick <- function(comm_matrix, reps = 999) {
  comm_matrix <- as.matrix(comm_matrix)
  comm_matrix[comm_matrix > 0] <- 1  # Presence/absence
  
  n_sites <- nrow(comm_matrix)
  species_pool <- colSums(comm_matrix)
  species_prob <- species_pool / sum(species_pool)
  species_names <- colnames(comm_matrix)
  
  rc_matrix <- matrix(NA, n_sites, n_sites)
  rownames(rc_matrix) <- colnames(rc_matrix) <- rownames(comm_matrix)
  
  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {
      sp_i <- comm_matrix[i, ]
      sp_j <- comm_matrix[j, ]
      
      alpha_i <- sum(sp_i)
      alpha_j <- sum(sp_j)
      
      shared_obs <- sum(sp_i * sp_j)
      
      shared_null <- numeric(reps)
      for (k in 1:reps) {
        sim_i <- sample(species_names, alpha_i, prob = species_prob, replace = FALSE)
        sim_j <- sample(species_names, alpha_j, prob = species_prob, replace = FALSE)
        shared_null[k] <- length(intersect(sim_i, sim_j))
      }
      
      # Calculate p-value (two-tailed)
      p_upper <- mean(shared_null >= shared_obs)
      p_lower <- mean(shared_null <= shared_obs)
      
      # Scale to -1 to 1
      rc <- if (shared_obs > mean(shared_null)) {
        1 - p_upper
      } else {
        -(1 - p_lower)
      }
      
      rc_matrix[i, j] <- rc
      rc_matrix[j, i] <- rc
    }
  }
  
  return(as.dist(rc_matrix))
}



format_raup <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute raup
  rc_dist <- custom_raup_crick(df_wide, reps = 999)
  
  # format raup
  raup_crick_matrix <- as.matrix(rc_dist)
  consecutive_dissimilarity <- sapply(1:(nrow(raup_crick_matrix) - 1), function(i) {
    raup_crick_matrix[i, i + 1]
  })
  names(consecutive_dissimilarity) <- paste0("T", 1:(length(consecutive_dissimilarity)), "_vs_T", 2:(length(consecutive_dissimilarity)+1))
  consecutive_dissimilarity <- (consecutive_dissimilarity-.5)*2

  result <- data.frame(
    unique_id = unique(df$unique_id),
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))], 
                      sort(unique(df$time))[-1], sep = " - "),
    raup_dissimilarity = consecutive_dissimilarity
  )
  
  return(result)
}


rc_dist <- custom_raup_crick(srs_sub_wider, reps = 999)
raup_crick_matrix <- as.matrix(rc_dist)
consecutive_dissimilarity <- sapply(1:(nrow(srs_sub_wider) - 1), function(i) {
  raup_crick_matrix[i, i + 1]
})
names(consecutive_dissimilarity) <- paste0("T", 1:(length(consecutive_dissimilarity)), "_vs_T", 2:(length(consecutive_dissimilarity)+1))
consecutive_dissimilarity <- (consecutive_dissimilarity-.5)*2

consecutive_dissimilarity <- as.data.frame(consecutive_dissimilarity)
consecutive_dissimilarity <- consecutive_dissimilarity %>%
  rownames_to_column("time")



prepare_comm_matrix <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
}

srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() 

raup_results <- srs_data_split %>%
  lapply(format_raup)

raup_results <- raup_results %>%
  bind_rows() %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) # converting to numeric

raup_results %>%
  #mutate(raup_dissimilarity = (raup_dissimilarity-.5)*2) %>%
  #filter(!raup_dissimilarity>-0.99999) %>%
  ggplot(aes(time2, raup_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_jitter(size = 2, alpha = 0.7) + 
  #geom_line(aes(time2, raup_dissimilarity, group = patch_type), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2")# +
  #stat_smooth(method = "lm", se = F, linewidth = 3)


# applying function to data
#raup_results <- srs_data_split %>% 
 # lapply(prepare_comm_matrix) %>%
 # lapply(custom_raup_crick)

#raup_results_matrix <- raup_results %>%
#  lapply(as.matrix) %>%
 # lapply(format_raup) 



################


# Simulated community data
comm <- matrix(c(
  1, 0, 1, 0, 1,
  0, 1, 1, 0, 0,
  1, 1, 0, 1, 0,
  1, 0, 0, 1, 1
), nrow = 4, byrow = TRUE)
rownames(comm) <- paste0("T", 1:4)
colnames(comm) <- paste0("Sp", 1:5)

# Run the custom Raup-Crick
rc_dist <- sapply(1:(nrow(comm) - 1), function(i) {
  vegan::vegdist(comm[i:(i+1), ], method = "raup")
})
rc_dist <- vegan::vegdist(comm, method = "raup")
# View dissimilarity matrix
as.matrix(rc_dist)













beta_raup <- raupcrick (decostand(srs_sub_wider, "pa"), 
                        nsimul = 9999, chase = F)
raup_crick_matrix <- as.matrix(beta_raup)
consecutive_dissimilarity <- sapply(1:(nrow(srs_sub_wider) - 1), function(i) {
  raup_crick_matrix[i, i + 1]
})
# Name the output for clarity
names(consecutive_dissimilarity) <- paste0("T", 1:(length(consecutive_dissimilarity)), "_vs_T", 2:(length(consecutive_dissimilarity)+1))

# Print results
print(consecutive_dissimilarity)
consecutive_dissimilarity <- (consecutive_dissimilarity-.5)*2
consecutive_dissimilarity


beta_raup <- sapply(1:(nrow(srs_sub_wider) - 1), function(i) {
  vegan::raupcrick(srs_sub_wider[i:(i+1), ], nsimul = 9999, chase = F)
})
beta_raup <- (beta_raup-.5)*2
beta_raup




# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# applying function to data
raup_results <- srs_data_split %>% 
  lapply(compute_raup) %>%
  bind_rows() # putting together into a dataframe

raup_results <- raup_results %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) # converting to numeric


raup_results %>%
  filter(unique_id == "53S-B-connected")


raup_results %>%
 # mutate(raup_dissimilarity = (raup_dissimilarity-.5)*2) %>%
  #filter(!raup_dissimilarity>-0.99999) %>%
  ggplot(aes(time2, raup_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_jitter(size = 2, alpha = 0.7) + 
  #geom_line(aes(time2, raup_dissimilarity, group = patch_type), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") #+
  #stat_smooth(method = "lm", se = F, linewidth = 3)


  



compute_raup <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute jaccard dissimilarity iteratively between consecutive years
  raup_values <- sapply(1:(nrow(df_wide) - 1), function(i) {
    vegan::vegdist(df_wide[i:(i+1), ], method = "raup")
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id),
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))], 
                      sort(unique(df$time))[-1], sep = " - "),
    raup_dissimilarity = raup_values
  )
  
  return(result)
}





