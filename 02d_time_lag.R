librarian::shelf(tidyverse, vegan, glmmTMB) # Install missing packages and load needed libraries


#### computing time lag jaccard values ####
compute_time_lag <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # Calculate pairwise dissimilarities
  beta_dist <- vegdist(df_wide, method = "jaccard") 
  
  # Create time lag matrix
  n <- nrow(df_wide)
  time_lags <- abs(outer(1:n, 1:n, "-"))
  time_lag_values <- time_lags[lower.tri(time_lags)]
  beta_values <- beta_dist[lower.tri(as.matrix(beta_dist))]
  
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id), # unique id
    time_lags = time_lag_values,
    beta_values = beta_values
    )
  
  return(result)
}


### applying function
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)



#### ALL SPECIES ####
# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# applying function to data
time_lag_results <- srs_data_split %>% 
  lapply(compute_time_lag) %>%
  bind_rows() # putting together into a dataframe

time_lag_results <- time_lag_results %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  filter(patch_type != "center")

m1 <- glmmTMB(beta_values ~ patch_type*time_lags + (1|block/patch),
              data = time_lag_results)
summary(m1)

time_lag_results %>%
  ggplot() +
  geom_point(aes(time_lags, beta_values, color = patch_type)) +
  geom_smooth(aes(time_lags, beta_values, color = patch_type), method = "lm") +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2")
