# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, vegan)

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# writing function to calculate jaccard's dissimilarity iteratively between consectutive years
compute_jaccard <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute jaccard dissimilarity, nestedness, and turnover iteratively between consecutive years
  jaccard_values <- sapply(1:(nrow(df_wide) - 1), function(i) {
    nestedbetajac(df_wide[i:(i+1), ])
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id), # unique id
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))],  # year pairs
                      sort(unique(df$time))[-1], sep = " - "),
    jaccard_dissimilarity = jaccard_values[row.names(jaccard_values) %in% c("jaccard"),], # jaccard values
    turnover_values = jaccard_values[row.names(jaccard_values) %in% c("turnover"),], # turnover values
    nestedness_values = jaccard_values[row.names(jaccard_values) %in% c("nestedness"),] # nestedness values
  )
  
  return(result)
}

# applying function to data
jaccard_results <- srs_data_split %>% 
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

# converting year pairs into a usable format for modeling/plotting
jaccard_results <- jaccard_results %>%
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  mutate(years_bw_surveys = time2 - time1)
  
jaccard_results <- jaccard_results %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) # separating unique ID into columns





# repeating excluding rare species (species with less than 10 occurances across all surveys)
common_spp <- srs_data %>%
  count(sppcode) %>%
  arrange(n) %>%
  filter(n > 10) %>%
  mutate(rare = 0) %>%
  select(!c("n"))

common_spp_jaccard <- srs_data %>%
  left_join(common_spp, by = c("sppcode")) %>%
  filter(rare == 0) %>%
  dplyr::select(!c("rare")) %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

common_spp_jaccard <- common_spp_jaccard %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  mutate(years_bw_surveys = time2 - time1)

