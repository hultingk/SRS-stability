# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, vegan) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)

# writing function to calculate jaccard's dissimilarity iteratively between consecutive years
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


#### ALL SPECIES ####
# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# applying function to data
jaccard_results <- srs_data_split %>% 
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe


#### COMMON SPECIES ####
# repeating excluding rare species (species with less than 10 occurrences across all surveys)

# figuring out which species are common, assigning a value to those to keep later 
common_spp <- srs_data %>% 
  count(sppcode) %>% # counting # of occurrences per species across all surveys
  filter(n > 10) %>% # removing species that had less than 10 occurrences
  mutate(rare = 0) %>% #  assigning a 0 to common species, rare will be NA later when joining
  dplyr::select(!c("n")) # removing counts of occurrences

commom_spp_data <- srs_data %>%
  left_join(common_spp, by = c("sppcode")) %>% # joining data on which species are common
  filter(rare == 0) %>% # only keeping common species
  dplyr::select(!c("rare")) # removing column

common_spp_jaccard <- commom_spp_data %>%
  dplyr::count(unique_id, time, sppcode) %>% # preparing data for function
  group_by(unique_id) %>%
  group_split() %>% # splitting dataframe into each patch as a list for function
  lapply(compute_jaccard) %>% # applying function to common species
  bind_rows() %>% # putting together into a dataframe
  rename(COMMON_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         COMMON_turnover = turnover_values, 
         COMMON_nestedness = nestedness_values)

#### LONGLEAF SPECIES (excluding rare) ####
# IGNORE FOR NOW 

#### DISPERSAL MODE (excluding rare) ####
## gravity dispersed
gravity_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Gravity") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(GRAVITY_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         GRAVITY_turnover = turnover_values, 
         GRAVITY_nestedness = nestedness_values)

## wind dispersed 
wind_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Wind") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(WIND_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         WIND_turnover = turnover_values, 
         WIND_nestedness = nestedness_values)

## animal dispersed 
animal_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Animal") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(ANIMAL_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         ANIMAL_turnover = turnover_values, 
         ANIMAL_nestedness = nestedness_values)




#### Temporal gamma diversity (spp pool size of each patch?)

# gamma diversity across time
gamma_diversity <- commom_spp_data %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(gamma_diversity = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, gamma_diversity)

# gamma of gravity dispersed
GRAVITY_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Gravity") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(GRAVITY_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, GRAVITY_gamma)

# gamma of wind dispersed
WIND_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Wind") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(WIND_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, WIND_gamma)

# gamma of animal dispersed
ANIMAL_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Animal") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(ANIMAL_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, ANIMAL_gamma)




#### joining all together ####
jaccard_results <- jaccard_results %>%
  left_join(common_spp_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(gravity_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(wind_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(animal_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) %>% # converting to numeric
  mutate(years_bw_surveys = time2 - time1) %>% # calculating # of years between surveys (most are annual)
  left_join(gamma_diversity, by = c("unique_id")) %>%
  left_join(GRAVITY_gamma, by = c("unique_id")) %>%
  left_join(WIND_gamma, by = c("unique_id")) %>%
  left_join(ANIMAL_gamma, by = c("unique_id"))

# adding soil moisture
soil_moisture <- srs_data %>%
  count(block, patch, soil_moisture) %>%
  dplyr::select(!n)

jaccard_results <- jaccard_results %>%
  left_join(soil_moisture, by = c("block", "patch"))

# writing summarized file
write_csv(jaccard_results, file = file.path("data", "L2_summarized", "patch_jaccard.csv"))
