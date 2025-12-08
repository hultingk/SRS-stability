# loading libraries
librarian::shelf(tidyverse, vegan) # Install missing packages and load needed libraries

### compute jaccard dissimilarity between patch pairs for each year
compute_convergence_jaccard <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # Ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
    dplyr::select(-block, -patch_type, -patch, -year, -time) %>%
    column_to_rownames("unique_id")
  
  jac <- vegdist(df_wide, method = "jaccard")
  jac_df <- as.data.frame(as.table(as.matrix(jac)))
  names(jac_df) <- c("unique_id1", "unique_id2", "jaccard")
  
  jac_df <- jac_df %>%
    separate(unique_id1, into = c("block1", "patch_rep1", "patch_type1", "time1"), remove = T) %>%
    separate(unique_id2, into = c("block2", "patch_rep2", "patch_type2", "time2"), remove = T) %>%
    mutate(unique_id1 = paste(block1, patch_rep1, patch_type1, sep = "-")) %>%
    mutate(unique_id2 = paste(block2, patch_rep2, patch_type2, sep = "-")) %>%
    filter(!(unique_id1 == unique_id2)) %>%
    filter(time1 == time2) %>%
    mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
    dplyr::select(-block2, -time2, -patch_rep1, -patch_rep2, -patch_type1, -patch_type2) %>%
    rename(block = block1, time = time1) %>%
    mutate(time = as.numeric(time)) %>%
    mutate(patch_pair = dplyr::case_when(
      patch_pair %in% c("Connected-Rectangular", "Rectangular-Connected") ~ "Connected-Rectangular",
      patch_pair %in% c("Connected-Winged", "Winged-Connected") ~ "Connected-Winged",
      patch_pair %in% c("Rectangular-Winged", "Winged-Rectangular") ~ "Rectangular-Winged",
      patch_pair %in% c("Center-Connected", "Connected-Center") ~ "Center-Connected",
      patch_pair %in% c("Center-Rectangular", "Rectangular-Center") ~ "Center-Rectangular",
      patch_pair %in% c("Center-Winged", "Winged-Center") ~ "Center-Winged",
      .default = patch_pair
    ))
  
  return(jac_df)
}


#### partition changes in spatial beta diversity into changes due to colonization or extinction ####
pairwise_ecopart <- function(nested_list, components) {
  
  # make sure times are sorted numerically
  time_keys <- sort(as.numeric(names(nested_list)))
  results <- map2(time_keys[-length(time_keys)], time_keys[-1], function(t1, t2) {
    comp_names <- intersect(names(nested_list[[as.character(t1)]]),
                            names(nested_list[[as.character(t2)]]))
    
    comp_res <- map(comp_names, function(cmp) {
      mat1 <- nested_list[[as.character(t1)]][[cmp]]
      mat2 <- nested_list[[as.character(t2)]][[cmp]]
      # drop metadata columns, keep only species
      mat1 <- mat1 %>% dplyr::select(-block, -unique_id, -time, -unique_id_time)
      mat2 <- mat2 %>% dplyr::select(-block, -unique_id, -time, -unique_id_time)
      # make sure same species order
      common_species <- intersect(names(mat1), names(mat2))
      mat1 <- mat1[, common_species, drop = FALSE]
      mat2 <- mat2[, common_species, drop = FALSE]
      
      ecopart::ecopart.pair(mat1, mat2, components = components, index = "jaccard")
    })
    
    names(comp_res) <- comp_names
    comp_res
  })
  
  names(results) <- paste(time_keys[-length(time_keys)], time_keys[-1], sep = "_vs_")
  results
}


#### convert pairwise_ecopart results into dataframe ####
results_to_df <- function(results) { # taking a nested list and flattening it - two metrics nested within patch comparison nested within year comparison
  map_dfr(names(results), function(year_pair) { # looping over each year pair
    comps <- results[[year_pair]]
    
    map_dfr(names(comps), function(comp) { # looping over each patch pair
      vals <- comps[[comp]] # getting metrics
      
      tibble(
        year_pair = year_pair,
        comparison = comp,
        metric = names(vals),
        value = as.numeric(vals)
      )
    })
  })
}


#### calculate raw composition change between consecutive years ####
compute_composition_change <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
    dplyr::select(-block, -patch_type, -patch, -year) %>% # remove unneeded columns
    column_to_rownames("unique_id")
  
  # calculate losses and gains between consecutive years
  turnover <- df_wide %>%
    arrange(time) %>% # ensure years are sorted properly
    mutate(
      gains = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 0)), na.rm = TRUE),     # calculate sp gains, where 0 -> 1
      losses = rowSums(across(-time, ~ (.x == 0 & lag(.x) == 1)), na.rm = TRUE),    # calculate sp losses, 1 -> 0
      stayed_present = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 1)), na.rm = TRUE),  # calculate sp staying consistent, 1 -> 1
      changed_total = gains + losses # total turnover
    ) %>%
    dplyr::select(time, gains, losses, stayed_present, changed_total)
  
  
  return(turnover)
}
