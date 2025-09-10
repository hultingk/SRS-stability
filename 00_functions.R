# loading libraries
librarian::shelf(tidyverse, vegan) # Install missing packages and load needed libraries

### compute jaccard dissimilarity between patch pairs for each year
compute_convergence_jaccard <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # Ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
    select(-block, -patch_type, -patch, -year, -time) %>%
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
      patch_pair %in% c("connected-rectangle", "rectangle-connected") ~ "Connected-Rectangular",
      patch_pair %in% c("connected-wing", "wing-connected") ~ "Connected-Winged",
      patch_pair %in% c("rectangle-wing", "wing-rectangle") ~ "Rectangular-Winged",
      .default = patch_pair
    ))
  
  return(jac_df)
}




#### computing jaccard values ####
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
  
  # calculate alpha diversity (species richness) for each year
  alpha_diversity <- apply(df_wide, 1, function(row) sum(row > 0))
  # calculate mean alpha diversity for each year pair
  mean_alpha <- sapply(1:(length(alpha_diversity) - 1), function(i) {
    mean(c(alpha_diversity[i], alpha_diversity[i + 1]))
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id), # unique id
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))],  # year pairs
                      sort(unique(df$time))[-1], sep = " - "),
    jaccard_dissimilarity = jaccard_values[row.names(jaccard_values) %in% c("jaccard"),], # jaccard values
    turnover_values = jaccard_values[row.names(jaccard_values) %in% c("turnover"),], # turnover values
    nestedness_values = jaccard_values[row.names(jaccard_values) %in% c("nestedness"),], # nestedness values
    mean_alpha_diversity = mean_alpha
  )
  
  return(result)
}




#### computing sorensen values ####
# writing function to calculate sorensen's dissimilarity iteratively between consecutive years
compute_sorensen <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute sorensen dissimilarity, nestedness, and turnover iteratively between consecutive years
  sorensen_values <- sapply(1:(nrow(df_wide) - 1), function(i) {
    nestedbetasor(df_wide[i:(i+1), ])
  })
  
  # calculate alpha diversity (species richness) for each year
  alpha_diversity <- apply(df_wide, 1, function(row) sum(row > 0))
  # calculate mean alpha diversity for each year pair
  mean_alpha <- sapply(1:(length(alpha_diversity) - 1), function(i) {
    mean(c(alpha_diversity[i], alpha_diversity[i + 1]))
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id), # unique id
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))],  # year pairs
                      sort(unique(df$time))[-1], sep = " - "),
    sorensen_dissimilarity = sorensen_values[row.names(sorensen_values) %in% c("sorensen"),], # jaccard values
    turnover_values = sorensen_values[row.names(sorensen_values) %in% c("turnover"),], # turnover values
    nestedness_values = sorensen_values[row.names(sorensen_values) %in% c("nestedness"),], # nestedness values
    mean_alpha_diversity = mean_alpha
  )
  
  return(result)
}



#### raup crick function ####
compute_raup <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute raup
  rc_dist <- RC.pc(df_wide, weighted = F, taxo.metric = "jaccard")
  
  # format raup
  raup_crick_matrix <- as.matrix(rc_dist[["index"]])
  consecutive_dissimilarity <- sapply(1:(nrow(raup_crick_matrix) - 1), function(i) {
    raup_crick_matrix[i, i + 1]
  })
  
  
  result <- data.frame(
    unique_id = unique(df$unique_id),
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))], 
                      sort(unique(df$time))[-1], sep = " - "),
    raup_dissimilarity = consecutive_dissimilarity
  )
  
  return(result)
}



