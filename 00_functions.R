# loading libraries
librarian::shelf(tidyverse, vegan) # Install missing packages and load needed libraries


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


