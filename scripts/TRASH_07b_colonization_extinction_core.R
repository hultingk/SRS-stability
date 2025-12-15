librarian::shelf(tidyverse, vegan, purrr)
remotes::install_github("communityecologist/ecopart")

# loading data
srs_data_core <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_core.csv"))

srs_data_core <- srs_data_core %>% # removing experimentally planted species 
  filter(transplant != TRUE)# %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%

srs_data_core_wider <- srs_data_core %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
  arrange(unique_id, time) %>%  # Ensure years are sorted properly
  mutate(unique_id_time = paste(unique_id, time, sep = "-")) %>%
  dplyr::select(-patch_type, -patch, -year) 

block_list <- split(srs_data_core_wider, srs_data_core_wider$block)


block_time_list <- lapply(block_list, function(block_df) {
  time_list <- split(block_df, block_df$time)
  
  time_patch_pairs <- lapply(time_list, function(time) {
    patches <- unique(time$unique_id)
    combs <- combn(patches, 2, simplify = FALSE)
    
    # Create list of data.frames for each site pair
    pair_list <- lapply(combs, function(pair) {
      subset(time, unique_id %in% pair)
    })
    names(pair_list) <- sapply(combs, function(x) paste(x, collapse = "_vs_"))
    return(pair_list)
  })
}) ## redo this so it nested block:patch combo:time?


block_08 <- block_time_list[["08"]]
block_10 <- block_time_list[["10"]]
block_52 <- block_time_list[["52"]]
block_57 <- block_time_list[["57"]]
block_53N <- block_time_list[["53N"]]
block_53S <- block_time_list[["53S"]]
block_54N <- block_time_list[["54N"]]
block_54S <- block_time_list[["54S"]]
block_75E <- block_time_list[["75E"]]
block_75W <- block_time_list[["75W"]]


# 8
block_08_ecopart_four <- pairwise_ecopart(block_08, components = "four")
block_08_ecopart_four_df <- results_to_df(block_08_ecopart_four)
block_08_ecopart_two <- pairwise_ecopart(block_08, components = "two")
block_08_ecopart_two_df <- results_to_df(block_08_ecopart_two)

# 10
block_10_ecopart_four <- pairwise_ecopart(block_10, components = "four")
block_10_ecopart_four_df <- results_to_df(block_10_ecopart_four)
block_10_ecopart_two <- pairwise_ecopart(block_10, components = "two")
block_10_ecopart_two_df <- results_to_df(block_10_ecopart_two)

# 52
block_52_ecopart_four <- pairwise_ecopart(block_52, components = "four")
block_52_ecopart_four_df <- results_to_df(block_52_ecopart_four)
block_52_ecopart_two <- pairwise_ecopart(block_52, components = "two")
block_52_ecopart_two_df <- results_to_df(block_52_ecopart_two)

# 57
block_57_ecopart_four <- pairwise_ecopart(block_57, components = "four")
block_57_ecopart_four_df <- results_to_df(block_57_ecopart_four)
block_57_ecopart_two <- pairwise_ecopart(block_57, components = "two")
block_57_ecopart_two_df <- results_to_df(block_57_ecopart_two)

# 53S
block_53S_ecopart_four <- pairwise_ecopart(block_53S, components = "four")
block_53S_ecopart_four_df <- results_to_df(block_53S_ecopart_four)
block_53S_ecopart_two <- pairwise_ecopart(block_53S, components = "two")
block_53S_ecopart_two_df <- results_to_df(block_53S_ecopart_two)

# 53N
block_53N_ecopart_four <- pairwise_ecopart(block_53N, components = "four")
block_53N_ecopart_four_df <- results_to_df(block_53N_ecopart_four)
block_53N_ecopart_two <- pairwise_ecopart(block_53N, components = "two")
block_53N_ecopart_two_df <- results_to_df(block_53N_ecopart_two)

# 54S
block_54S_ecopart_four <- pairwise_ecopart(block_54S, components = "four")
block_54S_ecopart_four_df <- results_to_df(block_54S_ecopart_four)
block_54S_ecopart_two <- pairwise_ecopart(block_54S, components = "two")
block_54S_ecopart_two_df <- results_to_df(block_54S_ecopart_two)

# 54N
block_54N_ecopart_four <- pairwise_ecopart(block_54N, components = "four")
block_54N_ecopart_four_df <- results_to_df(block_54N_ecopart_four)
block_54N_ecopart_two <- pairwise_ecopart(block_54N, components = "two")
block_54N_ecopart_two_df <- results_to_df(block_54N_ecopart_two)

# 75E
block_75E_ecopart_four <- pairwise_ecopart(block_75E, components = "four")
block_75E_ecopart_four_df <- results_to_df(block_75E_ecopart_four)
block_75E_ecopart_two <- pairwise_ecopart(block_75E, components = "two")
block_75E_ecopart_two_df <- results_to_df(block_75E_ecopart_two)

# 75W
block_75W_ecopart_four <- pairwise_ecopart(block_75W, components = "four")
block_75W_ecopart_four_df <- results_to_df(block_75W_ecopart_four)
block_75W_ecopart_two <- pairwise_ecopart(block_75W, components = "two")
block_75W_ecopart_two_df <- results_to_df(block_75W_ecopart_two)


# combining all together
all_ecopart_four <- rbind(
  block_08_ecopart_four_df, block_10_ecopart_four_df, block_52_ecopart_four_df,
  block_57_ecopart_four_df, block_53S_ecopart_four_df, block_53N_ecopart_four_df, 
  block_54S_ecopart_four_df, block_54N_ecopart_four_df, block_75E_ecopart_four_df, 
  block_75W_ecopart_four_df
)

all_ecopart_two <- rbind(
  block_08_ecopart_two_df, block_10_ecopart_two_df, block_52_ecopart_two_df,
  block_57_ecopart_two_df, block_53S_ecopart_two_df, block_53N_ecopart_two_df, 
  block_54S_ecopart_two_df, block_54N_ecopart_two_df, block_75E_ecopart_two_df, 
  block_75W_ecopart_two_df
)

# wrangling
all_ecopart_four <- all_ecopart_four %>%
  separate(comparison, into = c("unique_id1", "unique_id2"), sep = "_vs_", remove = F) %>%
  separate(unique_id1, into = c("block1", "patch1", "patch_type1"), sep = "-", remove = F) %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-", remove = F) %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
  separate(year_pair, into = c("time1", "time2"), sep = "_vs_", remove = F) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Winged-Rectangular",
    patch_pair %in% c("connected-rectangle") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing") ~ "Connected-Winged",
    patch_pair %in% c("center-connected") ~ "Center-Connected",
    patch_pair %in% c("center-rectangle") ~ "Center-Rectangular",
    patch_pair %in% c("center-wing") ~ "Center-Winged",
    .default = patch_pair
  )) %>%
  mutate(center_peripheral = dplyr::case_when(
    patch_pair %in% c("Winged-Rectangular", "Connected-Rectangular", "Connected-Winged") ~ "Peripheral Comparison",
    patch_pair %in% c("Center-Connected", "Center-Rectangular", "Center-Winged") ~ "Center Comparison",
    .default = patch_pair
  ))

all_ecopart_two <- all_ecopart_two %>%
  separate(comparison, into = c("unique_id1", "unique_id2"), sep = "_vs_", remove = F) %>%
  separate(unique_id1, into = c("block1", "patch1", "patch_type1"), sep = "-", remove = F) %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-", remove = F) %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
  separate(year_pair, into = c("time1", "time2"), sep = "_vs_", remove = F) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Winged-Rectangular",
    patch_pair %in% c("connected-rectangle") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing") ~ "Connected-Winged",
    patch_pair %in% c("center-connected") ~ "Center-Connected",
    patch_pair %in% c("center-rectangle") ~ "Center-Rectangular",
    patch_pair %in% c("center-wing") ~ "Center-Winged",
    .default = patch_pair
  )) %>%
  mutate(center_peripheral = dplyr::case_when(
    patch_pair %in% c("Winged-Rectangular", "Connected-Rectangular", "Connected-Winged") ~ "Peripheral Comparison",
    patch_pair %in% c("Center-Connected", "Center-Rectangular", "Center-Winged") ~ "Center Comparison",
    .default = patch_pair
  ))




ecopart_plot_four <- all_ecopart_four %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  scale_fill_brewer(palette = "BrBG", name = NULL, labels = c(expression(paste("Colonization differentiation (", Delta, beta["C+"], ")")),
                                                              expression(paste("Colonization homogenization (", Delta, beta["C-"], ")")),
                                                              expression(paste("Extinction differentiation (", Delta, beta["E+"], ")")),
                                                              expression(paste("Extinction homogenization (", Delta, beta["E-"], ")")))) +
  scale_color_brewer(palette = "BrBG", name = NULL, labels = c(expression(paste("Colonization differentiation (", Delta, beta["C+"], ")")),
                                                               expression(paste("Colonization homogenization (", Delta, beta["C-"], ")")),
                                                               expression(paste("Extinction differentiation (", Delta, beta["E+"], ")")),
                                                               expression(paste("Extinction homogenization (", Delta, beta["E-"], ")")))) +
  ylab(expression(paste("Change in ", beta , " diversity"))) +
  xlab("Time since site creation (years)") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 4)) +
  theme(panel.spacing = unit(1.5, "lines"))
ecopart_plot_four


ecopart_plot_two <- all_ecopart_two %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), name = NULL, labels = c(expression(paste("Colonization component (", Delta, beta["C"], ")")),
                                                                              expression(paste("Extinction component (", Delta, beta["E"], ")")))) +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"), name = NULL, labels = c(expression(paste("Colonization component (", Delta, beta["C"], ")")),
                                                                               expression(paste("Extinction component (", Delta, beta["E"], ")")))) +
  ylab(expression(paste("Change in ", beta , " diversity"))) +
  xlab(NULL) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 2)) +
  theme(panel.spacing = unit(1.5, "lines"))
ecopart_plot_two
