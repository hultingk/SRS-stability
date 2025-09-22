librarian::shelf(tidyverse, vegan, purrr, ecopart)
#remotes::install_github("communityecologist/ecopart")

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")


srs_data_wider <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
  arrange(unique_id, time) %>%  # Ensure years are sorted properly
   mutate(unique_id_time = paste(unique_id, time, sep = "-")) %>%
  dplyr::select(-patch_type, -patch, -year) 


block_list <- split(srs_data_wider, srs_data_wider$block)
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

results_to_df <- function(results) {
  map_dfr(names(results), function(year_pair) {
    comps <- results[[year_pair]]
    
    map_dfr(names(comps), function(comp) {
      vals <- comps[[comp]]
      
      tibble(
        year_pair = year_pair,
        comparison = comp,
        metric = names(vals),
        value = as.numeric(vals)
      )
    })
  })
}

# 8
block_08_ecopart_four <- pairwise_ecopart(block_08, components = "four")
block_08_ecopart_four_df <- results_to_df(block_08_ecopart_four)
block_08_ecopart_two <- pairwise_ecopart(block_08, components = "two")
block_08_ecopart_two_df <- results_to_df(block_08_ecopart_two)
block_08_ecopart_six <- pairwise_ecopart(block_08, components = "six")
block_08_ecopart_six_df <- results_to_df(block_08_ecopart_six)

# 10
block_10_ecopart_four <- pairwise_ecopart(block_10, components = "four")
block_10_ecopart_four_df <- results_to_df(block_10_ecopart_four)
block_10_ecopart_two <- pairwise_ecopart(block_10, components = "two")
block_10_ecopart_two_df <- results_to_df(block_10_ecopart_two)
block_10_ecopart_six <- pairwise_ecopart(block_10, components = "six")
block_10_ecopart_six_df <- results_to_df(block_10_ecopart_six)

# 52
block_52_ecopart_four <- pairwise_ecopart(block_52, components = "four")
block_52_ecopart_four_df <- results_to_df(block_52_ecopart_four)
block_52_ecopart_two <- pairwise_ecopart(block_52, components = "two")
block_52_ecopart_two_df <- results_to_df(block_52_ecopart_two)
block_52_ecopart_six <- pairwise_ecopart(block_52, components = "six")
block_52_ecopart_six_df <- results_to_df(block_52_ecopart_six)

# 57
block_57_ecopart_four <- pairwise_ecopart(block_57, components = "four")
block_57_ecopart_four_df <- results_to_df(block_57_ecopart_four)
block_57_ecopart_two <- pairwise_ecopart(block_57, components = "two")
block_57_ecopart_two_df <- results_to_df(block_57_ecopart_two)
block_57_ecopart_six <- pairwise_ecopart(block_57, components = "six")
block_57_ecopart_six_df <- results_to_df(block_57_ecopart_six)

# 53S
block_53S_ecopart_four <- pairwise_ecopart(block_53S, components = "four")
block_53S_ecopart_four_df <- results_to_df(block_53S_ecopart_four)
block_53S_ecopart_two <- pairwise_ecopart(block_53S, components = "two")
block_53S_ecopart_two_df <- results_to_df(block_53S_ecopart_two)
block_53S_ecopart_six <- pairwise_ecopart(block_53S, components = "six")
block_53S_ecopart_six_df <- results_to_df(block_53S_ecopart_six)

# 53N
block_53N_ecopart_four <- pairwise_ecopart(block_53N, components = "four")
block_53N_ecopart_four_df <- results_to_df(block_53N_ecopart_four)
block_53N_ecopart_two <- pairwise_ecopart(block_53N, components = "two")
block_53N_ecopart_two_df <- results_to_df(block_53N_ecopart_two)
block_53N_ecopart_six <- pairwise_ecopart(block_53N, components = "six")
block_53N_ecopart_six_df <- results_to_df(block_53N_ecopart_six)

# 54S
block_54S_ecopart_four <- pairwise_ecopart(block_54S, components = "four")
block_54S_ecopart_four_df <- results_to_df(block_54S_ecopart_four)
block_54S_ecopart_two <- pairwise_ecopart(block_54S, components = "two")
block_54S_ecopart_two_df <- results_to_df(block_54S_ecopart_two)
block_54S_ecopart_six <- pairwise_ecopart(block_54S, components = "six")
block_54S_ecopart_six_df <- results_to_df(block_54S_ecopart_six)

# 54N
block_54N_ecopart_four <- pairwise_ecopart(block_54N, components = "four")
block_54N_ecopart_four_df <- results_to_df(block_54N_ecopart_four)
block_54N_ecopart_two <- pairwise_ecopart(block_54N, components = "two")
block_54N_ecopart_two_df <- results_to_df(block_54N_ecopart_two)
block_54N_ecopart_six <- pairwise_ecopart(block_54N, components = "six")
block_54N_ecopart_six_df <- results_to_df(block_54N_ecopart_six)

# 75E
block_75E_ecopart_four <- pairwise_ecopart(block_75E, components = "four")
block_75E_ecopart_four_df <- results_to_df(block_75E_ecopart_four)
block_75E_ecopart_two <- pairwise_ecopart(block_75E, components = "two")
block_75E_ecopart_two_df <- results_to_df(block_75E_ecopart_two)
block_75E_ecopart_six <- pairwise_ecopart(block_75E, components = "six")
block_75E_ecopart_six_df <- results_to_df(block_75E_ecopart_six)

# 75W
block_75W_ecopart_four <- pairwise_ecopart(block_75W, components = "four")
block_75W_ecopart_four_df <- results_to_df(block_75W_ecopart_four)
block_75W_ecopart_two <- pairwise_ecopart(block_75W, components = "two")
block_75W_ecopart_two_df <- results_to_df(block_75W_ecopart_two)
block_75W_ecopart_six <- pairwise_ecopart(block_75W, components = "six")
block_75W_ecopart_six_df <- results_to_df(block_75W_ecopart_six)

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

all_ecopart_six <- rbind(
  block_08_ecopart_six_df, block_10_ecopart_six_df, block_52_ecopart_six_df,
  block_57_ecopart_six_df, block_53S_ecopart_six_df, block_53N_ecopart_six_df, 
  block_54S_ecopart_six_df, block_54N_ecopart_six_df, block_75E_ecopart_six_df, 
  block_75W_ecopart_six_df
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
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Winged-Rectangular",
    patch_pair %in% c("connected-rectangle") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing") ~ "Connected-Winged",
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
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Winged-Rectangular",
    patch_pair %in% c("connected-rectangle") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing") ~ "Connected-Winged",
    .default = patch_pair
  ))

all_ecopart_six <- all_ecopart_six %>%
  separate(comparison, into = c("unique_id1", "unique_id2"), sep = "_vs_", remove = F) %>%
  separate(unique_id1, into = c("block1", "patch1", "patch_type1"), sep = "-", remove = F) %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-", remove = F) %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
  separate(year_pair, into = c("time1", "time2"), sep = "_vs_", remove = F) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("wing-rectangle", "rectangle-wing") ~ "Winged-Rectangular",
    patch_pair %in% c("connected-rectangle") ~ "Connected-Rectangular",
    patch_pair %in% c("connected-wing") ~ "Connected-Winged",
    .default = patch_pair
  ))


ecopart_plot_four <- all_ecopart_four %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
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

# pdf(file = file.path("plots", "ecopart_plot.pdf"), width = 13, height = 8)
# ecopart_plot
# dev.off()

ecopart_plot_two <- all_ecopart_two %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
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




ecopart_plot_six <- all_ecopart_six %>%
  filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  #scale_fill_manual(values = c("#E1BE6A", "#40B0A6"), name = NULL, labels = c(expression(paste("Colonization component (", Delta, beta["C"], ")")),
 #                                                                             expression(paste("Extinction component (", Delta, beta["E"], ")")))) +
 # scale_color_manual(values = c("#E1BE6A", "#40B0A6"), name = NULL, labels = c(expression(paste("Colonization component (", Delta, beta["C"], ")")),
 #                                                                              expression(paste("Extinction component (", Delta, beta["E"], ")")))) +
  ylab(expression(paste("Change in ", beta , " diversity"))) +
  xlab(NULL) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 2)) +
  theme(panel.spacing = unit(1.5, "lines"))
ecopart_plot_six

pdf(file = file.path("plots", "ecopart_plot.pdf"), width = 14, height = 9.5)
cowplot::plot_grid(ecopart_plot_two, ecopart_plot_four, nrow = 2, axis = "tblr", align = "hv", 
                   labels = c("(A)", "(B)"),  label_x = 0.03, label_y = 1, label_size = 20)
dev.off()

# all_ecopart %>%
#   filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
#   ggplot(aes(metric, value, fill = patch_pair)) +
#   geom_boxplot()
# 
# all_ecopart %>%
#   group_by(patch_pair, metric) %>%
#   summarize(avg_value = mean(value))

# taking the absolute value for comparison
all_ecopart_two <- all_ecopart_two %>%
  mutate(abs_value = abs(value)) %>%
  filter(!patch_pair %in% c("rectangle-rectangle", "wing-wing"))

# dividing into patch pair types for modeling
c_r_ecopart <- all_ecopart_two %>%
  filter(patch_pair == "Connected-Rectangular")
c_w_ecopart <- all_ecopart_two %>%
  filter(patch_pair == "Connected-Winged")
w_r_ecopart <- all_ecopart_two %>%
  filter(patch_pair == "Winged-Rectangular")

# models
m.c_r_ecopart <- glmmTMB(abs_value ~ metric*time1 + (1|block1),
              data = c_r_ecopart)
summary(m.c_r_ecopart)

m.c_w_ecopart <- glmmTMB(abs_value ~ metric*time1 + (1|block1),
                         data = c_w_ecopart)
summary(m.c_w_ecopart)

m.w_r_ecopart <- glmmTMB(abs_value ~ metric*time1  + (1|block1),
                         data = w_r_ecopart)
summary(m.w_r_ecopart)


