########
## SCRIPT NAME: 07_colonization_extinction.R
## AUTHOR: Katherine Hulting
## PURPOSE: Partition change in spatial beta diversity over time into changes due to colonization and extinction
## PRODUCTS: FigureS5.pdf
#########

# load libraries
librarian::shelf(tidyverse, vegan, purrr, ecopart)
#remotes::install_github("communityecologist/ecopart")

source(here::here(file.path("scripts", "00_functions.R"))) # loading functions

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "Center")

# pivot to wider format
srs_data_wider <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
  arrange(unique_id, time) %>%  # Ensure years are sorted properly
  mutate(unique_id_time = paste(unique_id, time, sep = "-")) %>%
  dplyr::select(-patch_type, -patch, -year) 

# split up by blocks 
block_list <- split(srs_data_wider, srs_data_wider$block)

# need to get into the correct format for ecopart
# create nested list of community matrices - each matrix is the community composition of the two patches being compared in one year
block_time_list <- lapply(block_list, function(block_df) {
  time_list <- split(block_df, block_df$time) # splitting each block in the list into a further list of one matrix per year
  
  time_patch_pairs <- lapply(time_list, function(time) { # apply to all years within a block
    patches <- unique(time$unique_id) # get unique patches
    patch_combinations <- combn(patches, 2, simplify = FALSE) # get all combinations of patches
    
    # create list of dataframes for each site pair
    pair_list <- lapply(patch_combinations, function(pair) {
      subset(time, unique_id %in% pair) # subsetting the patches in each pair for each year
    })
    names(pair_list) <- sapply(patch_combinations, function(x) paste(x, collapse = "_vs_")) # naming lists by what patches they are comparing
    return(pair_list)
  })
}) 

# dividing up by block - probably a more concise way to do this
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


# applying pairwise_ecopart function to every block, calculating two, four, or six components - definitely a more concise way to do this
# 8
block_08_ecopart_two <- pairwise_ecopart(block_08, components = "two")
block_08_ecopart_two_df <- results_to_df(block_08_ecopart_two)
block_08_ecopart_four <- pairwise_ecopart(block_08, components = "four")
block_08_ecopart_four_df <- results_to_df(block_08_ecopart_four)
block_08_ecopart_six <- pairwise_ecopart(block_08, components = "six")
block_08_ecopart_six_df <- results_to_df(block_08_ecopart_six)

# 10
block_10_ecopart_two <- pairwise_ecopart(block_10, components = "two")
block_10_ecopart_two_df <- results_to_df(block_10_ecopart_two)
block_10_ecopart_four <- pairwise_ecopart(block_10, components = "four")
block_10_ecopart_four_df <- results_to_df(block_10_ecopart_four)
block_10_ecopart_six <- pairwise_ecopart(block_10, components = "six")
block_10_ecopart_six_df <- results_to_df(block_10_ecopart_six)

# 52
block_52_ecopart_two <- pairwise_ecopart(block_52, components = "two")
block_52_ecopart_two_df <- results_to_df(block_52_ecopart_two)
block_52_ecopart_four <- pairwise_ecopart(block_52, components = "four")
block_52_ecopart_four_df <- results_to_df(block_52_ecopart_four)
block_52_ecopart_six <- pairwise_ecopart(block_52, components = "six")
block_52_ecopart_six_df <- results_to_df(block_52_ecopart_six)

# 57
block_57_ecopart_two <- pairwise_ecopart(block_57, components = "two")
block_57_ecopart_two_df <- results_to_df(block_57_ecopart_two)
block_57_ecopart_four <- pairwise_ecopart(block_57, components = "four")
block_57_ecopart_four_df <- results_to_df(block_57_ecopart_four)
block_57_ecopart_six <- pairwise_ecopart(block_57, components = "six")
block_57_ecopart_six_df <- results_to_df(block_57_ecopart_six)

# 53S
block_53S_ecopart_two <- pairwise_ecopart(block_53S, components = "two")
block_53S_ecopart_two_df <- results_to_df(block_53S_ecopart_two)
block_53S_ecopart_four <- pairwise_ecopart(block_53S, components = "four")
block_53S_ecopart_four_df <- results_to_df(block_53S_ecopart_four)
block_53S_ecopart_six <- pairwise_ecopart(block_53S, components = "six")
block_53S_ecopart_six_df <- results_to_df(block_53S_ecopart_six)

# 53N
block_53N_ecopart_two <- pairwise_ecopart(block_53N, components = "two")
block_53N_ecopart_two_df <- results_to_df(block_53N_ecopart_two)
block_53N_ecopart_four <- pairwise_ecopart(block_53N, components = "four")
block_53N_ecopart_four_df <- results_to_df(block_53N_ecopart_four)
block_53N_ecopart_six <- pairwise_ecopart(block_53N, components = "six")
block_53N_ecopart_six_df <- results_to_df(block_53N_ecopart_six)

# 54S
block_54S_ecopart_two <- pairwise_ecopart(block_54S, components = "two")
block_54S_ecopart_two_df <- results_to_df(block_54S_ecopart_two)
block_54S_ecopart_four <- pairwise_ecopart(block_54S, components = "four")
block_54S_ecopart_four_df <- results_to_df(block_54S_ecopart_four)
block_54S_ecopart_six <- pairwise_ecopart(block_54S, components = "six")
block_54S_ecopart_six_df <- results_to_df(block_54S_ecopart_six)

# 54N
block_54N_ecopart_two <- pairwise_ecopart(block_54N, components = "two")
block_54N_ecopart_two_df <- results_to_df(block_54N_ecopart_two)
block_54N_ecopart_four <- pairwise_ecopart(block_54N, components = "four")
block_54N_ecopart_four_df <- results_to_df(block_54N_ecopart_four)
block_54N_ecopart_six <- pairwise_ecopart(block_54N, components = "six")
block_54N_ecopart_six_df <- results_to_df(block_54N_ecopart_six)

# 75E
block_75E_ecopart_two <- pairwise_ecopart(block_75E, components = "two")
block_75E_ecopart_two_df <- results_to_df(block_75E_ecopart_two)
block_75E_ecopart_four <- pairwise_ecopart(block_75E, components = "four")
block_75E_ecopart_four_df <- results_to_df(block_75E_ecopart_four)
block_75E_ecopart_six <- pairwise_ecopart(block_75E, components = "six")
block_75E_ecopart_six_df <- results_to_df(block_75E_ecopart_six)

# 75W
block_75W_ecopart_two <- pairwise_ecopart(block_75W, components = "two")
block_75W_ecopart_two_df <- results_to_df(block_75W_ecopart_two)
block_75W_ecopart_four <- pairwise_ecopart(block_75W, components = "four")
block_75W_ecopart_four_df <- results_to_df(block_75W_ecopart_four)
block_75W_ecopart_six <- pairwise_ecopart(block_75W, components = "six")
block_75W_ecopart_six_df <- results_to_df(block_75W_ecopart_six)

# combining all together
# dataframe of results of two components
all_ecopart_two <- rbind(
  block_08_ecopart_two_df, block_10_ecopart_two_df, block_52_ecopart_two_df,
  block_57_ecopart_two_df, block_53S_ecopart_two_df, block_53N_ecopart_two_df, 
  block_54S_ecopart_two_df, block_54N_ecopart_two_df, block_75E_ecopart_two_df, 
  block_75W_ecopart_two_df
)
# dataframe of results of four components
all_ecopart_four <- rbind(
  block_08_ecopart_four_df, block_10_ecopart_four_df, block_52_ecopart_four_df,
  block_57_ecopart_four_df, block_53S_ecopart_four_df, block_53N_ecopart_four_df, 
  block_54S_ecopart_four_df, block_54N_ecopart_four_df, block_75E_ecopart_four_df, 
  block_75W_ecopart_four_df
)
# dataframe of results of six components
all_ecopart_six <- rbind(
  block_08_ecopart_six_df, block_10_ecopart_six_df, block_52_ecopart_six_df,
  block_57_ecopart_six_df, block_53S_ecopart_six_df, block_53N_ecopart_six_df, 
  block_54S_ecopart_six_df, block_54N_ecopart_six_df, block_75E_ecopart_six_df, 
  block_75W_ecopart_six_df
)

# wrangling each into correct format
all_ecopart_two <- all_ecopart_two %>%
  separate(comparison, into = c("unique_id1", "unique_id2"), sep = "_vs_", remove = F) %>% # separating comparison into two columns 
  separate(unique_id1, into = c("block1", "patch1", "patch_type1"), sep = "-", remove = F) %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-", remove = F) %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
  separate(year_pair, into = c("time1", "time2"), sep = "_vs_", remove = F) %>% # separating the two years
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("Winged-Rectangular", "Rectangular-Winged") ~ "Winged-Rectangular",
    .default = patch_pair
  ))

all_ecopart_four <- all_ecopart_four %>%
  separate(comparison, into = c("unique_id1", "unique_id2"), sep = "_vs_", remove = F) %>%
  separate(unique_id1, into = c("block1", "patch1", "patch_type1"), sep = "-", remove = F) %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-", remove = F) %>%
  mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>%
  separate(year_pair, into = c("time1", "time2"), sep = "_vs_", remove = F) %>%
  mutate(time1 = as.numeric(time1)) %>%
  mutate(time2 = as.numeric(time2)) %>%
  mutate(patch_pair = dplyr::case_when(
    patch_pair %in% c("Winged-Rectangular", "Rectangular-Winged") ~ "Winged-Rectangular",
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
    patch_pair %in% c("Winged-Rectangular", "Rectangular-Winged") ~ "Winged-Rectangular",
    .default = patch_pair
  ))


## plotting four components
ecopart_plot_four <- all_ecopart_four %>%
  filter(!patch_pair %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  scale_fill_brewer(palette = "BrBG", name = NULL, labels = c(expression(paste("Colonization differentiation (", Delta, beta["C+"], ")")),
                                                              expression(paste("Colonization homogenization (", Delta, beta["C-"], ")")),
                                                              expression(paste("Extinction differentiation (", Delta, beta["E+"], ")")),
                                                              expression(paste("Extinction homogenization (", Delta, beta["E-"], ")")))) +
  scale_color_brewer(palette = "BrBG", name = NULL, labels = c(expression(paste("Colonization differentiation (", Delta, beta["C+"], ")")),
                                                               expression(paste("Colonization homogenization (", Delta, beta["C-"], ")")),
                                                               expression(paste("Extinction differentiation (", Delta, beta["E+"], ")")),
                                                               expression(paste("Extinction homogenization (", Delta, beta["E-"], ")")))) +
  ylab(expression(paste("Change in ", beta , " diversity"))) +
  xlab("Years since site creation") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 4)) +
  theme(panel.spacing = unit(1.5, "lines"))
ecopart_plot_four



## plotting two components
ecopart_plot_two <- all_ecopart_two %>%
  filter(!patch_pair %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
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



# plotting six components
ecopart_plot_six <- all_ecopart_six %>%
  filter(!patch_pair %in% c("Winged-Winged", "Rectangular-Rectangular")) %>%
  ggplot(aes(time2, value, color = metric, fill = metric)) +
  geom_point(alpha = 0.09, size = 4) +
  stat_smooth(size = 3, method = "lm", formula = y ~ x + I(x^2)) +
  facet_wrap(~patch_pair) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  ylab(expression(paste("Change in ", beta , " diversity"))) +
  xlab(NULL) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 2)) +
  theme(panel.spacing = unit(1.5, "lines"))
ecopart_plot_six


# combining the plot of two and four components into one plot
figureS4 <- cowplot::plot_grid(ecopart_plot_two, ecopart_plot_four, nrow = 2, axis = "tblr", align = "hv", 
                   labels = c("(A)", "(B)"),  label_x = 0.03, label_y = 1, label_size = 20)
figureS4

# exporting
# pdf(file = file.path("plots", "figureS4.pdf"), width = 15, height = 9.5)
# figureS4
# dev.off()







# all_ecopart %>%
#   filter(!patch_pair %in% c("wing-wing", "rectangle-rectangle")) %>%
#   ggplot(aes(metric, value, fill = patch_pair)) +
#   geom_boxplot()
# 
# all_ecopart %>%
#   group_by(patch_pair, metric) %>%
#   summarize(avg_value = mean(value))

# taking the absolute value for comparison
# all_ecopart_two <- all_ecopart_two %>%
#   mutate(abs_value = abs(value)) %>%
#   filter(!patch_pair %in% c("rectangle-rectangle", "wing-wing"))
# 
# # dividing into patch pair types for modeling
# c_r_ecopart <- all_ecopart_two %>%
#   filter(patch_pair == "Connected-Rectangular")
# c_w_ecopart <- all_ecopart_two %>%
#   filter(patch_pair == "Connected-Winged")
# w_r_ecopart <- all_ecopart_two %>%
#   filter(patch_pair == "Winged-Rectangular")
# 
# # models
# m.c_r_ecopart <- glmmTMB(abs_value ~ metric*time1 + (1|block1),
#               data = c_r_ecopart)
# summary(m.c_r_ecopart)
# 
# m.c_w_ecopart <- glmmTMB(abs_value ~ metric*time1 + (1|block1),
#                          data = c_w_ecopart)
# summary(m.c_w_ecopart)
# 
# m.w_r_ecopart <- glmmTMB(abs_value ~ metric*time1  + (1|block1),
#                          data = w_r_ecopart)
# summary(m.w_r_ecopart)
# 
# 
# all_ecopart_two$s.time <- as.numeric(scale(all_ecopart_two$time1))
# 
# colonizations <- all_ecopart_two %>%
#   filter(metric == "Colonization component")
# extinctions <- all_ecopart_two %>%
#   filter(metric == "Extinction component")
# 
# m1 <- glmmTMB(abs_value ~ metric + (1|block1),
#               data = s) 
# summary(m1)
# 
# m1.posthoc <- emmeans(m1, ~s.time*patch_pair)
# pairs(m1.posthoc, simple = "patch_pair")
# 
# 
# x <- all_ecopart_two %>%
#   select(-value) %>%
#   pivot_wider(names_from = "metric", values_from = abs_value) %>%
#   mutate(total = `Extinction component` - `Colonization component`)
# 
# 
# 
# 
# 
# 
# 
# 
# 
