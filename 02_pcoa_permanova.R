
librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB)

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(patch_type != "center") %>%
  #filter(!block %in% c("75W", "75E")) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  ))

# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode, soil_moisture, year_since_fire) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.factor(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)
srs_data_wider$year_since_fire <- as.numeric(srs_data_wider$year_since_fire)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year, soil_moisture, year_since_fire)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year", "soil_moisture", "year_since_fire"))


# Jaccard distance matrix
jaccard_dist <- vegdist(sp_info, method = "jaccard")
jaccard_dist_all <- as.matrix(jaccard_dist)
jaccard_dist_all_df <- as.data.frame(jaccard_dist_all)

jaccard_dist_all_df <- cbind(patch_info, jaccard_dist_all_df) # merge with patch info

# pcoa
pcoa_all <- pcoa(jaccard_dist)
pcoa_all_cmd <- cmdscale(jaccard_dist, eig=TRUE, add=FALSE)  ## Another way of doing the pcoa that gives same results#


# add axes to patch and time
pcoa_axes <- pcoa_all$vectors[,c(1,2)]
pcoa_axes <- cbind(patch_info, pcoa_axes)


#### species scores ####
# pcoa_all_cmd_scores <- add.spec.scores(pcoa_all_cmd, sp_info, method = 'cor.scores')
# 
# spscores<-as.data.frame(pcoa_all_cmd_scores$cproj)
# names<-rownames(spscores)
# rownames(spscores) <- NULL
# 
# spscores_df <- cbind(names,spscores)


#### permanova BY YEAR ####
# srs_wider_permanova <- srs_data_wider %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
#   filter(patch_rep %in% c("B", "C", "D", "E")) #%>%
#   #filter(!block %in% c("52", "57"))
# 
# 
# permanova_01_03 <- srs_wider_permanova %>%
#   filter(year %in% c("2001", "2002", "2003"))
# permanova_05_06 <- srs_wider_permanova %>%
#   filter(year %in% c("2005", "2006"))
# permanova_07 <- srs_wider_permanova %>%
#   filter(year %in% c("2007"))
# permanova_08_12 <- srs_wider_permanova %>%
#   filter(year %in% c("2008", "2009", "2010", "2011", "2012"))
# permanova_13 <- srs_wider_permanova %>%
#   filter(year %in% c("2013"))
# permanova_14_15 <- srs_wider_permanova %>%
#   filter(year %in% c("2014", "2015"))
# permanova_16 <- srs_wider_permanova %>%
#   filter(year %in% c("2016"))
# permanova_17 <- srs_wider_permanova %>%
#   filter(year %in% c("2017"))
# permanova_18_24 <- srs_wider_permanova %>%
#   filter(year %in% c("2018", "2019", "2020", "2021", "2022", "2023", "2024"))
# 
# loop_permanova <- function(df_wide) {
#   results <- adonis2(vegdist(df_wide[,8:332], method = "jaccard") ~ patch_type + soil_moisture + block, data = df_wide, by = "margin", method = "jaccard")
#   return(data.frame(results))
# }
# 
# 
# permanova_results_01_03 <- permanova_01_03 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2001:2003, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_05_06 <- permanova_05_06 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2005:2006, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_07 <- permanova_07 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2007, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_08_12 <- permanova_08_12 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2008:2012, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_13 <- permanova_13 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2013, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_14_15 <- permanova_14_15 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2014:2015, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_16 <- permanova_16 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2016, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_17 <- permanova_17 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2017, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results_18_24 <- permanova_18_24 %>%
#   group_by(year) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(year = rep(2018:2024, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(year, variable, R2)
# 
# permanova_results <- rbind(
#   permanova_results_01_03, permanova_results_05_06, permanova_results_07, permanova_results_08_12, permanova_results_13, permanova_results_14_15, permanova_results_16, permanova_results_17, permanova_results_18_24
# )
# 
# r2_permanova_plot <- permanova_results %>%
#   filter(!variable %in% c("Residual", "Total")) %>%
#   ggplot(aes(year, R2)) +
#   geom_point(size = 3) + 
#   geom_smooth(color = "black") +
#   facet_grid(cols = vars(variable)) +
#   theme_minimal(base_size = 22) +
#   xlab("Year") +
#   theme(panel.spacing = unit(1.5, "lines"),
#         plot.margin = margin(12,24,12,12)) +
#   theme(axis.text.x = element_text(angle = 60,  hjust=1))
# r2_permanova_plot
# 
# pdf(file = "r2_permanova_plot.pdf", width = 11, height = 7)
# r2_permanova_plot
# dev.off()


#### permanova BY TIME ####
time_permanova <- srs_data_wider %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  filter(patch_rep %in% c("B", "C", "D", "E")) #%>%
 # filter(!block %in% c("75E", "75W", "52", "57", "54N"))



time_permanova_0 <- time_permanova %>%
  filter(time %in% c("0"))
time_permanova_1_3 <- time_permanova %>%
  filter(time %in% c("1", "2", "3"))
time_permanova_4 <- time_permanova %>%
  filter(time %in% c("4"))
time_permanova_5 <- time_permanova %>%
  filter(time %in% c("5"))
time_permanova_6 <- time_permanova %>%
  filter(time %in% c("6"))
time_permanova_7_12 <- time_permanova %>%
  filter(time %in% c("7", "8", "9", "10", "11", "12"))
time_permanova_13 <- time_permanova %>%
  filter(time %in% c("13"))
time_permanova_14_15 <- time_permanova %>%
  filter(time %in% c("14", "15"))
time_permanova_16 <- time_permanova %>%
  filter(time %in% c("16"))
time_permanova_17 <- time_permanova %>%
  filter(time %in% c("17"))
time_permanova_18_24 <- time_permanova %>%
  filter(time %in% c("18", "19", "20", "21", "22", "23", "24"))


loop_permanova <- function(df_wide) {
  results <- adonis2(vegdist(df_wide[,9:327], method = "jaccard") ~ patch_type + soil_moisture + block, data = df_wide, by = "margin", method = "jaccard")
  return(data.frame(results))
}


## EXCLUDING time 0 and 4 -- only 57 and 52 sampled at those time points
# time_permanova_results_0 <- time_permanova_0 %>%
#   group_by(time) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(time = rep(0, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(time, variable, R2)

time_permanova_results_1_3 <- time_permanova_1_3 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(1:3, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)

# time_permanova_results_4 <- time_permanova_4 %>%
#   group_by(time) %>%
#   group_split() %>%
#   lapply(loop_permanova) %>%
#   bind_rows() %>%
#   rownames_to_column("variable") %>%
#   mutate(time = rep(4, each = 5)) %>%
#   mutate(variable = dplyr::case_when(
#     str_detect(variable, "patch_type") ~ "Patch Type",
#     str_detect(variable, "block") ~ "Block",
#     str_detect(variable, "soil_moisture") ~ "Soil Moisture",
#     str_detect(variable, "Residual") ~ "Residual",
#     str_detect(variable, "Total") ~ "Total"
#   )) %>%
#   select(time, variable, R2)


time_permanova_results_5 <- time_permanova_5 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(5, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)

time_permanova_results_6 <- time_permanova_6 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(6, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)

time_permanova_results_7_12 <- time_permanova_7_12 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(7:12, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)


time_permanova_results_13 <- time_permanova_13 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(13, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)


time_permanova_results_14_15 <- time_permanova_14_15 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(14:15, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)


time_permanova_results_16 <- time_permanova_16 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(16, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)



time_permanova_results_17 <- time_permanova_17 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(17, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)



time_permanova_results_18_24 <- time_permanova_18_24 %>%
  group_by(time) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(time = rep(18:24, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(time, variable, R2)

### putting it all together
permanova_results_time <- rbind(time_permanova_results_1_3, 
  time_permanova_results_5, time_permanova_results_6, time_permanova_results_7_12, 
  time_permanova_results_13, time_permanova_results_14_15, time_permanova_results_16,
  time_permanova_results_17, time_permanova_results_18_24
)

# removing residual and total -- don't care about those
permanova_results_time <- permanova_results_time %>%
  filter(!variable %in% c("Residual", "Total"))

# plotting
r2_permanova_time_plot <- permanova_results_time %>%
  ggplot(aes(time, R2)) +
  geom_point(size = 3) + 
  geom_smooth(color = "black") +
  facet_grid(cols = vars(variable)) +
  theme_minimal(base_size = 24) +
  xlab("Time since site creation (years)") +
  ylab("Explained variation") +
  theme(panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(12,24,12,12)) +
  theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.text.x= element_text(size = 24))
r2_permanova_time_plot

# pdf(file = file.path("plots", "r2_permanova.pdf"), width = 11, height = 7)
# r2_permanova_time_plot
# dev.off()


# patch_permanova <- permanova_results_time %>%
#   filter(variable == "Patch Type")
# block_permanova <- permanova_results_time %>%
#   filter(variable == "Block")
# block_permanova <- permanova_results_time %>%
#   filter(variable == "Block")
# 
# m1 <- glmmTMB(R2 ~ time*variable,
#               data = permanova_results_time)
# summary(m1)


#### plotting ####
pcoa_axes_plot <- pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  ))
pcoa_axes_plot$time <- as.numeric(pcoa_axes_plot$time)

plot_pcoa <- pcoa_axes_plot %>%
  filter(patch_type != "center") %>%
  filter(patch_rep %in% c("B", "C", "D")) %>%
  ggplot(aes(Axis.1, Axis.2, color = time)) +
  geom_point(size = 1.5) +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = time), linewidth = 1) +
  facet_grid(block~patch_type) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x= element_text(size = 24),
        strip.text.y= element_text(size = 24),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(color='Time',
       x = "PCoA 1",
       y = "PCoA 2") 
plot_pcoa

pdf(file = file.path("plots", "pcoa_plot.pdf"), width = 10, height =15)
plot_pcoa
dev.off()

# pdf(file = file.path("plots", "time9_53N_pcoa_plot.pdf"), width = 12, height =6)
# plot_pcoa
# dev.off()

plot_pcoa_53N <- pcoa_axes_plot %>%
  filter(block %in% c("53N")) %>%
  filter(patch_type != "center") %>%
  filter(patch_rep %in% c("B", "C", "D")) %>%
  ggplot(aes(Axis.1, Axis.2, color = time)) +
  geom_point(size = 2) +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = time), linewidth = 1) +
  facet_grid(block~patch_type) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x= element_text(size = 26),
        strip.text.y= element_text(size = 26),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(color='Time',
       x = "PCoA 1",
       y = "PCoA 2") 
plot_pcoa_53N



