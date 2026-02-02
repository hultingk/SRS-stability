librarian::shelf(tidyverse, summarytools)

srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

# loading data from EDI 
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
edi_srs_data <- read.csv(file = srs_data_url)
edi_srs_data2 <- edi_srs_data %>%
  rowwise() %>% # grouping by rows
  mutate(occurance = sum(across(starts_with("X")), na.rm = T))%>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  select(!X1:X8) %>% # removing plot info - redundant
  filter(occurance == 1)



# comparing species lists - taking out transplants and center patch to be consistent with 2019 data
new_spp_list <- srs_data %>%
  filter(patch_type != "Center") %>%
  filter(transplant != TRUE) %>%
  filter(year <= 2017) %>%
  filter(!(block == "54N" & year %in% c("2015", "2017"))) %>%
  count(sppcode, transplant) %>%
  mutate(data1 = "2024")


edi_spp_list <- edi_srs_data %>%
  count(SppCode) %>%
  mutate(data2 = "2019")
edi_spp_list2 <- edi_srs_data2 %>%
  count(SppCode) %>%
  mutate(data2 = "2019")

edi_spp_list %>%
  anti_join(edi_spp_list2, by = "SppCode")

all_spp_list <- new_spp_list %>%
  full_join(edi_spp_list2, by = c("sppcode" = "SppCode")) %>%
  dplyr::select(sppcode, data1, data2) %>%
  filter(is.na(data1) | is.na(data2))

write.csv(all_spp_list, file = file.path("data", "L0_original", "spp_compare.csv"))

# blocks and years
edi_blocks <- edi_srs_data %>%
  count(EU, Year)
edi_blocks$EU <- substring(edi_blocks$EU, 3)

new_blocks <-srs_data %>%
  filter(patch_type != "Center") %>%
  filter(transplant != TRUE) %>%
  filter(year <= 2017) %>%
  count(block, year)

new_blocks %>%
  anti_join(edi_blocks, by = c("block" = "EU",
                               "year" = "Year")) %>%
  View()
