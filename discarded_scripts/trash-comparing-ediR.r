librarian::shelf(tidyverse, summarytools)

# loading cleaned 2024 data 
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

# loading data from EDI 
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
edi_srs_data <- read.csv(file = srs_data_url)

# formating EDI data same way as 2024 data (each row is occurrence of a sp in one patch in one year)
edi_srs_data2 <- edi_srs_data %>%
  rowwise() %>% # grouping by rows
  mutate(occurance = sum(across(starts_with("X")), na.rm = T))%>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  select(!X1:X8) %>% # removing plot info - redundant
  filter(occurance == 1)
edi_srs_data2$EU <- substring(edi_srs_data2$EU, 3)



# comparing species lists
# making list of species in new data
new_spp_list <- srs_data %>%
  filter(patch_type != "Center") %>% # taking out transplants and center patch to be consistent with EDI data
  filter(transplant != TRUE) %>%
  filter(year <= 2017) %>% # only comparing same data as EDI - nothing past 2017
  filter(!(block == "54N" & year %in% c("2015", "2017"))) %>% # removing 2015 and 2017 surveys in 54N - not in EDI data
  count(sppcode, dispersal_mode) %>% 
  rename(new_dispersal = dispersal_mode) %>%
  mutate(data1 = "new data")

# making list of species in new data
edi_spp_list2 <- edi_srs_data2 %>%
  count(SppCode) %>%
  mutate(data2 = "edi data")

# joining lists together and filtering for species that occur in one list but not the other
all_spp_list <- new_spp_list %>%
  full_join(edi_spp_list2, by = c("sppcode" = "SppCode")) %>%
  dplyr::select(sppcode, data1, data2) %>%
  filter(is.na(data1) | is.na(data2))

# write.csv(all_spp_list, file = file.path("data", "L0_original", "spp_compare.csv"))

# comparing blocks and years
# EDI blocks
edi_blocks <- edi_srs_data %>%
  count(EU, Year)
edi_blocks$EU <- substring(edi_blocks$EU, 3)

# new data blocks
new_blocks <-srs_data %>%
  filter(patch_type != "Center") %>%
  filter(transplant != TRUE) %>%
  filter(year <= 2017) %>%
  count(block, year)

# joining and comparing
new_blocks %>%
  anti_join(edi_blocks, by = c("block" = "EU",
                               "year" = "Year")) %>%
  View()




#### comparing dispersal modes ####
# importing edi data
edi_dispersal_url  <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/8de4a490a6ac6b05d2406c975d25b649" 
edi_dispersal <- read.csv(file = edi_dispersal_url)
edi_dispersal <- edi_dispersal %>%
  rename(EDI_dispersal = DispersalMode)

disperal_compare <- new_spp_list %>%
  full_join(edi_dispersal, by = c("sppcode" = "SppCode")) %>%
  dplyr::select(sppcode, new_dispersal, EDI_dispersal) 

# looking at where one dispersal mode isn't in the other -- these correspond to the spp list above except for HOULON, LECLEG, and PRESER, which are present in the EDI data but do not occur in any patches
disperal_compare %>%
  filter(is.na(new_dispersal) | is.na(EDI_dispersal)) %>%
  View()

# no observations where dispersal modes are inconsistent
disperal_compare %>%
  filter(new_dispersal != EDI_dispersal)


#### comparing actual data ####
# making object of species that we know are going to differ between datasets
inconsistent_sppcode <- all_spp_list$sppcode

# calculating # of species per year per patch, excluding inconsistent species
new_spp_rich <- srs_data %>%
  filter(patch_type != "Center") %>% # taking out transplants and center patch to be consistent with EDI data
  filter(transplant != TRUE) %>%
  filter(year <= 2017) %>% # only comparing same data as EDI - nothing past 2017
  filter(!(block == "54N" & year %in% c("2015", "2017"))) %>%
  filter(!sppcode %in% inconsistent_sppcode) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("Connected") ~ "C", 
    patch_type %in% c("Rectangular") ~ "R", 
    patch_type %in% c("Winged") ~ "W", 
    .default = patch_type
  )) %>%
  mutate(unique_id2 = paste(block, patch_type, sep = "-")) %>%
  count(unique_id2, year) %>%
  rename(new_rich = n)

# calculating # of species per year per patch, excluding inconsistent species
edi_spp_rich <- edi_srs_data2 %>%
  filter(!SppCode %in% inconsistent_sppcode) %>%
  mutate(unique_id = paste(EU, PatchType, sep = "-")) %>%
  count(unique_id, Year) %>%
  rename(edi_rich = n)

# checking to see when the species numbers aren't the same
new_spp_rich %>%
  full_join(edi_spp_rich, by = c("unique_id2" = "unique_id",
                                 "year" = "Year")) %>%
  filter(new_rich != edi_rich) %>%
  View()


