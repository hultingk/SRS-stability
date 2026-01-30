

srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

# loading data from EDI 
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
edi_srs_data <- read.csv(file = srs_data_url)
edi_srs_data2 <- edi_srs_data %>%
  rowwise() %>% # grouping by rows
  mutate(occurance = sum(across(starts_with("X")), na.rm = T))%>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  select(!X1:X8) %>% # removing plot info - redundant
  filter(occurance == 1)

# loading species info from EDI
srs_dispersal_url <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.414.1&entityid=8de4a490a6ac6b05d2406c975d25b649"
edi_srs_dispersal <- read.csv(file = srs_dispersal_url)


# comparing species lists
new_spp_list <- srs_data %>%
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
  full_join(edi_spp_list, by = c("sppcode" = "SppCode")) %>%
  dplyr::select(sppcode, data1, data2, transplant) %>%
  filter(is.na(data1) | is.na(data2))



# blocks and years
edi_blocks <- edi_srs_data %>%
  count(EU, Year)
edi_blocks$EU <- substring(edi_blocks$EU, 3)
new_blocks <- srs_data %>%
  count(block, year)

new_blocks %>%
  anti_join(edi_blocks, by = c("block" = "EU",
                               "year" = "Year")) %>%
  View()
