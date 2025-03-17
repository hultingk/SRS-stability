# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, vegan)

# loading data from EDI
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
srs_data <- read.csv(file = srs_data_url)

srs_data <- srs_data %>%
  rowwise() %>%
  mutate(occurance = sum(across(starts_with("X")), na.rm = T))

srs_data <- srs_data %>%
  select(!X1:X8) %>%
  filter(occurance == 1)

