library(tidyverse)
library(here)

files_to_combine <- read_csv(here("sim_output/to_combine"))

all_dat <- data.frame() 

for(f in files_to_combine){
  this_file <- read_csv(here(paste0("sim_output/",f)))
  all_dat <- bind_rows(all_dat, this_file)
}

write_csv(all_dat,file = here("data/variability_data.csv"))