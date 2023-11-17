#Sit assessment specific script
#Read in latest linelist file and set target dates for modelling

#currently manually specifying date range and linelist file name.
#Update to automatically read in latest linelist file and assign target dates
source("R/summarise_linelist.R")

linelist_file <- "data-raw/processed_linelist_20231109.rds"
linelist <- readRDS(linelist_file)

local_summary <- summarise_linelist(linelist,
                                    import_status_option = 'local')

# make target dates
target_dates <- as.character(
    seq.Date(as.Date("2023-06-01"),
             as.Date("2023-11-08"),
             by = "day"))

