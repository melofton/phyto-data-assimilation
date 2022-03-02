##DATA WRANGLING TO COLLATE OBSERVATION INPUT FILE FOR REGULARIZED REGRESSION
##Author: Mary Lofton
##Date: 02MAR22

##To-do in this script:####
#1. Pull in all relevant cleaned data files
#2. Collate data
#3. Write to file

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

##READ IN DATA####
phytos <- read_csv("./1_Data_wrangling/collated_phyto_data.csv")
wt <- read_csv("./1_Data_wrangling/collated_water_temp_data.csv") %>%
  select(Temp_C)
met <- read_csv("./1_Data_wrangling/collated_met_data.csv") %>%
  select(-Date,-daily_WindDir_degrees,-daily_PAR_Average_umol_s_m2)
inf <- read_csv("./1_Data_wrangling/collated_inflow_data.csv") %>%
  select(Flow_cms)

##COLLATE DATA####
dat <- bind_cols(phytos, wt, met, inf)

##WRITE TO FILE####
write.csv(dat, "./1_Data_wrangling/collated_obs_data.csv",row.names = FALSE)
