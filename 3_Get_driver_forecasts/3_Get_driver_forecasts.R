#Get driver forecasts of meteorology, water temperature, and inflow
#Author: Mary Lofton
#Date: 30MAR22

#Tasks to do in this script: ####

#1. Retrieve NOAA GEFS ensemble forecasts for 00:00 hrs from 2021-01-01 to 2021-12-31
#   These will come from the flare-forecast.org S3 bucket.
#   a. Aggregate each variable to a daily average
#   b. Format for use in Jake's EnKF code (see load_noaa_forecast function developed
#      by Tadhg for EDDIE and use as baseline; this may even do the daily averaging)
#2. Retrieve and format water temperature forecasts from FLARE for FCR from 2021-01-01
#   to 2021-12-31; have pinged Quinn to see where these are currently stored as of 30MAR22
#3. Create discharge ensemble forecasts; this is going to take a bit of detective work.
#   Whitney did this for her chl-a forecasting paper. The function is written in the
#   create_discharge_forecast_ensembles.R script but it is poorly annotated. Your best
#   bet is to see the source code where Whitney applies the function
#   (https://github.com/wwoelmer/FLARE_AR_CHLA/blob/EcoApps_MS_Feb2022/Rscripts/run_arima_any_timestep.R)
#   and try to back out from there what each of the arguments should be. I would guess
#   this will take you the better part of a day.

#SET-UP####

#clear environment
rm(list = ls())

#load packages
#install.packages('pacman')
pacman::p_load(tidyverse, lubridate, aws.s3, ncdf4)

#set environment to retrieve from proper S3 bucket
Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

#source helper functions
source("./0_Function_library/format_noaa_forecast.R")
source("./0_Function_library/create_discharge_forecast_ensembles.R")

#RETRIEVE AND FORMAT NOAA FORECASTS####
#set series of dates over which you want to retrieve forecasts
model_dates = seq.Date(from = as.Date("2021-01-21"), to = as.Date("2021-12-31"), by = "days")


#NOTE to self (MEL): you are storing the .nc files of each ensemble member on your
#external hard drive due to limited space on your PC laptop; however, formatted
#output has been copied over to the phyto-data-assimilation RProject on your laptop
for (i in seq_along(model_dates)){
  
  #assign bucket prefix for each date
  prefix <- paste0("noaa/NOAAGEFS_1hr/fcre/",model_dates[i],"/00")
  
  #retrieve filenames of ensemble member files
  tryCatch({
    fc.files <- get_bucket_df(bucket = "drivers", prefix = prefix)

  for(j in 1:nrow(fc.files)){
    
    #download and save each ensemble member file
    #MEL using external hard drive here b/c my laptop is old and sad
    save_object(
      object = fc.files[j,1],
      bucket = "drivers",
      file = file.path(paste0("F:/phyto_DA/",fc.files[j,1])))
  }

  #once all ensembles are downloaded, load them in and re-format to daily timestep 
  fc <- format_noaa_forecast(siteID = "fcre",
                           start_date = as.character(model_dates[i]),
                           my_directory = "F:/phyto_DA/")
  #I *think* the result of the this function should be directly compatible as input
  #to the get_drivers function of in my most recent version of the module7 Rmd
  
  write.csv(fc, file = file.path(paste0("F:/phyto_DA/formatted_noaa/met_",model_dates[i],".csv")),row.names = FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

#note that three days are missing: 21-22 JAN and one more...
#placeholder - may decide to come back here and aggregate the daily files to one file

#2. RETRIEVE AND FORMAT WATER TEMPERATURE FORECASTS####

#Waiting for the team to find the 2021 FLARE water temperature forecasts at FCR.

#3. USE NOAA FORECASTS TO CREATE DISCHARGE FORECASTS####

#set up arguments for create_discharge_forecast_ensembles file
met_file_names <- list.files("./00_Data_files/formatted_noaa") 
fc_dates = seq.Date(from = as.Date("2021-01-01"), to = as.Date("2022-02-04"), by = "days")


working_directory <- "./00_Data_files/discharge_forecasts/"

#don't need to run this once it's written to file
# inflow_file1 <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
#   filter(Date %in% fc_dates) %>%
#   select(Date, Flow_cms) %>%
#   rename(time = Date,
#          FLOW = Flow_cms)
# write.csv(inflow_file1,file = "./00_Data_files/inflow_for_forecasts.csv",row.names = FALSE)


create_discharge_forecast_ensembles(working_directory = working_directory,
                                    start_forecast_step = 2,
                                    inflow_file1 = "./00_Data_files/inflow_for_forecasts.csv",
                                    met_file_names = met_file_names,
                                    inflow_process_uncertainty = TRUE)

#check how they look
inf_fc_files <- list.files("./00_Data_files/discharge_forecasts")
test <- read_csv(file.path(paste0("./00_Data_files/discharge_forecasts/",inf_fc_files[361])))
