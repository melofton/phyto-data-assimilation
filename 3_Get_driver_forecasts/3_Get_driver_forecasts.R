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
install.packages('pacman')
pacman::p_load(tidyverse, lubridate, aws.s3, ncdf4)

#set environment to retrieve from proper S3 bucket
Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

#source helper functions
source("./0_Function_library/format_noaa_forecast.R")

#set series of dates over which you want to retrieve forecasts
model_dates = seq.Date(from = as.Date("2021-01-21"), to = as.Date("2021-12-31"), by = "days")

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
