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
pacman::p_load(tidyverse, lubridate, aws.s3, ncdf4, ncdf4.helpers)

#set environment to retrieve from proper S3 bucket
Sys.setenv("AWS_DEFAULT_REGION" = "s3",
           "AWS_S3_ENDPOINT" = "flare-forecast.org")

#source helper functions
source("./0_Function_library/format_noaa_forecast.R")
source("./0_Function_library/create_discharge_forecast_ensembles.R")

#define simple functions
split_params = function (x, sep, n) {
  # Splits string into list of substrings separated by 'sep'.
  # Returns nth substring.
  x = strsplit(x, sep)[[1]][n]
  
  return(x)
}

#1. RETRIEVE AND FORMAT NOAA FORECASTS####
#set series of dates over which you want to retrieve forecasts
model_dates = seq.Date(from = as.Date("2021-01-01"), to = as.Date("2021-12-31"), by = "days")


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
model_dates = seq.Date(from = as.Date("2021-01-01"), to = as.Date("2021-12-31"), by = "days")

for (i in seq_along(model_dates)){
  
  tryCatch({ #because we know there are 4 days missing so this will bonk 4 times
    
  #create object key name within bucket
  object = paste0("fcre/fcre-",model_dates[i],"-marylofton.nc")
  
  #create filename where object will be written locally
  filename = paste0("C:/Users/Mary Lofton/Downloads/fcre-",model_dates[i],"-marylofton.nc")
  
  #retrieve object from bucket and write to file locally
  save_object(
    object = object,
    bucket = "forecasts",
    file = file.path(filename))
  
  #open nc file
  nc <- nc_open(filename)
  
  #retrieve water temperature forecast
  #the dim names of the temp variable are time, depth, ensemble (18,28,100)
  #time includes yesterday, today, and 16 days into the future
  #depth includes 28 depths for FCR, ~ every third of a meter
  #there are 100 ensemble members
  #we are indexing to retrieve the forecast starting from today (time dimension = 2),
  #for a depth of 1.6 m (depth dimension = 6), for all ensemble members
  #that's how we're getting to start = c(2,6,1) and count = c(-1,1,-1)
  wt <- ncvar_get(nc, attributes(nc$var)$names[1], start = c(2,6,1), count = c(-1,1,-1))
  
  #retrieve and format forecast dates
  t <- ncdf4::ncvar_get(nc,'time')
  full_time <- as.POSIXct(t,
                          origin = '1970-01-01 00:00.00 EST',
                          tz = "UTC")
  full_time_day <- lubridate::as_date(full_time)
  
  #format forecast output into long form to match other driver forecasts
  fc <- data.frame(wt)
  fc$fc_date <- full_time_day[2:18]
  ens <- c(1:100)
  colnames(fc)[1:100] <- paste0("ens_",ens)
  
  fc <- fc %>%
    gather(ens_1:ens_100, key = "ensemble_member", value = "Temp_C") %>%
    mutate(ensemble_member = sapply(
      X = ensemble_member,
      FUN = split_params,
      sep = '_',
      n = 2)) %>%
    add_column(issue_date = model_dates[i])
  fc <- fc[,c(4,1,2,3)]
  
  #write forecast to file
  write.csv(fc, file = file.path(paste0("./00_Data_files/FLARE_forecasts/wt_",model_dates[i],".csv")),row.names = FALSE)
  
  #close nc file
  ncdf4::nc_close(nc)
  
  #delete local nc file
  file.remove(filename)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #print error if occurs
  
} #end loop

test <- read_csv("./00_Data_files/FLARE_forecasts/wt_2021-01-01.csv")
#looks good

check <- list.files("./00_Data_files/FLARE_forecasts/")
#ha ok we are now down to 360 days...hmmm...ok whatever I am not too concerned
#we can circle back to this once we are writing it all up if that happens


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
