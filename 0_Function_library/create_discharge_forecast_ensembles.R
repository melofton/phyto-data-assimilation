#' Function to create inflow discharge ensemble forecast
#' Adapted from https://github.com/wwoelmer/FLARE_AR_CHLA/blob/EcoApps_MS_Feb2022/Rscripts/create_discharge_forecast_ensembles.R
#' My understanding of authorship provenance is the R. Quinn Thomas originally developed the function for FLARE,
#' it was subsequently adapted by Whitney Woelmer for her chl-a ARIMA forecasting project,
#' and now MEL is adapting for use in phyto community forecasting using regularized regression
#' Last updated: 05APR22
#' Major changes to accommodate that we don't really need all the GLM file formatting for use
#' in regularized regression, which simplifies things somewhat
#' Also don't need forecasted inflow water temperature, just flow
#' 
#' @param working_directory where discharge forecast ensemble output is written
#' @param start_forecast_step should be set to 2 b/c the first day in each forecast file
#' vector is "today" so we want to start forecasts at second timestep
#' @param inflow_file1 .csv file with time and FLOW; time is the date in yyyy-mm-dd;
#' FLOW is discharge in cms
#' @param met_file_names names of met files you will read in to convert to discharge forecast
#' for this version of the function, this file is a re-formatted version of NOAA GEFS 35-day 
#' forecasts with the following columns:
#' [1] "fc_date": forecast date                     "air_pressure"                             
#' [3] "air_temperature"                           "cloud_area_fraction"                      
#' [5] "precipitation_flux"                        "relative_humidity"                        
#' [7] "specific_humidity"                         "surface_downwelling_longwave_flux_in_air" 
#' [9] "surface_downwelling_shortwave_flux_in_air" "wind_speed"                               
#' [11] "issue_date": date forecast was issued    "ensemble_member" : ranges from 1-30
#' all met vars are in NOAA units aggregated to daily
#' @param inflow_process_uncertainty TRUE or FALSE as to whether to include process
#' uncertainty; the distribution for process uncertainty is specified w/in function
#' 


create_discharge_forecast_ensembles <- function(full_time_day_local,
                                       working_directory,
                                       start_forecast_step,
                                       inflow_file1,
                                       met_file_names,
                                       inflow_process_uncertainty){
  
  #read in inflow file
  inflow <- read_csv(inflow_file1)

  #setting null object for output of for-loop below
  curr_all_days <- NULL
  
  #loop through all met files, read in, convert to GLM units for regression equation,
  #and collate into a single data frame
  for(m in 1:length(met_file_names)){ 
    
    curr_met_daily <- read_csv(file.path(paste0("./00_Data_files/formatted_noaa/",met_file_names[m]))) %>%
      #select columns we need
      select(fc_date, precipitation_flux, issue_date, ensemble_member) %>%
      #convert NOAA units to GLM units
      mutate(Rain = precipitation_flux * (60 * 60 * 24)/1000) %>%
      group_by(ensemble_member) %>%
      #get lags for regression equation below
      mutate(Rain_lag1 = lag(Rain, 1)) %>%
      #renaming fc_date column to be consistent w/ inflow file for joining purposes
      rename(time = fc_date) %>%
      #again, select only columns we need
      select(time, Rain_lag1, issue_date, ensemble_member) %>%
      ungroup()
    
    #other NOAA to GLM unit conversion for reference in case want to build out again
    # GLM$AirTemp <- NOAA$AirTemp - 273.15
    # GLM$RelHum <- NOAA$RelHum * 100
    
    
    #binding "today" to the other days; yeah so here it seems like we could use this loop
    #over multiple days but then why are we naming ensemble "m"???? makes no sense to me
    curr_all_days <- rbind(curr_all_days,curr_met_daily)
  }
  
  #creating combined file where daily rain lags is joined to daily mean inflow
  tmp <- left_join(curr_all_days, inflow, by = "time")

  
  #determining whether to populate output file with a forecast (below) or observed inflow data
  #the issue date for each forecast file should be populated with observed inflow from that day
  #which will be used as the starting point in the autoregressive equation below to generate
  #an inflow forecast
  tmp <- tmp %>%
    mutate(forecast = ifelse(time == issue_date, 0, 1),
           FLOW = ifelse(forecast == 1, NA, FLOW))
  
  #generating error to add to regression output if want to take into account process error
  #of the model
  if(inflow_process_uncertainty == TRUE){
    inflow_error <- rnorm(nrow(tmp), 0, 0.00965)
    #temp_error <- rnorm(nrow(tmp), 0, 0.943) #just in case we want to build this out again
  }else{
    inflow_error <- rep(0.0, nrow(tmp))
  }
  
  #so reminder here that each row of tmp is one ensemble member from one day of a forecast
  #issued on a particular day
  print("iterating through tmp df")
  for(i in 1:nrow(tmp)){
    #checking to see if we want a forecast for that day or not; if not,
    #and the observed inflow is missing
    #we subset and summarize the inflow data file using doy from other years
    #and populate the row with mean inflow on that doy from other years so we have a starting
    #point for the forecast
    # if(tmp$forecast[i] == 0 & is.na(tmp$FLOW[i])){
    #   list(tmp[i, c("FLOW")]  <- inflow %>% 
    #     filter(time < tmp[i,"time"]) %>% 
    #     mutate(doy = yday(time)) %>% 
    #     filter(doy == yday(tmp$time[i])) %>% 
    #     summarize_at(.vars = c("FLOW"), mean, na.rm = TRUE)) %>% 
    #     unlist()
    # }
    #if we do want a forecast, we use the regression equation below (+ error) to generate
    #forecasted inflow and inflow temperature
    if(tmp$forecast[i] == 1){
      tmp$FLOW[i] = 0.0010803 + 0.9478724 * tmp$FLOW[i - 1] +  0.3478991 * tmp$Rain_lag1[i] + inflow_error[i]
      #keeping the equation below just in case
      #tmp$TEMP[i] = 0.20291 +  0.94214 * tmp$TEMP[i-1] +  0.04278 * tmp$AirTemp_lag1[i] + temp_error[i]
      
    }
  }
  print("finished iterating thru tmp df")
  #if FLOW is negative, set it to 0 - makes sense
  tmp <- tmp %>% 
    mutate(FLOW = ifelse(FLOW < 0.0, 0.0, FLOW))
  
  #ok now we are iterating through forecast issue dates, rounding the estimated inflow,
  #and writing to file based on issue date
  
  #this is a little redundant b/c we could just write the collated file, and I may come
  #back and edit the function to do just that, but for now keeping as is to match
  #how the met forecasts for this project are organized
  issue_dates <- unique(tmp$issue_date)
  
  for(i in 1:n_distinct(tmp$issue_date)){

    tmp2 <- tmp %>% 
      filter(issue_date == issue_dates[i]) %>%
      rename(fc_date = time) %>%
      dplyr::select(fc_date, issue_date, ensemble_member,FLOW) %>% 
      mutate(FLOW = round(FLOW,4))
   
      file_name_discharge <- paste0('inflow_', issue_dates[i])
   
    #writing to file
    write.csv(tmp2,file = file.path(paste0(working_directory,"/", file_name_discharge, '.csv')),row.names = FALSE)

  }
  #yup
  
}
