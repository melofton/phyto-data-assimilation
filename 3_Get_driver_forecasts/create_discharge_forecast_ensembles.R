#' Function to create inflow discharge ensemble forecast
#' 
#' @param full_time_day_local vector of dates over which you want to create discharge
#' forecast ensembles; format is chr string in yyyy-mm-dd
#' @param working_directory where forecast ensemble output is written
#' @param input_file_tz time zone of the met file inputs (?) DON'T CHANGE THIS
#' @param start_forecast_step should be set to 2 b/c the first day in the full_time_day_local
#' vector is "today" so we want to start forecasts at second timestep
#' @param inflow_file1 .csv file with time, FLOW, TEMP and time is the date in yyyy-mm-dd;
#' FLOW is discharge in cms; TEMP is water temperature of inflow in C
#' @param local_tzone is the same as input_file_tz
#' @param met_file_names GLM compatible met files with columns and units accordingly; each
#' file is the forecast from one NOAA ensemble member and you need a separate file for each 
#' ensemble member/date combination that you want discharge forecasts for; what it 
#' ultimately uses is AirTemp and Rain for linear model of discharge
#' @param forecast_days number of days in the future that you want discharge forecasts for
#' @param inflow_process_uncertainty TRUE or FALSE as to whether to include process
#' uncertainty; the distribution for process uncertainty is specified w/in function
#' 
#' If need to look in more detail, run the automated forecast scripts (6) in order first
#' to set everything up and then use the run_arima_any_timestep when you get to that script
#' to step through the function

create_discharge_forecast_ensembles <- function(full_time_day_local,
                                       working_directory,
                                       input_file_tz = 'EST5EDT', 
                                       start_forecast_step,
                                       inflow_file1,
                                      # inflow_file2,
                                      # outflow_file1,
                                       #chemistry_file,
                                       local_tzone,
                                       met_file_names,
                                       forecast_days,
                                       inflow_process_uncertainty){
  
  full_time_day_local <- as_date(full_time_day_local)
  
  inflow <- read_csv(inflow_file1)
  #wetland <- read_csv(inflow_file2)

  if(include_wq){
    wq_names_tmp <- wq_names[which(wq_names %in% names(inflow))]
  }else{
    wq_names_tmp <- NULL
  }
  
  curr_all_days <- NULL
  
  col_types <- cols(
    time = col_datetime(format = ""),
    ShortWave = col_double(),
    LongWave = col_double(),
    AirTemp = col_double(),
    RelHum = col_double(),
    WindSpeed = col_double(),
    Rain = col_double(),
    Snow = col_double())
  
  for(m in 1:length(met_file_names)){
    curr_met_daily <- read_csv(met_file_names[m],
                               col_types = col_types) %>% 
      mutate(time = as_date(time)) %>% 
      group_by(time) %>% 
      dplyr::summarize(Rain = mean(Rain, na.rm = TRUE),
                AirTemp = mean(AirTemp, na.rm = TRUE)) %>% 
      mutate(ensemble = m) %>% 
      mutate(AirTempMean = roll_mean(AirTemp, n = 5, align = "right",fill=NA),
             RainMean = roll_mean(Rain, n = 5, align = "right",fill=NA),
             AirTemp_lag1 = lag(AirTemp, 1),
             Rain_lag1 = lag(Rain, 1))
    
    curr_all_days <- rbind(curr_all_days,curr_met_daily)
  }
  
  
  tmp <- left_join(curr_all_days, inflow, by = "time")
  
  forecasts_days <- full_time_day_local[start_forecast_step:length(full_time_day_local)]
  if(use_future_inflow == FALSE || start_forecast_step == length(full_time_day_local)){
    forecasts_days <- NULL 
  }
  
  tmp <- tmp %>%
    mutate(forecast = ifelse(time %in% forecasts_days, 1, 0),
           TEMP = ifelse(forecast == 1, NA, TEMP),
           FLOW = ifelse(forecast == 1, NA, FLOW))
  
  if(inflow_process_uncertainty == TRUE){
    inflow_error <- rnorm(nrow(tmp), 0, 0.00965)
    temp_error <- rnorm(nrow(tmp), 0, 0.943)
  }else{
    inflow_error <- rep(0.0, nrow(tmp))
    temp_error <- rep(0.0, nrow(tmp))
  }
  
  for(i in 1:nrow(tmp)){
    if(tmp$forecast[i] == 0 & is.na(tmp$FLOW[i]) & !include_wq){
      list(tmp[i, c("FLOW", "TEMP",wq_names_tmp)]  <- inflow %>% 
        filter(time < full_time_day_local[start_forecast_step]) %>% 
        mutate(doy = yday(time)) %>% 
        filter(doy == yday(tmp$time[i])) %>% 
        summarize_at(.vars = c("FLOW", "TEMP", wq_names_tmp), mean, na.rm = TRUE)) %>% 
        unlist()
    }
    
    if(tmp$forecast[i] == 1){
      tmp$FLOW[i] = 0.0010803 + 0.9478724 * tmp$FLOW[i - 1] +  0.3478991 * tmp$Rain_lag1[i] + inflow_error[i]
      tmp$TEMP[i] = 0.20291 +  0.94214 * tmp$TEMP[i-1] +  0.04278 * tmp$AirTemp_lag1[i] + temp_error[i]
      if(include_wq){
        tmp[i, c(wq_names_tmp)] <- inflow %>% 
          filter(time < full_time_day_local[start_forecast_step]) %>% 
          mutate(doy = yday(time)) %>% 
          filter(doy == yday(tmp$time[i])) %>% 
          summarize_at(.vars = c(wq_names_tmp), mean, na.rm = TRUE) %>% 
          unlist()
      }
    }
  }
  
  tmp <- tmp %>% 
    mutate(FLOW = ifelse(FLOW < 0.0, 0.0, FLOW))
  

  for(i in 1:n_distinct(tmp$ensemble)){
    tmp2 <- tmp %>% 
      filter(ensemble == i) %>% 
      mutate(SALT = 0.0) %>% 
      dplyr::select(time, FLOW, TEMP, SALT) %>% 
      mutate_at(vars(c("FLOW", "TEMP", "SALT", wq_names_tmp)), funs(round(., 4)))
      
   if(i < 10){
      file_name_discharge <- paste0('inflow_forecast_ensemble_', "0", i)
    }else{
      file_name_discharge <- paste0('inflow_forecast_ensemble_', i) 
    }
    write_csv(x = tmp2,
              path = paste0(working_directory,"/", file_name_discharge, '.csv'),
              quote_escape = "none")

    tmp2 <- tmp2 %>% 
     dplyr::select(time, FLOW)
    

  }
  
  return(list(tmp2 = tmp2))
}
