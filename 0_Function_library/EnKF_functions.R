#Source code for functions to run EnKF
#Author: adapted by Mary Lofton using code from Jake Zwart provided to 2019 GLEON ecological forecasting workshop
#Date: 15APR22

rm(list = ls())
pacman::p_load(tidyverse, lubridate, PerformanceAnalytics)

##' retreive the model time steps based on start and stop dates and time step ----
#'
#' @param model_start model start date in date class
#' @param model_stop model stop date in date class
#' @param time_step model time step, defaults to daily timestep
get_model_dates = function(model_start, model_stop, time_step = 'days'){
  
  model_dates = seq.Date(from = as.Date(model_start), to = as.Date(model_stop), by = time_step)
  
  return(model_dates)
}

model_dates <- get_model_dates(model_start = "2021-01-01",
                               model_stop = "2021-12-31",
                               time_step = 'days')

##' get observational data ----
#' function to read in data 
#'
#' @param obs_file filepath and filename containing observational data
#' @param assimilation_expt choose from c("FP","FP_EXO","FP_CTD","all");
#' determines which data streams to include
#' @param missing_data_expt choose from TRUE or FALSE; determines whether or
#' not to drop every other FP observation

get_obs_df <- function(obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                       assimilation_expt = "all",
                       missing_data_expt = FALSE){
  
  obs <- read_csv(obs_file) 
  
  #check for skewness and log-transform if needed
  for(i in 2:ncol(obs)){
    if(abs(skewness(obs[,i],na.rm = TRUE)) > abs(skewness(log(obs[,i]+0.0001),na.rm = TRUE))){
      obs[,i] <- log(obs[,i]+0.0001)
    }
  }
  
  obs <- obs %>% 
    filter(year(Date) %in% c(2021)) %>%
    rename(datetime = Date)
  obs <- obs[,c(1,25,28,31,34,8,9,2:4)]
  
  if(assimilation_expt == "FP"){
    obs[,c("daily_EXOChla_ugL_1", "daily_EXOBGAPC_ugL_1", "CTDChla_ugL_above", "CTDChla_ugL", "CTDChla_ugL_below")] <- NA
  }
  
  if(assimilation_expt == "FP_EXO"){
    obs[,c("daily_EXOChla_ugL_1", "daily_EXOBGAPC_ugL_1")] <- NA
  }
  
  if(assimilation_expt == "FP_CTD"){
    obs[,c("CTDChla_ugL_above", "CTDChla_ugL", "CTDChla_ugL_below")] <- NA
  }
  
  if(missing_data_expt == TRUE){
    #place holder for later but honestly I don't think I'm going to get to this by JASM
  }
  
  return(obs)
  
}

obs_df <- get_obs_df(obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                     assimilation_expt = "all",
                     missing_data_expt = FALSE)

##' turn observation dataframe into matrix ----
#'
#' @param obs_df observation data frame
#' @param model_dates dates over which you're modeling
#' @param n_step number of model time steps
#' @param n_states number of states we're updating in data assimilation routine
#' @param states character string vector of state names in obs_file

get_obs_matrix = function(obs_df = obs_df, 
                          model_dates = model_dates, 
                          n_step = length(model_dates), 
                          n_states = 4, 
                          states = c("GreenAlgae_ugL",
                                     "Bluegreens_ugL",
                                     "BrownAlgae_ugL",
                                     "MixedAlgae_ugL")){
  
  # need to know location and time of observation
  
  obs_df_filtered = obs_df %>%
    dplyr::filter(as.Date(datetime) %in% model_dates) %>%
    mutate(date = as.Date(datetime)) %>%
    mutate(date_step = which(model_dates %in% date))
  
  obs_df_filtered <- obs_df_filtered[,c("date",states,"date_step")]
  
  obs_matrix = array(NA, dim = c(n_states, 1, n_step))
  
  for(i in 1:n_states){
    for(j in obs_df_filtered$date_step){
      obs_matrix[i, 1, j] = dplyr::filter(obs_df_filtered,
                                          date_step == j) %>%
        pull(states[i])
    }}
  
  return(obs_matrix)
}

obs <- get_obs_matrix(obs_df = obs_df, 
                      model_dates = model_dates, 
                      n_step = length(model_dates), 
                      n_states = 4, 
                      states = c("GreenAlgae_ugL",
                                 "Bluegreens_ugL",
                                 "BrownAlgae_ugL",
                                 "MixedAlgae_ugL"))

ancillary_obs <- get_obs_matrix(obs_df = obs_df, 
                                model_dates = model_dates, 
                                n_step = length(model_dates), 
                                n_states = 5, 
                                states = c("daily_EXOChla_ugL_1",
                                           "daily_EXOBGAPC_ugL_1",
                                           "CTDChla_ugL_above",
                                           "CTDChla_ugL",
                                           "CTDChla_ugL_below"))



##' observation error matrix, should be a square matrix where ----
#'   col & row = the number of states and params for which you have observations
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_step number of model timesteps
#' @param FP_sd vector of state observation standard deviation; assuming sd is constant through time; in our case sd of FP observations on log scale
#' @param param_sd vector of parameter observation standard deviation; assuming sd is constant through time
#' @param EXO_sd vector of standard deviations of EXO model predictions for each state
#' @param CTD_sd vector of standard deviations of CTD model predictions for each state
#' @param EXO_CTD_sd vector of standard deviations of EXO + CTD model predictions for each state

get_obs_error_matrix = function(n_states = 4, 
                                n_params_obs = 0, 
                                n_step = length(model_dates), 
                                FP_var = c(0.009,0.01,0.04,0.90), #remember log scale; estimated this using var(log(rnorm(1000,x,0.5))),
                                #where x is the mean of each PG from 2014-2020; for green algae x = 5.6, for bluegreens x = 4.8, 
                                #for brown algae x = 2.7, for mixed algae x = 0.74 ugL; had to use var(log(abs(rnorm(1000,0.74,0.5)))) for mixed algae due to negative values
                                param_var = NULL,
                                EXO_var = c(1.5,7.7,6.3,19.5), #see 2B_DA_model_selection.R for calculation of these for each model
                                CTD_var = c(1.8,19.4,8.5,20.8),
                                EXO_CTD_var = c(0.3,4.3,7.2,18.9)){
  
  R = array(0, dim = c(n_states + n_params_obs, n_states + n_params_obs, n_step))
  
  state_var <- list()
  
    for(t in 1:n_step){
      if(all(is.na(obs[,,t]))){
        
        if(all(is.na(ancillary_obs[c(1:2),,t]))){
          state_var[[t]] = CTD_var
        } else if(all(is.na(ancillary_obs[c(3:5),,t]))){
          state_var[[t]] = EXO_var
        } else {
          state_var[[t]] = EXO_CTD_var
        }
        
      }
        else{
          state_var[[t]] = FP_var
        }
    }
    
  
  if(n_params_obs > 0){
    all_var = c(state_var, param_var)
  }else{
    all_var = state_var
  }
  
  
  for(i in 1:n_step){
    # variance is the same for each depth and time step; could make dynamic or varying by time step if we have good reason to do so
    R[,,i] = diag(all_var[[i]], n_states + n_params_obs, n_states + n_params_obs)
  }
  
  return(R)
}

R <- get_obs_error_matrix(n_states = 4, 
                          n_params_obs = 0, 
                          n_step = length(model_dates), 
                          FP_var = c(0.009,0.01,0.04,0.90), #remember log scale; estimated this using var(log(rnorm(1000,x,0.5))),
                          #where x is the mean of each PG from 2014-2020; for green algae x = 5.6, for bluegreens x = 4.8, 
                          #for brown algae x = 2.7, for mixed algae x = 0.74 ugL; had to use var(log(abs(rnorm(1000,0.74,0.5)))) for mixed algae due to negative values
                          param_var = NULL,
                          EXO_var = c(1.5,7.7,6.3,19.5), #see 2B_DA_model_selection.R for calculation of these for each model
                          CTD_var = c(1.8,19.4,8.5,20.8),
                          EXO_CTD_var = c(0.3,4.3,7.2,18.9))                                    

## ' Measurement operator matrix saying 1 if there is observation data available, 0 otherwise ----
#' for this project we need to tweak according to which data we are assimilating
#' for example, if we are only assimilating FP data, the function could stay as is
#' but if we want to assimilate both FP and EXO, then there should be a 1 if EITHER
#' of those is available, so this leads to the assimilation_expt parameter
#' and ancillary_obs inputs
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param obs observation matrix created with get_obs_matrix function
#' @param ancillary_obs ancillary observation matrix created with get_obs_matrix function
#' @param assimilation_expt choose from c("FP","FP_EXO","FP_CTD","all")
get_obs_id_matrix = function(n_states = 4, 
                             n_params_obs = 0, 
                             n_params_est = 2, 
                             n_step = length(model_dates), 
                             obs = obs,
                             ancillary_obs = ancillary_obs,
                             assimilation_expt = "all"){
  
  H = array(0, dim=c(n_states + n_params_obs, n_states + n_params_est, n_step))
  
  # order goes 1) states, 2)params for which we have obs, 3) params for which we're estimating but don't have obs
  
  if(assimilation_expt == "FP"){
  for(t in 1:n_step){
    H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
  }
  }
  
  if(assimilation_expt == "FP_EXO"){
    for(t in 1:n_step){
      if(all(is.na(ancillary_obs[c(1:2),,t]))){
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
      } else {
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(1, n_states + n_params_obs, n_states + n_params_obs)
      }
    }
  }
  
  if(assimilation_expt == "FP_CTD"){
    for(t in 1:n_step){
      if(all(is.na(ancillary_obs[c(3:5),,t]))){
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
      } else {
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(1, n_states + n_params_obs, n_states + n_params_obs)
      }
    }
  }
  
  if(assimilation_expt == "all"){
    for(t in 1:n_step){
      if(all(is.na(ancillary_obs[,,t]))){
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
      } else {
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(1, n_states + n_params_obs, n_states + n_params_obs)
    }
    }
  }
  
  return(H)
}

H <- get_obs_id_matrix(n_states = 4, 
                       n_params_obs = 0, 
                       n_params_est = 2, 
                       n_step = length(model_dates), 
                       obs = obs,
                       ancillary_obs = ancillary_obs,
                       assimilation_expt = "all")

##' vector for holding states and parameters for updating ----
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param n_en number of ensembles
get_Y_vector = function(n_states = 4, 
                        n_params_est = 2, 
                        n_step = length(model_dates), 
                        n_en = 30){
  
  Y = array(dim = c(n_states + n_params_est, n_step, n_en))
  
  return(Y)
}

Y <- get_Y_vector(n_states = 4, 
                  n_params_est = 2, 
                  n_step = length(model_dates), 
                  n_en = 30)

####' initialize Y vector with draws from distribution of obs ----
#' for our project, this will only be happening on 2021-01-01
#' I think it should be estimated the same way we estimate IC for the 
#' lagged FP observations in the driver data functions
#' params will be initialized using lambda.min and lambda.1se from glmnet output
#' @param Y Y vector
#' @param obs observation matrix
#' @param init_params initial values for estimated parameters (lambda.min and lambda.1se from glmnet output)
#' @param init_states values at which you'd like to initialize states
#' @param n_states_est number of states being estimated
#' @param n_params_est number of parameters being estimated
#' @param n_params_obs number of parameters being observed
#' @param n_step number of model timesteps
initialize_Y = function(Y = Y, 
                        init_params = c(0.35,1.07), 
                        init_states = c(6.15,2.26,1.93,3.60),
                        n_states_est = 4, 
                        n_params_est = 2, 
                        n_params_obs = 0, 
                        n_step = length(model_dates), 
                        n_en = 30, 
                        state_sd = rep(0.5,4), 
                        param_sd = c(0.05,0.05)){
  
  temp = array(dim = c(n_states_est + n_params_est, 1, n_en))
  
  temp[c(1:n_states_est),,] = array(log(rnorm(n = n_en * (n_states_est),
                                   mean = c(init_states),
                                   sd = c(state_sd))),
                         dim = c(c(n_states_est), n_en))
  temp[c((n_states_est+1):(n_states_est + n_params_est)),,] = array(rnorm(n = n_en * (n_params_est),
                                   mean = c(init_params),
                                   sd = c(param_sd)),
                         dim = c(c(n_params_est), n_en))
  
  
  Y[ , 1, ] = temp
  
  return(Y)
}

Y <- initialize_Y(Y = Y, 
                  init_params = c(0.35,1.07), 
                  init_states = c(6.15,2.26,1.93,3.60),
                  n_states_est = 4, 
                  n_params_est = 2, 
                  n_params_obs = 0, 
                  n_step = length(model_dates), 
                  n_en = 30, 
                  state_sd = rep(0.5,4), 
                  param_sd = c(0.05,0.05))

##' Kalman filter function ----
##' @param Y vector for holding states and parameters you're estimating
##' @param R observation error matrix
##' @param obs observations at current timestep
##' @param H observation identity matrix
##' @param n_en number of ensembles
##' @param cur_step current model timestep
kalman_filter = function(Y, R, obs, H, n_en, cur_step){
  
  cur_obs = obs[ , , cur_step]
  
  cur_obs = ifelse(is.na(cur_obs), 0, cur_obs) # setting NA's to zero so there is no 'error' when compared to estimated states
  
  ###### estimate the spread of your ensembles #####
  Y_mean = matrix(apply(Y[ , cur_step, ], MARGIN = 1, FUN = mean), nrow = length(Y[ , 1, 1])) # calculating the mean of each temp and parameter estimate
  delta_Y = Y[ , cur_step, ] - matrix(rep(Y_mean, n_en), nrow = length(Y[ , 1, 1])) # difference in ensemble state/parameter and mean of all ensemble states/parameters
  
  ###### estimate Kalman gain #########
  K = ((1 / (n_en - 1)) * delta_Y %*% t(delta_Y) %*% t(H[, , cur_step])) %*%
    qr.solve(((1 / (n_en - 1)) * H[, , cur_step] %*% delta_Y %*% t(delta_Y) %*% t(H[, , cur_step]) + R[, , cur_step]))
  
  ###### update Y vector ######
  for(q in 1:n_en){
    Y[, cur_step, q] = Y[, cur_step, q] + K %*% (cur_obs - H[, , cur_step] %*% Y[, cur_step, q]) # adjusting each ensemble using kalman gain and observations
  }
  return(Y)
}


#### Function to create driver data dataframe ----

#' Need to:
#' 1. get 30 ensemble members of FLARE
#' 2. pull in NOAA forecasts
#' 3. pull in inf forecasts
#' 4. pull in previous day's FP forecast for lag term or specify IC for lag term
#' 5. combine all the above in a data frame with the following columns:
#' date, column for each driver, ensemble member
#' 
#' @param issue_date issue date of forecast
#' @param start_of_run is this the first forecast you are running? T/F to tell 
#' the function how to specify the lagged term in the GLM
#' @param ic_lag vector of four values corresponding to log concentrations of
#' green algae, cyanobacteria, brown algae, and mixed algae in that order that you will use
#' to initialize the lag term if start_of_run == TRUE
#' @param sd_obs standard deviation of FP observations
#' @param fp_prev previous day's forecast of FP groups as matrix? dataframe?
#' with 
#' @param n_en number of ensemble members (will only work up to 30 due to NOAA)
#' @param forecast_horizon how many days into the future do you want to forecast? (max of 16 due to FLARE)

collate_driver_forecasts <- function(issue_date = model_dates[i],
                                     start_of_run = TRUE,
                                     ic_lag = c(6.15,2.26,1.93,3.60),
                                     sd_obs = 0.5,
                                     n_en = 30,
                                     forecast_horizon = 10){


#assign days for which you are making a forecast
fc_dates <- get_model_dates(model_start = issue_date, model_stop = (as.Date(issue_date) + forecast_horizon), time_step = 'days') 

#assign rows to pull from for FLARE ensemble
erow <- sample.int(100, n_en)

#'read in NOAA fc for start date
#'air_temperature = deg C
#'precipitation_flux = avg mm/sec I think?? convert to daily sum
#'relative_humidity = 0-1, convert to percent
#'surface_downwelling_longwave_flux_in_air = W/m2
#'surface_downwelling_shortwave_flux_in_air = W/m2
#'wind_speed = m/s

noaa <- read_csv(file.path(paste0("./00_Data_files/formatted_noaa/met_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates) %>%
  select(-issue_date,-air_pressure,-cloud_area_fraction,-specific_humidity) %>%
  mutate(precipitation_flux = precipitation_flux*60*60*24,#to convert to mm/day
         relative_humidity = relative_humidity*100) #to convert to percent

#read in inf fc for start date
inf <- read_csv(file.path(paste0("./00_Data_files/discharge_forecasts/inflow_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates) %>%
  select(-issue_date)

#read in FLARE fc for start date
wt <- read_csv(file.path(paste0("./00_Data_files/FLARE_forecasts/wt_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates & ensemble_member %in% erow) %>%
  select(-issue_date) 
wt$ensemble_member <- rep(1:n_en,each = length(fc_dates))

#specify data frame for lagged FP observations depending on whether this is
#the first forecast run or not
if(start_of_run == TRUE){
  fp <- data.frame("GreenAlgae_ugL_lag" = rnorm(n_en,ic_lag[1],sd_obs),
                   "Bluegreens_ugL_lag" = rnorm(n_en,ic_lag[2],sd_obs),
                   "BrownAlgae_ugL_lag" = rnorm(n_en,ic_lag[3],sd_obs),
                   "MixedAlgae_ugL_lag" = rnorm(n_en,ic_lag[4],sd_obs),
                   "ensemble_member" = rep(c(1:n_en),times = length(fc_dates)),
                   "fc_date" = rep(fc_dates, each = n_en))
} else {
  fp <- data.frame("GreenAlgae_ugL_lag" = rep(NA, n_en),
                   "Bluegreens_ugL_lag" = rep(NA, n_en),
                   "BrownAlgae_ugL_lag" = rep(NA, n_en),
                   "MixedAlgae_ugL_lag" = rep(NA, n_en),
                   "ensemble_member" = rep(c(1:n_en),times = length(fc_dates)),
                   "fc_date" = rep(fc_dates, each = n_en))
}

drivers <- left_join(noaa, inf, by = c("fc_date","ensemble_member")) %>%
  left_join(wt, by = c("fc_date","ensemble_member")) %>%
  left_join(fp, by = c("fc_date","ensemble_member")) %>%
  rename(daily_AirTemp_C = air_temperature,
         daily_RH_percent = relative_humidity,
         daily_Rain_Total_mm = precipitation_flux,
         daily_WindSpeed_Average_m_s = wind_speed,
         daily_ShortwaveRadiationDown_Average_W_m2 = surface_downwelling_shortwave_flux_in_air,
         daily_InfraredRadiationDown_Average_W_m2 = surface_downwelling_longwave_flux_in_air,
         Flow_cms = FLOW,
         date = fc_date)
drivers <- drivers[,c(1,10,2,4,3,7,6,5,9,11:14,8)] #making sure driver cols are in same order as expected by 'predict' function based on model selection script

log_vars <- read_csv("./2_Model_selection/logged_vars.csv")

for(i in 1:ncol(drivers)){
  if(colnames(drivers)[i] %in% log_vars$log_vars){
    drivers[,i] <- log(drivers[,i]+0.0001)
  }
}

fc_conv <- drivers
return(fc_conv)
}

fc_conv <- collate_driver_forecasts(issue_date = model_dates[i],
                                    start_of_run = TRUE,
                                    ic_lag = c(6.15,2.26,1.93,3.60),
                                    sd_obs = 0.5,
                                    n_en = 30,
                                    forecast_horizon = 10)
#### Function to create driver data matrix ----

#' matrix for holding driver data
#'
#' @param fc_conv dataframe which holds all the driver data; output of collate_driver_forecasts function 
#' @param n_drivers number of model drivers 
#' @param driver_colnames column names of the drivers in the driver dataframe 
#' @param issue_date issue date of forecast
#' @param n_en number of ensemble members
#' @param forecast_horizon forecast horizon in days
get_drivers = function(fc_conv = fc_conv, 
                       n_drivers = 12, 
                       driver_colnames = colnames(fc_conv)[2:13], 
                       issue_date = model_dates[i], 
                       n_en = 30,
                       forecast_horizon = 10){
  
  #assign days for which you are making a forecast
  fc_dates <- get_model_dates(model_start = issue_date, model_stop = (as.Date(issue_date) + forecast_horizon), time_step = 'days') 
  
  drivers_out = array(NA, dim = c(length(fc_dates), n_drivers, n_en))
  
  for(i in 1:n_drivers){
    for(j in 1:length(fc_dates)){
      fc1 <- fc_conv %>% filter(date == fc_dates[j])
      drivers_out[j,i,] = unlist(fc1[,driver_colnames[i]])
      
    }
  }
  
  return(drivers_out) 
}

forecast_drivers <- get_drivers(fc_conv = fc_conv, 
                       n_drivers = 12, 
                       driver_colnames = colnames(fc_conv)[2:13], 
                       issue_date = model_dates[i], 
                       n_en = 30,
                       forecast_horizon = 10)

####Function to aggregate driver ensembles for all forecast issue dates----
#' @param model_dates vector of dates for which you'd like to produce forecasts
#' @param n_en number of ensemble members
#' @param forecast_horizon forecast horizon in days
#' @param ic_lag vector of four values corresponding to log concentrations of
#' green algae, cyanobacteria, brown algae, and mixed algae in that order that you will use
#' to initialize the lag term if start_of_run == TRUE
#' @param sd_obs standard deviation of FP observations
#' #' @param n_drivers number of model drivers 
#' @param driver_colnames column names of the drivers in the driver dataframe 
#' BEST ESTIMATE THIS WILL TAKE ~15 MINUTES TO RUN
get_driver_list <- function(model_dates = model_dates,
                            n_en = 30,
                            forecast_horizon = 10,
                            ic_lag = c(6.15,2.26,1.93,3.60),
                            sd_obs = 0.5,
                            n_drivers = 12,
                            driver_colnames = c("Temp_C",                                   
                                                "daily_AirTemp_C",                          
                                                "daily_RH_percent",                         
                                                "daily_Rain_Total_mm" ,                     
                                                "daily_WindSpeed_Average_m_s" ,             
                                                "daily_ShortwaveRadiationDown_Average_W_m2",
                                                "daily_InfraredRadiationDown_Average_W_m2", 
                                                "Flow_cms" ,                                
                                                "GreenAlgae_ugL_lag"  ,                     
                                                "Bluegreens_ugL_lag"  ,                     
                                                "BrownAlgae_ugL_lag" ,                      
                                                "MixedAlgae_ugL_lag")){
  
  driver_list <- list()
  
  for(i in 1:length(model_dates)){
    if(file.exists(file.path(paste0("./00_Data_files/FLARE_forecasts/wt_",model_dates[i],".csv")))){
      if(i == 1){
        fc_conv <- collate_driver_forecasts(issue_date = model_dates[i],
                                            start_of_run = TRUE,
                                            ic_lag = ic_lag,
                                            sd_obs = sd_obs,
                                            n_en = n_en,
                                            forecast_horizon = forecast_horizon)
      } else {
        fc_conv <- collate_driver_forecasts(issue_date = model_dates[i],
                                            start_of_run = FALSE,
                                            ic_lag = ic_lag,
                                            sd_obs = sd_obs,
                                            n_en = n_en,
                                            forecast_horizon = forecast_horizon)
      }
      
      forecast_drivers <- get_drivers(fc_conv = fc_conv, 
                                      n_drivers = n_drivers, 
                                      driver_colnames = driver_colnames, 
                                      issue_date = model_dates[i], 
                                      n_en = n_en,
                                      forecast_horizon = forecast_horizon)
      
      date_last_file_found <- model_dates[i]
      
      
    } else{ #so this is the case where a NOAA forecast is missing
      
      missing_data_forecast_horizon = 16
      
      fc_conv <- collate_driver_forecasts(issue_date = date_last_file_found,
                                          start_of_run = FALSE,
                                          ic_lag = ic_lag,
                                          sd_obs = sd_obs,
                                          n_en = n_en,
                                          forecast_horizon = missing_data_forecast_horizon)
      
      missing_data_fc_dates = get_model_dates(model_start = model_dates[i], model_stop = (as.Date(model_dates[i]) + forecast_horizon), time_step = 'days')
      
      fc_conv <- fc_conv %>%
        filter(date %in% missing_data_fc_dates)
      
      forecast_drivers <- get_drivers(fc_conv = fc_conv, 
                                      n_drivers = n_drivers, 
                                      driver_colnames = driver_colnames, 
                                      issue_date = model_dates[i], 
                                      n_en = n_en,
                                      forecast_horizon = forecast_horizon)
      
    }
    
    driver_list[[i]] <- forecast_drivers
    
  } #end of for-loop along model_dates
  
  return(driver_list)
  
} #end of function

driver_list <- get_driver_list(model_dates = model_dates,
                               n_en = 30,
                               forecast_horizon = 10,
                               ic_lag = c(6.15,2.26,1.93,3.60),
                               sd_obs = 0.5,
                               n_drivers = 12,
                               driver_colnames = c("Temp_C",                                   
                                                   "daily_AirTemp_C",                          
                                                   "daily_RH_percent",                         
                                                   "daily_Rain_Total_mm" ,                     
                                                   "daily_WindSpeed_Average_m_s" ,             
                                                   "daily_ShortwaveRadiationDown_Average_W_m2",
                                                   "daily_InfraredRadiationDown_Average_W_m2", 
                                                   "Flow_cms" ,                                
                                                   "GreenAlgae_ugL_lag"  ,                     
                                                   "Bluegreens_ugL_lag"  ,                     
                                                   "BrownAlgae_ugL_lag" ,                      
                                                   "MixedAlgae_ugL_lag"))

####Function to run glmnet model ----
#' This function needs to:
#' 1. read in model file
#' 2. use driver and various draws from a runif bounded by lambda.min and lambda.1se
#' to generate a forecast with a horizon of 10 days
#' 
#' Notes: lambda.min and lambda.1se will be tuned so need to remember to re-populate them
#' This will need to be looped b/c there are 30 ensemble members for each element of driver_list
#' and 'predict' can't handle that; will need to have a loop of n_en
#' @param lambda.min lambda.min from glmnet output
#' @param lambda.1se lambda.1se from glmnet output
#' @param model cv.glmnet object
#' @param drivers array of driver data; dim = c(forecast_horizon, n_drivers, n_en)
#' @param n_en number of ensemble members
#' 
model <- readRDS("./2_Model_selection/fit_1.RData")

predict_phytos <- function(lambda.min = 0.35,
                           lambda.1se = 1.07,
                           model = model,
                           drivers = driver_list[[i]],
                           n_en = 30,
                           forecast_horizon = 10,
                           n_states = 4){
  
  y <- array(dim = c((forecast_horizon+1),n_states,n_en))
  
  lambda <- runif(n_en,lambda.min,lambda.1se)
  
  for(n in 1:n_en){
    y[,,n] <- predict(model, s = lambda[n], drivers[,,n])
  }
  
  return(list(y = y, lambda = lambda)) #come back to this and make sure y is in the right format to be plugged into drivers
  
}

driver_list[[2]][,c(9:12),] <- y #do we want to do this AFTER we've updated with obs from today or no?
#need to think about indexing here more; we can update first lag with today's data
#to run tmrw's forecast but how to update for the +2,+3...etc. day forecasts?

##Function to populate obs with predictions from DA models----
#' This function needs to:
#' 1. Determine which obs are available for today if FP is not available
#' 2. Apply the appropriate model to generate predictions of today's FP groups
#' 3. Uncertainty doesn't need to be applied here b/c it is treated like
#' obs uncertainty and is included in the obs variance matrix
#' The format needs to be inputtable to the obs matrix
#' 
#' @param x matrix of ancillary observations available for today
#' @param EXO.model lm object using EXO to predict FP
#' @param CTD.model lm object using CTD to predict FP
#' @param EXO.CTD.model lm object using EXO and CTD to predict FP
#' 

EXO.model <- readRDS("./2_Model_selection/lm_EXO.RData")
CTD.model <- readRDS("./2_Model_selection/lm_EXO.RData")
EXO.CTD.model <- readRDS("./2_Model_selection/lm_EXO.RData")

predict_FP_obs <- function(predictors = ancillary_obs[,,t],
                           EXO.model = EXO.model,
                           CTD.model = CTD.model,
                           EXO.CTD.model = EXO.CTD.model){
  
  if(all(is.na(predictors[c(1:2)]))){
    
    x = data.frame(CTD_chla_ugL_above = unlist(predictors[3]), CTD_chla_ugL = unlist(predictors[4]), CTD_chla_ugL_below = unlist(predictors[5]))
    
    FP.pred <- predict(CTD.model,x)
    
  } else if(all(is.na(predictors[c(3:5)]))){
    
    x = data.frame(daily_EXOChla_ugL_1 = unlist(predictors[1]), daily_EXOBGAPC_ugL_1 = unlist(predictors[2]))
    
    FP.pred <- predict(EXO.model,x) 
    
  } else {
    
    x = data.frame(CTD_chla_ugL_above = unlist(predictors[3]), CTD_chla_ugL = unlist(predictors[4]), CTD_chla_ugL_below = unlist(predictors[5]), daily_EXOChla_ugL_1 = unlist(predictors[1]), daily_EXOBGAPC_ugL_1 = unlist(predictors[2]))
    
    FP.pred <- predict(EXO.CTD.model,x)
  }
  
  return(FP.pred)
  
}

obs[,,t] <- FP.pred #will populate obs with these predictions


##' wrapper for running EnKF ----
#' 
#' @param n_en number of model ensembles 
#' @param start start date of model run 
#' @param stop date of model run
#' @param time_step model time step, defaults to days 
#' @param obs_file observation file 
#' @param driver_file driver data file 
#' @param n_states_est number of states we're estimating 
#' @param n_params_est number of parameters we're estimating 
#' @param n_params_obs number of parameters for which we have observations 
#' @param decay_init initial decay rate of DOC 
#' @param obs_cv coefficient of variation of observations 
#' @param param_cv coefficient of variation of parameters 
#' @param driver_cv coefficient of variation of driver data for DOC Load, Discharge out, and Lake volume, respectively 
#' @param init_cond_cv initial condition CV (what we're )
EnKF = function(n_en = 30, 
                start = "2021-01-01",
                stop = "2021-01-17", 
                time_step = 'days', 
                obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                driver_file = 'A_EcoForecast/Data/lake_c_data.rds',
                n_states_est = 1, 
                n_params_est = 1,
                n_params_obs = 0, 
                decay_init = 0.005, 
                obs_cv = 0.1,
                param_cv = 0.25,
                driver_cv = c(0.2, 0.2, 0.2),
                init_cond_cv = 0.1){
  
  
  n_en = n_en
  start = as.Date(start)
  stop = as.Date(stop)
  time_step = 'days' 
  dates = get_model_dates(model_start = start, model_stop = stop, time_step = time_step)
  n_step = length(dates)
  
  # get observation matrix
  obs_df = readRDS(obs_file) %>% 
    select(datetime, doc_lake) 
  
  drivers_df = readRDS(driver_file) %>% 
    select(datetime, doc_load, water_out, lake_vol) 
  
  n_states_est = n_states_est # number of states we're estimating 
  
  n_params_est = n_params_est # number of parameters we're calibrating
  
  n_params_obs = n_params_obs # number of parameters for which we have observations
  
  decay_init = decay_init # Initial estimate of DOC decay rate day^-1 
  
  doc_init = obs_df$doc_lake[min(which(!is.na(obs_df$doc_lake)))]
  
  state_cv = obs_cv #coefficient of variation of DOC observations 
  state_sd = state_cv * doc_init
  init_cond_sd = init_cond_cv * doc_init
  
  param_cv = param_cv #coefficient of variation of DOC decay 
  param_sd = param_cv * decay_init
  
  # driver data coefficient of variation for DOC Load, Discharge out, and Lake volume, respectively 
  driver_cv = driver_cv 
  
  
  # setting up matrices
  # observations as matrix
  obs = get_obs_matrix(obs_df = obs_df,
                       model_dates = dates,
                       n_step = n_step,
                       n_states = n_states_est)
  
  # Y vector for storing state / param estimates and updates
  Y = get_Y_vector(n_states = n_states_est,
                   n_params_est = n_params_est,
                   n_step = n_step,
                   n_en = n_en)
  
  # observation error matrix
  R = get_obs_error_matrix(n_states = n_states_est,
                           n_params_obs = n_params_obs,
                           n_step = n_step,
                           state_sd = state_sd,
                           param_sd = param_sd)
  
  # observation identity matrix
  H = get_obs_id_matrix(n_states = n_states_est,
                        n_params_obs = n_params_obs,
                        n_params_est = n_params_est,
                        n_step = n_step,
                        obs = obs)
  
  # initialize Y vector
  Y = initialize_Y(Y = Y, obs = obs, init_params = decay_init, n_states_est = n_states_est,
                   n_params_est = n_params_est, n_params_obs = n_params_obs,
                   n_step = n_step, n_en = n_en, state_sd = init_cond_sd, param_sd = param_sd)
  
  # get driver data with uncertainty - dim = c(n_step, driver, n_en) 
  drivers = get_drivers(drivers_df = drivers_df, 
                        model_dates = dates,
                        n_drivers = 3, 
                        driver_colnames = c('doc_load', 'water_out', 'lake_vol'), 
                        driver_cv = driver_cv, 
                        n_step = n_step, 
                        n_en = n_en) 
  
  # start modeling
  # remember will need to build in functionality here to populate each forecast's driver data
  # with the previous day's forecast of FP groups
  for(t in 2:n_step){
    for(n in 1:n_en){
      
      # run model; 
      model_output = predict_lake_doc(doc_load = drivers[t-1, 1, n], 
                                      doc_lake = Y[1, t-1, n], 
                                      lake_vol = drivers[t-1, 3, n], 
                                      water_out = drivers[t-1, 2, n],
                                      decay = Y[2, t-1, n])
      #' so would need to have something here where the model_output is for many days ahead of time (forecast horizon)
      #' but only the prediction for tmrw is stored in Y b/c that's what'll be used if no data is assimilated
      #' the full model output should be written to file each time so have a forecast for each day
      
      Y[1 , t, n] = model_output$doc_predict # store in Y vector
      Y[2 , t, n] = model_output$decay
    }
    # check if there are any observations to assimilate 
    ##this is where you need to edit if you've got obs from 
    ##ancillary sources to run the GLM and get predictions of FP
    if(any(!is.na(obs[ , , t]))){
      Y = kalman_filter(Y = Y,
                        R = R,
                        obs = obs,
                        H = H,
                        n_en = n_en,
                        cur_step = t) # updating params / states if obs available
    }
  }
  out = list(Y = Y, dates = dates, drivers = drivers, R = R, obs = obs, state_sd = state_sd)
  
  return(out)
}


