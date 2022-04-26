#Source code for functions to run EnKF
#Author: adapted by Mary Lofton using code from Jake Zwart provided to 2019 GLEON ecological forecasting workshop
#Date: 15APR22

rm(list = ls())
pacman::p_load(tidyverse, lubridate, mvtnorm)

##' retreive the model time steps based on start and stop dates and time step ----
#'
#' @param model_start model start date in date class
#' @param model_stop model stop date in date class
#' @param time_step model time step, defaults to daily timestep
get_model_dates = function(model_start, model_stop, time_step = 'days'){
  
  model_dates = seq.Date(from = as.Date(model_start), to = as.Date(model_stop), by = time_step)
  
  return(model_dates)
}

# model_dates <- get_model_dates(model_start = "2021-04-27",model_stop = "2021-11-09",time_step = 'days')

##' get observational data ----
#' function to read in data 
#'
#' @param obs_file filepath and filename containing observational data
#' @param assimilation_expt choose from c("FP","FP_EXO_1","FP_EXO_2","FP_CTD","all");
#' determines which data streams to include
#' @param missing_data_expt choose from TRUE or FALSE; determines whether or
#' not to drop every other FP observation

get_obs_df <- function(obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                       assimilation_expt = "all",
                       missing_data_expt = FALSE){
  
  obs <- read_csv(obs_file) 
  
  obs <- obs %>% 
    filter(year(Date) %in% c(2021)) %>%
    rename(datetime = Date)
  obs <- obs[,c(1,25,28,31,34,8,9,2:4)]
  
  for(i in 2:ncol(obs)){
    obs[,i] <- log(obs[,i]+0.0001)
  }
  
  if(assimilation_expt == "FP"){
    obs[,c("daily_EXOChla_ugL_1", "daily_EXOBGAPC_ugL_1", "CTDChla_ugL_above", "CTDChla_ugL", "CTDChla_ugL_below")] <- NA
  }
  
  if(assimilation_expt == "FP_EXO_1"){
    obs[,c("CTDChla_ugL_above", "CTDChla_ugL", "CTDChla_ugL_below", "daily_EXOBGAPC_ugL_1")] <- NA
  }
  
  if(assimilation_expt == "FP_EXO_2"){
    obs[,c("CTDChla_ugL_above", "CTDChla_ugL", "CTDChla_ugL_below")] <- NA
  }
  
  if(assimilation_expt == "FP_CTD"){
    obs[,c("daily_EXOChla_ugL_1", "daily_EXOBGAPC_ugL_1")] <- NA
  }
  
  if(missing_data_expt == TRUE){
    #place holder for later but honestly I don't think I'm going to get to this by JASM
  }
  
  return(obs)
  
}

# obs_df <- get_obs_df(obs_file = "./1_Data_wrangling/collated_obs_data.csv",
#                      assimilation_expt = "all",
#                      missing_data_expt = FALSE)


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

# obs <- get_obs_matrix(obs_df = obs_df, 
#                       model_dates = model_dates, 
#                       n_step = length(model_dates), 
#                       n_states = 4, 
#                       states = c("GreenAlgae_ugL",
#                                  "Bluegreens_ugL",
#                                  "BrownAlgae_ugL",
#                                  "MixedAlgae_ugL"))
# 
# ancillary_obs <- get_obs_matrix(obs_df = obs_df, 
#                       model_dates = model_dates, 
#                       n_step = length(model_dates), 
#                       n_states = 5, 
#                       states = c("daily_EXOChla_ugL_1","daily_EXOBGAPC_ugL_1","CTDChla_ugL_above","CTDChla_ugL","CTDChla_ugL_below"))

##' observation error matrix, should be a square matrix where ----
#'   col & row = the number of states and params for which you have observations
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_params_obs number of parameters for which we have observations
#' @param n_step number of model timesteps
#' @param FP_var vector of state observation standard deviation; assuming sd is constant through time; in our case sd of FP observations on log scale
#' @param param_var vector of parameter observation standard deviation; assuming sd is constant through time
#' @param EXO_1_var vector of standard deviations of EXO model predictions for each state
#' @param EXO_2_var vector of standard deviations of EXO model predictions for each state
#' @param CTD_var vector of standard deviations of CTD model predictions for each state
#' @param EXO_CTD_var vector of standard deviations of EXO + CTD model predictions for each state
#' @param obs obs matrix
#' @param ancillary_obs ancillary obs matrix

get_obs_error_matrix = function(n_states = 4, 
                                n_params_obs = 0, 
                                n_step = length(model_dates), 
                                FP_var = c(0.007,0.012,0.02,0.93), #estimated assuming obs error of 0.5 ug/L, e.g., var(log(rnorm(1000,x,0.5))) where x is mean of biomass for a PG
                                param_var = NULL,
                                EXO_1_var = c(1.5,14.5,6.7,20.8), #see 2B_DA_model_selection.R for calculation of these for each model
                                EXO_2_var = c(1.5,7.7,6.3,19.5),
                                CTD_var = c(1.8,19.3,8.7,20.7),
                                EXO_CTD_var = c(0.3,4.3,7.7,17.9),
                                obs = obs,
                                ancillary_obs = ancillary_obs){
  
  R = array(0, dim = c(n_states + n_params_obs, n_states + n_params_obs, n_step))
  
  state_var <- list()
  
    for(t in 1:n_step){
      if(all(is.na(obs[,,t]))){
        
        if(all(is.na(ancillary_obs[c(1:2),,t]))){
          state_var[[t]] = CTD_var
        } else if(all(is.na(ancillary_obs[c(3:5),,t]))){
          if(is.na(ancillary_obs[2,,t])){
            state_var[[t]] = EXO_1_var
          } else{
            state_var[[t]] = EXO_2_var
          }
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

# R <- get_obs_error_matrix(n_states = 4, 
#                           n_params_obs = 0, 
#                           n_step = length(model_dates), 
#                           FP_var = c(0.007,0.012,0.02,0.93), #estimated assuming obs error of 0.5 ug/L, e.g., var(log(rnorm(1000,x,0.5))) where x is mean of biomass for a PG
#                           param_var = NULL,
#                           EXO_1_var = c(1.5,14.5,6.7,20.8), #see 2B_DA_model_selection.R for calculation of these for each model
#                           EXO_2_var = c(1.5,7.7,6.3,19.5),
#                           CTD_var = c(1.8,19.3,8.7,20.7),
#                           EXO_CTD_var = c(0.3,4.3,7.7,17.9),
#                           obs = obs,
#                           ancillary_obs = ancillary_obs)

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
                             n_params_est = 0, 
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
  
  if(assimilation_expt == "FP_EXO_1"){
    for(t in 1:n_step){
      if(all(is.na(ancillary_obs[1,,t]))){
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
      } else {
        H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(1, n_states + n_params_obs, n_states + n_params_obs)
      }
    }
  }
  
  if(assimilation_expt == "FP_EXO_2"){
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

# H <- get_obs_id_matrix(n_states = 4, 
#                        n_params_obs = 0, 
#                        n_params_est = 0, 
#                        n_step = length(model_dates), 
#                        obs = obs,
#                        ancillary_obs = ancillary_obs,
#                        assimilation_expt = "all")

##' vector for holding states and parameters for updating ----
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param n_en number of ensembles
get_Y_vector = function(n_states = 4, 
                        n_params_est = 0, 
                        n_step = length(model_dates), 
                        n_en = 100){
  
  Y = array(dim = c(n_states + n_params_est, n_step, n_en))
  
  return(Y)
}

# Y <- get_Y_vector(n_states = 4, 
#                   n_params_est = 0, 
#                   n_step = length(model_dates), 
#                   n_en = 100)

####' initialize Y vector with draws from distribution of obs ----
#' for our project, this will only be happening on 2021-01-01
#' I think it should be estimated the same way we estimate IC for the 
#' lagged FP observations in the driver data functions
#' params will be initialized using lambda.min and lambda.1se from glmnet output
#' @param Y Y vector
#' @param obs observation matrix
#' @param init_params initial values for estimated parameters (lambda.min and lambda.1se from glmnet output)
#' @param init_states values at which you'd like to initialize states
#' @param n_states number of states being estimated
#' @param n_params_est number of parameters being estimated
#' @param n_params_obs number of parameters being observed
#' @param n_step number of model timesteps
# initialize_Y = function(Y = Y, 
#                         init_params = c(0.35,1.07), 
#                         init_states = c(6.15,2.26,1.93,3.60),
#                         n_states = 4, 
#                         n_params_est = 2, 
#                         n_params_obs = 0, 
#                         n_step = length(model_dates), 
#                         n_en = 30, 
#                         state_sd = rep(0.5,4), 
#                         param_sd = c(0.05,0.05)){
#   
#   temp = array(dim = c(n_states + n_params_est, 1, n_en))
#   
#   temp[c(1:n_states),,] = array(log(rnorm(n = n_en * (n_states),
#                                    mean = c(init_states),
#                                    sd = c(state_sd))),
#                          dim = c(c(n_states), n_en))
#   temp[c((n_states+1):(n_states + n_params_est)),,] = array(rnorm(n = n_en * (n_params_est),
#                                    mean = c(init_params),
#                                    sd = c(param_sd)),
#                          dim = c(c(n_params_est), n_en))
#   
#   
#   Y[ , 1, ] = temp
#   
#   return(Y)
# }

# Y <- initialize_Y(Y = Y, 
#                   init_params = c(0.35,1.07), 
#                   init_states = c(6.15,2.26,1.93,3.60),
#                   n_states = 4, 
#                   n_params_est = 2, 
#                   n_params_obs = 0, 
#                   n_step = length(model_dates), 
#                   n_en = 30, 
#                   state_sd = rep(0.5,4), 
#                   param_sd = c(0.05,0.05))

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
#' @param n_en number of ensemble members (will only work up to 30 due to NOAA)
#' @param forecast_horizon how many days into the future do you want to forecast? (max of 16 due to FLARE)

collate_driver_forecasts <- function(issue_date = model_dates[i],
                                     start_of_run = TRUE,
                                     ic_lag = c(2,0.0001,3,0.0001), #based on obs on 2021-04-26
                                     sd_obs = 0.5,
                                     n_en = 100,
                                     forecast_horizon = 10){


#assign days for which you are making a forecast
fc_dates <- get_model_dates(model_start = issue_date, model_stop = (as.Date(issue_date) + forecast_horizon), time_step = 'days') 

#assign rows to pull from for FLARE ensemble
if(n_en > 100){
  erow <- sample.int(100,n_en,replace = TRUE)
} else{
  erow <- sample.int(100, n_en)
}
if(n_en > 30){
  wrow <- sample.int(30,n_en,replace = TRUE)
} else{
  wrow <- sample.int(30, n_en)
}

#'read in NOAA fc for start date
#'air_temperature = deg C
#'precipitation_flux = avg mm/sec I think?? convert to daily sum
#'relative_humidity = 0-1, convert to percent
#'surface_downwelling_longwave_flux_in_air = W/m2
#'surface_downwelling_shortwave_flux_in_air = W/m2
#'wind_speed = m/s

noaa <- read_csv(file.path(paste0("./00_Data_files/formatted_noaa/met_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates & ensemble_member %in% wrow) %>%
  select(-issue_date,-air_pressure,-cloud_area_fraction,-specific_humidity) %>%
  mutate(precipitation_flux = log(precipitation_flux*60*60*24+0.0001),#to convert to mm/day and log-transform
         relative_humidity = relative_humidity*100) #to convert to percent
new_noaa <- data.frame(fc_date = rep(fc_dates,times = n_en),
                       ensemble_member = rep(wrow,each = length(fc_dates)),
                       new_ensemble_member = rep(1:n_en, each = length(fc_dates)))
final_noaa <- left_join(new_noaa,noaa,by = c("fc_date","ensemble_member")) %>%
  select(-ensemble_member) %>%
  rename(ensemble_member = new_ensemble_member)

#read in inf fc for start date
inf <- read_csv(file.path(paste0("./00_Data_files/discharge_forecasts/inflow_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates & ensemble_member %in% wrow) %>%
  select(-issue_date) 
new_inf <- data.frame(fc_date = rep(fc_dates,times = n_en),
                       ensemble_member = rep(wrow,each = length(fc_dates)),
                       new_ensemble_member = rep(1:n_en, each = length(fc_dates)))
final_inf <- left_join(new_inf,inf,by = c("fc_date","ensemble_member")) %>%
  select(-ensemble_member) %>%
  mutate(FLOW = log(FLOW + 0.0001)) %>%
  rename(ensemble_member = new_ensemble_member)

#read in FLARE fc for start date
wt <- read_csv(file.path(paste0("./00_Data_files/FLARE_forecasts/wt_",issue_date,".csv"))) %>%
  filter(fc_date %in% fc_dates & ensemble_member %in% erow) %>%
  select(-issue_date) 
new_wt <- data.frame(fc_date = rep(fc_dates,times = n_en),
                      ensemble_member = rep(erow,each = length(fc_dates)),
                      new_ensemble_member = rep(1:n_en, each = length(fc_dates)))
final_wt <- left_join(new_wt,wt,by = c("fc_date","ensemble_member")) %>%
  select(-ensemble_member) %>%
  rename(ensemble_member = new_ensemble_member)

#specify data frame for lagged FP observations depending on whether this is
#the first forecast run or not
if(start_of_run == TRUE){
  fp <- data.frame("GreenAlgae_ugL_lag" = abs(rnorm(n_en,ic_lag[1],sd_obs)),
                   "Bluegreens_ugL_lag" = abs(rnorm(n_en,ic_lag[2],sd_obs)),
                   "BrownAlgae_ugL_lag" = abs(rnorm(n_en,ic_lag[3],sd_obs)),
                   "MixedAlgae_ugL_lag" = abs(rnorm(n_en,ic_lag[4],sd_obs)),
                   "ensemble_member" = rep(c(1:n_en),times = length(fc_dates)),
                   "fc_date" = rep(fc_dates, each = n_en))
  for(c in 1:(ncol(fp)-2)){
    fp[,c] <- log(fp[,c]+0.0001)
  }
} else {
  fp <- data.frame("GreenAlgae_ugL_lag" = rep(NA, n_en),
                   "Bluegreens_ugL_lag" = rep(NA, n_en),
                   "BrownAlgae_ugL_lag" = rep(NA, n_en),
                   "MixedAlgae_ugL_lag" = rep(NA, n_en),
                   "ensemble_member" = rep(c(1:n_en),times = length(fc_dates)),
                   "fc_date" = rep(fc_dates, each = n_en))
}

drivers <- left_join(final_noaa, final_inf, by = c("fc_date","ensemble_member")) %>%
  left_join(final_wt, by = c("fc_date","ensemble_member")) %>%
  left_join(fp, by = c("fc_date","ensemble_member")) %>%
  rename(daily_AirTemp_C = air_temperature,
         daily_RH_percent = relative_humidity,
         daily_Rain_Total_mm = precipitation_flux,
         daily_WindSpeed_Average_m_s = wind_speed,
         daily_ShortwaveRadiationDown_Average_W_m2 = surface_downwelling_shortwave_flux_in_air,
         daily_InfraredRadiationDown_Average_W_m2 = surface_downwelling_longwave_flux_in_air,
         Flow_cms = FLOW,
         date = fc_date)
drivers <- drivers[,c(1,10,3,5,4,8,7,6,9,11:14,2)] #making sure driver cols are in same order as expected by 'predict' function based on model selection script

fc_conv <- drivers
return(fc_conv)

}

# fc_conv <- collate_driver_forecasts(issue_date = model_dates[i],
#                                     start_of_run = TRUE,
#                                     ic_lag = c(2,0.0001,3,0.0001), #based on obs on 2021-04-26
#                                     sd_obs = 0.5,
#                                     n_en = 100,
#                                     forecast_horizon = 10)
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
                       n_en = 100,
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

# forecast_drivers <- get_drivers(fc_conv = fc_conv,
#                        n_drivers = 12,
#                        driver_colnames = colnames(fc_conv)[2:13],
#                        issue_date = model_dates[i],
#                        n_en = 100,
#                        forecast_horizon = 10)

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
                            n_en = 100,
                            forecast_horizon = 10,
                            ic_lag = c(2,0.0001,3,0.0001),
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
      
      
    } else{ #so this is the case where a driver forecast is missing
      
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

# driver_list <- get_driver_list(model_dates = model_dates,
#                                n_en = 100,
#                                forecast_horizon = 10,
#                                ic_lag = c(2,0.0001,3,0.0001),
#                                sd_obs = 0.5,
#                                n_drivers = 12,
#                                driver_colnames = c("Temp_C",
#                                                    "daily_AirTemp_C",
#                                                    "daily_RH_percent",
#                                                    "daily_Rain_Total_mm" ,
#                                                    "daily_WindSpeed_Average_m_s" ,
#                                                    "daily_ShortwaveRadiationDown_Average_W_m2",
#                                                    "daily_InfraredRadiationDown_Average_W_m2",
#                                                    "Flow_cms" ,
#                                                    "GreenAlgae_ugL_lag"  ,
#                                                    "Bluegreens_ugL_lag"  ,
#                                                    "BrownAlgae_ugL_lag" ,
#                                                    "MixedAlgae_ugL_lag"))
# 
# saveRDS(driver_list, file = file.path(paste0("./4_Data_assimilation/driver_list.RData")))

####Function to run glmnet model ----
#' This function needs to:
#' 1. read in model file
#' 2. use drivers 
#' to generate a forecast with a horizon of 10 days w/ proc error
#' 
#' @param model cv.glmnet object
#' @param drivers array of driver data; dim = c(forecast_horizon, n_drivers, n_en)
#' @param n_en number of ensemble members
#' @param forecast_horizon forecast horizon in days
#' @param n_states = 4
#' @param resid_matrix matrix of residuals to calculate sigma for proc error

# #new version with mlm
# model <- readRDS("./2_Model_selection/forecast_fit.RData")
# resid_matrix <- read_csv("./4_Data_assimilation/resid_matrix.csv")


predict_phytos <- function(model = model,
                           drivers = driver_list[[i]],
                           n_en = 100,
                           forecast_horizon = 10,
                           n_states = 4,
                           resid_matrix = resid_matrix){
  
  y <- array(dim = c((forecast_horizon+1),n_states,n_en))
  
  sigma <- cov(resid_matrix)
  
  for(d in 1:(forecast_horizon+1)){
    for(n in 1:n_en){
    
    x = data.frame(Temp_C = unlist(drivers[d,1,n]),
                   daily_AirTemp_C = unlist(drivers[d,2,n]),
                   daily_RH_percent = unlist(drivers[d,3,n]),
                   daily_Rain_Total_mm = unlist(drivers[d,4,n]),
                   daily_WindSpeed_Average_m_s = unlist(drivers[d,5,n]),
                   daily_ShortwaveRadiationDown_Average_W_m2 = unlist(drivers[d,6,n]),
                   daily_InfraredRadiationDown_Average_W_m2 = unlist(drivers[d,7,n]),
                   Flow_cms = unlist(drivers[d,8,n]),
                   GreenAlgae_ugL_lag = unlist(drivers[d,9,n]),
                   Bluegreens_ugL_lag = unlist(drivers[d,10,n]),
                   BrownAlgae_ugL_lag = unlist(drivers[d,11,n]),
                   MixedAlgae_ugL_lag= unlist(drivers[d,12,n]))
    
    preds <- predict(model,x)
    y[d,,n] <- rmvnorm(1,matrix(preds),sigma)
    
    }
    
    if(d != 11){
      drivers[d+1,c(9:12),] <- y[d,,]
    }
    
  }
  
  return(list(y = y)) #come back to this and make sure y is in the right format to be plugged into drivers
  
}

# pred <- predict_phytos(model = model,
#                        drivers = driver_list[[i]],
#                        n_en = 100,
#                        forecast_horizon = 10,
#                        n_states = 4,
#                        resid_matrix = resid_matrix)

##Function to populate obs with predictions from DA models----
#' This function needs to:
#' 1. Determine which obs are available for today if FP is not available
#' 2. Apply the appropriate model to generate predictions of today's FP groups
#' 3. Uncertainty doesn't need to be applied here b/c it is treated like
#' obs uncertainty and is included in the obs variance matrix
#' The format needs to be inputtable to the obs matrix
#' 
#' @param predictors matrix of ancillary observations available for today
#' @param EXO1.model lm object using EXO chla to predict FP
#' @param EXO2.model lm object using EXO to predict FP
#' @param CTD.model lm object using CTD to predict FP
#' @param EXO.CTD.model lm object using EXO and CTD to predict FP
#' 


predict_FP_obs <- function(predictors = ancillary_obs[,,t],
                           EXO1.model = EXO1.model,
                           EXO2.model = EXO2.model,
                           CTD.model = CTD.model,
                           EXO.CTD.model = EXO.CTD.model){
  
  if(all(is.na(predictors[c(1:2)]))){
    
    x = data.frame(CTDChla_ugL_above = unlist(predictors[3]), CTDChla_ugL = unlist(predictors[4]), CTDChla_ugL_below = unlist(predictors[5]))
    
    FP.pred <- predict(CTD.model,x)
    
  } else if(all(is.na(predictors[c(3:5)]))){
    
    if(is.na(predictors[2])){
      x = data.frame(daily_EXOChla_ugL_1 = unlist(predictors[1]))
      
      FP.pred <- predict(EXO1.model,x)
    } else {
    
    x = data.frame(daily_EXOChla_ugL_1 = unlist(predictors[1]), daily_EXOBGAPC_ugL_1 = unlist(predictors[2]))
    
    FP.pred <- predict(EXO2.model,x) 
    }
    
  } else {
    
    x = data.frame(CTDChla_ugL_above = unlist(predictors[3]), CTDChla_ugL = unlist(predictors[4]), CTDChla_ugL_below = unlist(predictors[5]), daily_EXOChla_ugL_1 = unlist(predictors[1]), daily_EXOBGAPC_ugL_1 = unlist(predictors[2]))
    
    FP.pred <- predict(EXO.CTD.model,x)
  }
  
  return(FP.pred)
  
}

# EXO1.model <- readRDS("./2_Model_selection/lm_EXO_1.RData")
# EXO2.model <- readRDS("./2_Model_selection/lm_EXO_2.RData")
# CTD.model <- readRDS("./2_Model_selection/lm_CTD.RData")
# EXO.CTD.model <- readRDS("./2_Model_selection/lm_EXO_CTD.RData")
# 
# FP.pred <- predict_FP_obs(predictors = ancillary_obs[,,t],
#                           EXO1.model = EXO1.model,
#                           EXO2.model = EXO2.model,
#                           CTD.model = CTD.model,
#                           EXO.CTD.model = EXO.CTD.model)

##' wrapper for running EnKF ----
#' For this project, we want to generate a 10-day hindcast, so for example on Jan. 1,
#' we are generating a hindcast for Jan. 1-11, which gets us 10 days in the future.
#' 
#' For BOTH IC and the lag term on Jan. 1, we are using draws from a distribution
#' based around the last observation in Dec. 2020.
#' 
#' AFTER we generate the Jan. 1 forecast, we update the prediction for Jan. 1
#' with any available observations using the Kalman filter, and use that to
#' populate the lag for Jan. 2. 
#' 
#' To break that down in a little more detail:
#' 1. Generate hindcast (where you are working with t as today)
#' 2. Write hindcast to file
#' 3. Query for observations in obs and ancillary_obs
#'  a. if find FP obs, run Kalman filter
#'  b. if find no FP obs but DO find ancillary_obs
#'    i. determine which ancillary data model to apply (run predict_FP_obs function)
#'    ii. insert those predicted FP obs into obs matrix
#'    iii. run Kalman filter
#'  c. if find no FP obs but DO find ancillary_obs, do nothing
#' 4. Populate the drivers[[i+1]] with either forecast corrected with "obs" or just forecast
#' 
#' So now that I think through this...
#' I don't think we need to "initialize Y" b/c it'll get populated on day 1
#' given the way I'm gonna write the function
#' Actually, Y is not that helpful to us except for tracking the parameters 
#' anyways, because we'll be writing each forecast to file.
#' 
#' The test here will be 2 days in June which have NO data of any kind.
#' In this case, we would just write the prediction for "today" to tmrw's 
#' file w/ no updating. I think Y might be a good place to store predictions
#' post-update b/c then we can retrieve those easily?? not sure
#' 
#' @param model_start first forecast issue date
#' @param model_stop last forecast issue date
#' @param time_step time step of forecast; in theory could vary but honestly
#' will probably bonk if set to anything but 'days'
#' @param obs_file file of lake observations
#' @param assimilation_expt choose from c("FP","FP_EXO","FP_CTD","all");
#' determines which data streams to include
#' @param missing_data_expt aspirational; doesn't have any functionality yet;
#' should be set to FALSE
#' @param n_states_est number of FP states (4 for 4 PGs)
#' @param n_ancillary_data_streams number of ancillary data streams that are assimilated (5
#' for EXO and CTD above and below and at depth of EXO)
#' @param states_est character vector of column names for states
#' @param ancillary_data_streams character vector of column names for ancillary data streams
#' @param n_params_obs number of parameters for which we have observations
#' @param FP_var vector of variances for each FP group
#' @param param_var vector of variances for each estimated parameter
#' @param EXO_var vector of residual variances for EXO model predictions of FP
#' @param CTD_var vector of residual variances for CTD model predictions of FP
#' @param EXO_CTD_var vector of residual variances for EXO + CTD model predictions of FP
#' @param n_params_est number of parameters being estimated
#' @param n_en number of ensemble members
#' @param ic_lag vector of values with which to initialize lag term
#' @param sd_obs vector of standard deviations with which to initialize lag term
#' @param forecast_horizon forecast horizon in days
#' @param n_drivers number of drivers in forecast model
#' @param driver_colnames vector of column names for forecast model driver variables
#' @param init_lambda.min vector of values to initialize lambda.min, equal in length to n_en
#' @param init_lambda.1se vector of values to initialize lambda.1se, equal in length to n_en
#' @param model RData object containing forecast model, in this case a cv.glmnet object
#' @param EXO.model RData object containing DA model to convert EXO obs to FP, in this case an lm object
#' @param CTD.model RData object containing DA model to convert CTD obs to FP, in this case an lm object
#' @param EXO.CTD.model RData object containing DA model to convert EXO + CTD obs to FP, in this case an lm object
#' @param resid_matrix matrix to initialize sigma (process error covariance matrix)

resid_matrix <- read_csv("./4_Data_assimilation/resid_matrix.csv")
model <- readRDS("./2_Model_selection/forecast_fit.RData")
EXO1.model <- readRDS("./2_Model_selection/lm_EXO_1.RData")
EXO2.model <- readRDS("./2_Model_selection/lm_EXO_2.RData")
CTD.model <- readRDS("./2_Model_selection/lm_CTD.RData")
EXO.CTD.model <- readRDS("./2_Model_selection/lm_EXO_CTD.RData")

EnKF = function(model_start = "2021-04-27",
                model_stop = "2021-11-09",
                time_step = 'days',
                obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                assimilation_expt = "all",
                missing_data_expt = FALSE,
                n_states_est = 4,
                n_ancillary_data_streams = 5,
                states_est = c("GreenAlgae_ugL","Bluegreens_ugL","BrownAlgae_ugL","MixedAlgae_ugL"),
                ancillary_data_streams = c("daily_EXOChla_ugL_1","daily_EXOBGAPC_ugL_1","CTDChla_ugL_above","CTDChla_ugL","CTDChla_ugL_below"),
                n_params_obs = 0,
                FP_var = c(0.007,0.012,0.02,0.93), #estimated assuming obs error of 0.5 ug/L, e.g., var(log(rnorm(1000,x,0.5))) where x is mean of biomass for a PG
                param_var = NULL,
                EXO_1_var = c(1.5,14.5,6.7,20.8), #see 2B_DA_model_selection.R for calculation of these for each model
                EXO_2_var = c(1.5,7.7,6.3,19.5),
                CTD_var = c(1.8,19.3,8.7,20.7),
                EXO_CTD_var = c(0.3,4.3,7.7,17.9),
                n_params_est = 0,
                n_en = 100,
                ic_lag = c(2,0.0001,3,0.0001),
                sd_obs = 0.5,
                forecast_horizon = 10,
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
                                    "MixedAlgae_ugL_lag"),
                model = model,
                EXO1.model = EXO1.model,
                EXO2.model = EXO2.model,
                CTD.model = CTD.model,
                EXO.CTD.model = EXO.CTD.model,
                hindcast_output_folder = "./4_Data_assimilation/hindcasts/all/",
                resid_matrix = resid_matrix){
  
  #get vector of forecast issue dates
  model_dates <- get_model_dates(model_start = model_start,
                                 model_stop = model_stop,
                                 time_step = time_step)
  #set number of model timesteps
  n_step = length(model_dates)
  
  #get observation dataframe
  obs_df <- get_obs_df(obs_file = "./1_Data_wrangling/collated_obs_data.csv",
                       assimilation_expt = assimilation_expt,
                       missing_data_expt = missing_data_expt)
  
  #get observation matrix and ancillary observation matrix
  obs <- get_obs_matrix(obs_df = obs_df,
                        model_dates = model_dates,
                        n_step = n_step,
                        n_states = n_states_est,
                        states = states_est)
  
  ancillary_obs <- get_obs_matrix(obs_df = obs_df,
                                  model_dates = model_dates,
                                  n_step = n_step,
                                  n_states = n_ancillary_data_streams,
                                  states = ancillary_data_streams)
  
  #get observation error matrix
  R <- get_obs_error_matrix(n_states = n_states_est,
                            n_params_obs = n_params_obs,
                            n_step = n_step,
                            FP_var = FP_var,
                            param_var = param_var,
                            EXO_1_var = EXO_1_var,
                            EXO_2_var = EXO_2_var,
                            CTD_var = CTD_var,
                            EXO_CTD_var = EXO_CTD_var,
                            obs = obs,
                            ancillary_obs = ancillary_obs)
  
  #get observation id matrix
  H <- get_obs_id_matrix(n_states = n_states_est,
                         n_params_obs = n_params_obs,
                         n_params_est = n_params_est,
                         n_step = n_step,
                         obs = obs,
                         ancillary_obs = ancillary_obs,
                         assimilation_expt = assimilation_expt)
  
  #get Y vector
  Y <- get_Y_vector(n_states = n_states_est,
                    n_params_est = n_params_est,
                    n_step = n_step,
                    n_en = n_en)
  
  #read in drivers
  driver_list <- readRDS("./4_Data_assimilation/driver_list.RData")
  
  # start modeling
  for(t in 1:n_step){
    
    #' 1. Generate hindcast (where you are working with t as today)
      pred <- predict_phytos(model = model,
                             drivers = driver_list[[t]],
                             n_en = n_en,
                             forecast_horizon = forecast_horizon,
                             n_states = n_states_est,
                             resid_matrix = resid_matrix)
      
      #store states and parameters in Y 
      Y[1 , t,  ] = pred$y[1,1,] # store in Y vector
      Y[2 , t,  ] = pred$y[1,2,]
      Y[3 , t,  ] = pred$y[1,3,]
      Y[4 , t,  ] = pred$y[1,4,]
      
    
    #' 2. Write hindcast to file
    
    #assign issue_date
    issue_date = model_dates[t]
    #assign fc_dates
    fc_dates <- get_model_dates(model_start = issue_date, model_stop = (as.Date(issue_date) + forecast_horizon), time_step = time_step) 
    
    #format hindcast
    gr <- reshape2::melt(pred$y[,1,],varnames = c("fc_date","ensemble_member"), value.name = "GreenAlgae_ugL")
    bg <- reshape2::melt(pred$y[,2,],varnames = c("fc_date","ensemble_member"), value.name = "Bluegreens_ugL")
    br <- reshape2::melt(pred$y[,3,],varnames = c("fc_date","ensemble_member"), value.name = "BrownAlgae_ugL")
    mx <- reshape2::melt(pred$y[,4,],varnames = c("fc_date","ensemble_member"), value.name = "MixedAlgae_ugL")
    
    formatted_pred <- left_join(gr,bg, by = c("fc_date","ensemble_member")) %>%
      left_join(br, by = c("fc_date","ensemble_member")) %>%
      left_join(mx, by = c("fc_date","ensemble_member"))
    
    formatted_pred$issue_date <- issue_date
    formatted_pred$fc_date <- rep(fc_dates, times = n_en)
    
    #write to file
    write.csv(formatted_pred,file = file.path(paste0(hindcast_output_folder,"phytos_",model_dates[t],".csv")),row.names = FALSE)
    
    #' 3. Query for observations in obs and ancillary_obs
    #'  a. if find FP obs, run Kalman filter
    #'  b. if find no FP obs but DO find ancillary_obs
    #'    i. determine which ancillary data model to apply (run predict_FP_obs function)
    #'    ii. insert those predicted FP obs into obs matrix
    #'    iii. run Kalman filter
    #'  c. if find no FP obs or ancillary_obs, do nothing
    
    if(any(!is.na(obs[ , , t]))){ #if there are FP obs
      
      #then apply kalman filter
      Y = kalman_filter(Y = Y,
                        R = R,
                        obs = obs,
                        H = H,
                        n_en = n_en,
                        cur_step = t) # updating params / states if obs available
      
    } else if(any(!is.na(ancillary_obs[ , , t]))){ #otherwise if there are ancillary obs
      
      #get predictions for FP using ancillary data
      FP.pred <- predict_FP_obs(predictors = ancillary_obs[,,t],
                                EXO1.model = EXO1.model,
                                EXO2.model = EXO2.model,
                                CTD.model = CTD.model,
                                EXO.CTD.model = EXO.CTD.model)
      
      #put those predictions in the obs matrix
      obs[, , t] <- FP.pred
      
      #then apply kalman filter
      Y = kalman_filter(Y = Y,
                        R = R,
                        obs = obs,
                        H = H,
                        n_en = n_en,
                        cur_step = t) # updating params / states if obs available
      
    }
    
    #' 4. Populate the drivers[[i+1]] with either forecast corrected with "obs" or just forecast
    if(t != length(model_dates)){
      
      for(k in 1:(forecast_horizon+1)){ #I am certain there is a more efficient way to do this :-/
        
        driver_list[[t+1]][k ,c(9:12), ] <- Y[c(1:4), t, ]
        
      }
    
    }
    
  } #end of n_step loop
  
  #return output
  out = list(Y = Y, dates = model_dates, drivers = driver_list, R = R, obs = obs)
  
  return(out)
}


