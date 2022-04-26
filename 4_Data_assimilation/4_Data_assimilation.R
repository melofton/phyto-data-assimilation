#4_Data_assimilation
#Author: Mary Lofton, adapted from code provided at the GLEON
#ecological forecasting workshop by Jake Zwart
#Date: 13APR22

##Tasks for this script: ####
#'1. Generate 10-day hindcasts from 2021-01-01 to 2021-12-31 at a daily timestep
#'for four PFGs in FCR. 
#'

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate, mvtnorm)
source("./0_Function_library/EnKF_functions.R")
scenarios <- c("FP","FP_EXO_1","FP_EXO_2","FP_CTD","all")

##RUN FOR loop of ENKF FUNCTION with DA scenarios####
for(i in 1:length(scenarios)){
out <- EnKF(model_start = "2021-04-27",
            model_stop = "2021-11-09",
            time_step = 'days',
            obs_file = "./1_Data_wrangling/collated_obs_data.csv",
            assimilation_expt = scenarios[i],
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
            hindcast_output_folder = paste0("./4_Data_assimilation/hindcasts/",scenarios[i],"/"),
            resid_matrix = resid_matrix)

##SAVE ENKF OUTPUT####
#Note that hindcasts are saved within the function to ./4_Data_assimilation/hindcasts

saveRDS(out, file = file.path(paste0("./4_Data_assimilation/EnKF_output_",scenarios[i],".RData")))

}
