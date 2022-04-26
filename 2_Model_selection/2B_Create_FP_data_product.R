#2C Create FP data product
#Author: Mary Lofton
#Date: 20APR22

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

##READ IN MODELS####
EXO1model <- readRDS("./2_Model_selection/lm_EXO_1.RData")
EXO2model <- readRDS("./2_Model_selection/lm_EXO_2.RData")
CTDmodel <- readRDS("./2_Model_selection/lm_CTD.RData")
EXOCTDmodel <- readRDS("./2_Model_selection/lm_EXO_CTD.RData")

##READ IN DATA####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(Date >= "2018-10-01" & Date <= "2020-12-02")
obs <- obs[,c(1,2:4,8:9,25,28,31,34)]
check <- subset(obs, is.na(obs$daily_EXOChla_ugL_1))

#log FP vars
for(i in 2:ncol(obs)){
  obs[,i] <- log(obs[,i]+0.0001)
}

##WRITE FP PREDICTION FUNCTION####
##Function to populate obs with predictions from DA models
#' This function needs to:
#' 1. Determine which obs are available for today if FP is not available
#' 2. Apply the appropriate model to generate predictions of today's FP groups
#' 3. Uncertainty doesn't need to be applied here b/c it is treated like
#' obs uncertainty and is included in the obs variance matrix
#' 
#' @param obs_df dataframe of phytoplankton observations
#' @param EXO1.model lm object using EXO chla to predict FP
#' @param EXO2.model lm object using EXO to predict FP
#' @param CTD.model lm object using CTD to predict FP
#' @param EXO.CTD.model lm object using EXO and CTD to predict FP
#' 
predict_FP_obs <- function(obs_df = obs,
                           EXO1.model = EXO1model,
                           EXO2.model = EXO2model,
                           CTD.model = CTDmodel,
                           EXO.CTD.model = EXOCTDmodel){
  for(i in 1:nrow(obs_df)){
    
    if(all(is.na(obs[i,c(7:10)]))){
  
  if(all(is.na(obs[i,c(5:6)]))){
    
    x = data.frame(CTDChla_ugL_above = unlist(obs[i,2]), CTDChla_ugL = unlist(obs[i,3]), CTDChla_ugL_below = unlist(obs[i,4]))
    
    FP.pred <- predict(CTD.model,x)
    
  } else if(all(is.na(obs[i,c(2:4)]))){
    
    if(is.na(obs[i,6])){
      x = data.frame(daily_EXOChla_ugL_1 = unlist(obs[i,5]))
      
      FP.pred <- predict(EXO1.model,x)
    } else {
      
      x = data.frame(daily_EXOChla_ugL_1 = unlist(obs[i,5]), daily_EXOBGAPC_ugL_1 = unlist(obs[i,6]))
      
      FP.pred <- predict(EXO2.model,x) 
    }
    
  } else {
    
    x = data.frame(CTDChla_ugL_above = unlist(obs[i,2]), CTDChla_ugL = unlist(obs[i,3]), CTDChla_ugL_below = unlist(obs[i,4]), daily_EXOChla_ugL_1 = unlist(obs[i,5]), daily_EXOBGAPC_ugL_1 = unlist(obs[i,6]))
    
    FP.pred <- predict(EXO.CTD.model,x)
  }
      
      obs_df[i,c(7:10)] <- FP.pred
      
    }
    
  }
  
  return(obs_df)
  
}

FP.data.product <- predict_FP_obs(obs_df = obs,
                                  EXO1.model = EXO1model,
                                  EXO2.model = EXO2model,
                                  CTD.model = CTDmodel,
                                  EXO.CTD.model = EXOCTDmodel)

hist(exp(FP.data.product$GreenAlgae_ugL))
hist(exp(FP.data.product$Bluegreens_ugL))
hist(exp(FP.data.product$BrownAlgae_ugL))
hist(exp(FP.data.product$MixedAlgae_ugL))

FP.data.product <- FP.data.product[,c(1,7:10)]

write.csv(FP.data.product, file = "./2_Model_selection/FPDataProduct.csv",row.names = FALSE)
