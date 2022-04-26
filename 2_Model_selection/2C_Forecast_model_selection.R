#Elastic net regularized model selection
#Author: Mary Lofton
#Date: 02MAR22

#Tasks to do in this script: ####

#1. Subset collated data for various test models
#   a. In all cases, target is FP data at
#   the depth of the EXO sonde in 2018-2020 (1 m before 2019-05-20, 1.6 m after);
#   for each subset, need to pare down to complete.cases only as glmnet cannot
#   handle NAs - CHECK
#   b. For each subset of data (a-f), check to see if each target vector and input
#   variable needs to be log-transformed (glmnet does standardization) - CHECK

#2. Fit models to each subset
#Notes for Task 2:
#
#Use elastic net (0<alpha<1)
#Use cv.glmnet to make sure cross-validation is done automatically
#Use code from https://uc-r.github.io/regularized_regression#predict to tune
# both alpha and lambda parameters

#3. Plot relevant calibration/tuning/validation output
#a. Lambda vs. coefficient values
#b. Lambda vs. MSE
#c. Input vars vs. coefficient value
#d. Alpha vs. MSE

#4. Select best-fit model for test cases (a-f)

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

##1. SUBSET COLLATED DATA FOR MODELS####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(Date >= "2018-10-01" & Date <= "2020-11-02")
obs <- obs[,c(1,48:55)]

#check for skewness and log-transform if needed
for(i in 2:ncol(obs)){
  hist(unlist(obs[,i]),main = colnames(obs)[i])
  hist(unlist(log(obs[,i])),main = colnames(obs)[i])
}

#based on this Flow_cms and daily_Rain_Total_mm should be logged
obs$daily_Rain_Total_mm <- log(obs$daily_Rain_Total_mm+0.0001)
obs$Flow_cms <- log(obs$Flow_cms+0.0001)


#read in FP data product
FP.data.product <- read_csv("./2_Model_selection/FPDataProduct.csv") %>%
  mutate(GreenAlgae_ugL_lag = lag(GreenAlgae_ugL,1),
         Bluegreens_ugL_lag = lag(Bluegreens_ugL,1),
         BrownAlgae_ugL_lag = lag(BrownAlgae_ugL,1),
         MixedAlgae_ugL_lag = lag(MixedAlgae_ugL,1))
FP.data.product <- FP.data.product[,c(1,6:9,2:5)]

obs <- left_join(obs,FP.data.product,by = "Date")

#a. Input: met, water temp, inflow, lagged FP @ 1.6 m
a <- obs[complete.cases(obs[ , ]), ]

combns <- list(a=a)
combns_x <- list()
combns_y <- list()

for (i in 1:length(combns)){
  combns_x[[i]] <- combns[[i]][,2:(length(combns[[i]])-4)]
  combns_y[[i]] <- combns[[i]][,(length(combns[[i]])-3):length(combns[[i]])]
}

##2. FIT MODELS
PFGs <- c("GreenAlgae","Bluegreens","BrownAlgae","MixedAlgae")
scenarios <- c("","EXO_1","EXO_2","CTD","EXO_CTD")

sink(file = file.path("./2_Model_Selection/Model_selection_output.txt"))

i=1
  
  #subset x and y variables
  x <- as.matrix(combns_x[[i]])
  y <- as.matrix(combns_y[[i]])
  model.dat <- data.frame(cbind(x,y))
    
  # fit model
  fit <- lm(cbind(GreenAlgae_ugL,Bluegreens_ugL,BrownAlgae_ugL,MixedAlgae_ugL) ~ Temp_C + daily_AirTemp_C + daily_RH_percent + daily_Rain_Total_mm + daily_WindSpeed_Average_m_s + daily_ShortwaveRadiationDown_Average_W_m2 + daily_InfraredRadiationDown_Average_W_m2 + Flow_cms + GreenAlgae_ugL_lag + Bluegreens_ugL_lag + BrownAlgae_ugL_lag + MixedAlgae_ugL_lag, data = model.dat)
  saveRDS(fit, file = file.path(paste0("./2_Model_selection/forecast_fit.RData")))
  
  print(paste("Summary of",scenarios[i],sep = " "))
  print(summary(fit))
  print("Variance for each PG")
  print(sigma(fit)^2)
  


sink()


##get residuals matrix
#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(Date >= "2020-11-03" & Date <= "2020-12-02")
obs <- obs[,c(1,48:55)]

#based on this Flow_cms and daily_Rain_Total_mm should be logged
obs$daily_Rain_Total_mm <- log(obs$daily_Rain_Total_mm+0.0001)
obs$Flow_cms <- log(obs$Flow_cms+0.0001)


#read in FP data product
FP.data.product <- read_csv("./2_Model_selection/FPDataProduct.csv") %>%
  mutate(GreenAlgae_ugL_lag = lag(GreenAlgae_ugL,1),
         Bluegreens_ugL_lag = lag(Bluegreens_ugL,1),
         BrownAlgae_ugL_lag = lag(BrownAlgae_ugL,1),
         MixedAlgae_ugL_lag = lag(MixedAlgae_ugL,1))
FP.data.product <- FP.data.product[,c(1,6:9,2:5)]

obs <- left_join(obs,FP.data.product,by = "Date")

#a. Input: met, water temp, inflow, lagged FP @ 1.6 m
a <- obs[complete.cases(obs[ , ]), ]

combns <- list(a=a)
combns_x <- list()
combns_y <- list()

for (i in 1:length(combns)){
  combns_x[[i]] <- combns[[i]][,2:(length(combns[[i]])-4)]
  combns_y[[i]] <- combns[[i]][,(length(combns[[i]])-3):length(combns[[i]])]
}

i=1

#subset x and y variables
x <- data.frame(as.matrix(combns_x[[i]]))
y <- as.matrix(combns_y[[i]])
model.dat <- data.frame(cbind(x,y))

pred <- predict(fit,x)
mean((y - pred)^2)
resid_matrix <- y - pred
write.csv(resid_matrix,"./4_Data_assimilation/resid_matrix.csv",row.names = FALSE)
