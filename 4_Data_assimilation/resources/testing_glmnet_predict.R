##Testing code
##30MAR22
##Author: Mary Lofton

##Let's see how the predict function from glmnet works

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate,PerformanceAnalytics, glmnet, broom)

##1. SUBSET COLLATED DATA FOR MODELS####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") 
colnames(obs)

#check for skewness and log-transform if needed

for(i in 2:ncol(obs)){
  if(abs(skewness(obs[,i],na.rm = TRUE)) > abs(skewness(log(obs[,i]+0.0001),na.rm = TRUE))){
    obs[,i] <- log(obs[,i]+0.0001)
  }
}

obs_test <- obs %>% filter(year(Date) %in% c(2021))
obs <- obs %>% filter(year(Date) %in% c(2018:2020))

b <- obs[complete.cases(obs[ ,c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]
b_test <- obs_test[complete.cases(obs_test[ ,c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]

combns <- list(b=b)
combns_x <- list()
combns_y <- list()

for (i in 1:length(combns)){
  combns_x[[i]] <- combns[[i]][,2:(length(combns[[i]])-4)]
  combns_y[[i]] <- combns[[i]][,(length(combns[[i]])-3):length(combns[[i]])]
}

for (i in 1:length(combns)){
  
  #subset x and y variables
  x <- as.matrix(combns_x[[i]])
  y <- as.matrix(combns_y[[i]])
  
  # fit model
  fit <- cv.glmnet(x, y, family = "mgaussian", alpha = 1)
  plot(fit)
  
}

combns_test <- list(b_test=b_test)
combns_test_x <- list()
combns_test_y <- list()

for (i in 1:length(combns_test)){
  combns_test_x[[i]] <- combns_test[[i]][c(1:4),2:(length(combns_test[[i]])-4)]
  combns_test_y[[i]] <- combns_test[[i]][c(1:4),(length(combns_test[[i]])-3):length(combns_test[[i]])]
}

test_x <- as.matrix(combns_test_x[[1]])
test_y <- as.matrix(combns_test_y[[1]])

lambda.test <- runif(3,fit$lambda.min,fit$lambda.1se)
hist(lambda.test)
pred <- predict(fit, s = fit$lambda.1se, test_x)

mean((test_y - pred[,,1])^2)
pred
as.matrix(pred)

test_y - pred[,,1] ## this is the matrix that you would use as basis for process
#error covariance matrix; it represents the residuals from four days of predictions (rows)
#for each of the four PFGs (columns) using the fitted model with lambda.1se
#Predictions are for the first four observed data points in 2021.
