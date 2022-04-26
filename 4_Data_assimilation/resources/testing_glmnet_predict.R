##Testing code
##30MAR22
##Author: Mary Lofton

##Let's see how the predict function from glmnet works

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate,PerformanceAnalytics, glmnet, broom)

##1. SUBSET COLLATED DATA FOR MODELS####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(Date >= "2018-10-01" & Date <= "2020-12-02")
obs <- obs[,c(1,48:55)]

#check for skewness and log-transform if needed
for(i in 2:ncol(obs)){
  hist(unlist(obs[,i]),main = colnames(obs)[i])
  hist(unlist(log(obs[,i])),main = colnames(obs)[i])
}

#based on this Flow_cms and daily_Rain_Total_mm should be logged
obs$daily_Rain_Total_mm <- log(obs$daily_Rain_Total_mm+0.0001)
obs$Flow_cms <- log(obs$Flow_cms+0.0001)

#check for skewness and log-transform if needed
# log_vars <- NULL
# for(i in 2:ncol(obs)){
#   if(abs(skewness(obs[,i],na.rm = TRUE)) > abs(skewness(log(obs[,i]+0.0001),na.rm = TRUE))){
#     obs[,i] <- log(obs[,i]+0.0001)
#     print(colnames(obs)[i])
#     log_vars <- c(log_vars,colnames(obs)[i])
#   }
# }
# log_vars <- data.frame(log_vars)
# write.csv(log_vars,"./2_Model_selection/logged_vars.csv",row.names = FALSE)

#read in FP data product
FP.data.product <- read_csv("./2_Model_selection/FPDataProduct.csv") %>%
  mutate(GreenAlgae_ugL_lag = lag(GreenAlgae_ugL,1),
         Bluegreens_ugL_lag = lag(Bluegreens_ugL,1),
         BrownAlgae_ugL_lag = lag(BrownAlgae_ugL,1),
         MixedAlgae_ugL_lag = lag(MixedAlgae_ugL,1))
FP.data.product <- FP.data.product[,c(1,6:9,2:5)]

obs <- left_join(obs,FP.data.product,by = "Date")

obs_test <- obs %>% filter(year(Date) %in% c(2020))

#a. Input: met, water temp, inflow, lagged FP @ 1.6 m
a <- obs[complete.cases(obs[ , ]), ]
a_test <- obs_test[complete.cases(obs_test[ , ]), ]


combns <- list(a=a)
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

combns_test <- list(a_test=a_test)
combns_test_x <- list()
combns_test_y <- list()
my.rows <- sample.int(303,size = 10)
for (i in 1:length(combns_test)){
  combns_test_x[[i]] <- combns_test[[i]][c(my.rows),2:(length(combns_test[[i]])-4)]
  combns_test_y[[i]] <- combns_test[[i]][c(my.rows),(length(combns_test[[i]])-3):length(combns_test[[i]])]
}

test_x <- as.matrix(combns_test_x[[1]])
test_y <- as.matrix(combns_test_y[[1]])

lambda.test <- runif(3,fit$lambda.min,fit$lambda.1se)
hist(lambda.test)
pred <- predict(fit, s = 1.07, test_x)

mean((test_y - pred[,,1])^2)
pred
as.matrix(pred)

resid_matrix <- test_y - pred[,,1] ## this is the matrix that you would use as basis for process
#error covariance matrix; it represents the residuals from four days of predictions (rows)
#for each of the four PFGs (columns) using the fitted model with lambda.1se
#Predictions are for the last four observed data points in 2020.
write.csv(resid_matrix,"./4_Data_assimilation/resid_matrix.csv",row.names = FALSE)
