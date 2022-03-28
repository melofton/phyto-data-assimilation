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
pacman::p_load(tidyverse, lubridate,PerformanceAnalytics, glmnet)

##1. SUBSET COLLATED DATA FOR MODELS####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(!year(Date) == 2021)
colnames(obs)

#check for skewness and log-transform if needed
for(i in 2:ncol(obs)){
  if(abs(skewness(obs[,i],na.rm = TRUE)) > abs(skewness(log(obs[,i]+0.0001),na.rm = TRUE))){
    obs[,i] <- log(obs[,i]+0.0001)
  }
}

#a. Input: met, water temp, inflow, CTD @ depth of EXO sonde
a <- obs[complete.cases(obs[ ,c(1,20:25,19,26,3,8,11,14,17)]), c(1,20:25,19,26,3,8,11,14,17)]

#b. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below
b <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,8,11,14,17)]), c(1,20:25,19,26,3,8,11,14,17)]

#c. Input: met, water temp, inflow, EXO chl-a
c <- obs[complete.cases(obs[ ,c(1,20:25,19,26,5,8,11,14,17)]), c(1,20:25,19,26,5,8,11,14,17)]

#d. Input: met, water temp, inflow, EXO chl-a and EXO phyco
d <- obs[complete.cases(obs[ ,c(1,20:25,19,26,5:6,8,11,14,17)]), c(1,20:25,19,26,5:6,8,11,14,17)]

#e. Input: met, water temp, inflow, FP data at 0.5 m above and below the EXO sonde
e <- obs[complete.cases(obs[ ,c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,7,9,10,12,13,15,16,18,8,11,14,17)]

#f. Input: met, water temp, inflow, CTD @ depth of EXO sonde, EXO chl-a
f <- obs[complete.cases(obs[ ,c(1,20:25,19,26,3,5,8,11,14,17)]), c(1,20:25,19,26,3,5,8,11,14,17)]

#g. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, EXO chl-a
g <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,5,8,11,14,17)]), c(1,20:25,19,26,2:4,5,8,11,14,17)]

#h. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, EXO phyco
h <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,6,8,11,14,17)]), c(1,20:25,19,26,2:4,6,8,11,14,17)]

#i. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, EXO chl-a and EXO phyco
ii <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,5:6,8,11,14,17)]), c(1,20:25,19,26,2:4,5:6,8,11,14,17)]

#j. Input: met, water temp, inflow, CTD @ depth of EXO sonde, FP data at 0.5 m above
#   and below the EXO sonde
j <- obs[complete.cases(obs[ ,c(1,20:25,19,26,3,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,3,7,9,10,12,13,15,16,18,8,11,14,17)]

#k. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, FP data at 0.5 m above and below the EXO sonde
k <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,2:4,7,9,10,12,13,15,16,18,8,11,14,17)]

#l. Input: met, water temp, inflow, EXO chl-a, FP data at 0.5 m above and below the 
#   EXO sonde
l <- obs[complete.cases(obs[ ,c(1,20:25,19,26,5,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,5,7,9,10,12,13,15,16,18,8,11,14,17)]

#m Input: met, water temp, inflow, EXO phyco, FP data at 0.5 m above and below the 
#   EXO sonde
m <- obs[complete.cases(obs[ ,c(1,20:25,19,26,6,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,6,7,9,10,12,13,15,16,18,8,11,14,17)]

#n. Input: met, water temp, inflow, EXO chl-a and EXO phyco, FP data at 0.5 m above 
#   and below the EXO sonde
n <- obs[complete.cases(obs[ ,c(1,20:25,19,26,5:6,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,5:6,7,9,10,12,13,15,16,18,8,11,14,17)]

#o. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, EXO chl-a and EXO phyco, FP data at 0.5 m above and below the EXO sonde
o <- obs[complete.cases(obs[ ,c(1,20:25,19,26,2:4,5:6,7,9,10,12,13,15,16,18,8,11,14,17)]), c(1,20:25,19,26,2:4,5:6,7,9,10,12,13,15,16,18,8,11,14,17)]

combns <- list(a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h, ii=ii, j=j, k=k, l=l, m=m, n=n, o=o)
combns_x <- list()
combns_y <- list()

for (i in 1:length(combns)){
  combns_x[[i]] <- combns[[i]][,2:(length(combns[[i]])-4)]
  combns_y[[i]] <- combns[[i]][,(length(combns[[i]])-3):length(combns[[i]])]
}

##2. FIT MODELS

#useful resources
#https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf


#a. first assess what level of alpha is best by tuning lambda and alpha 
#   simultaneously with elastic net

aggregate_tuning_grid <- NULL

for (i in 1:length(combns)){
  
  #subset x and y variables
  x <- as.matrix(combns_x[[i]])
  y <- as.matrix(combns_y[[i]])
  
  # maintain the same folds across all models
  fold_id <- sample(1:10, size = nrow(y), replace=TRUE)
  
  # search across a range of alphas
  tuning_grid <- tibble::tibble(
    alpha      = seq(0, 1, by = .1),
    mse_min    = NA,
    mse_1se    = NA,
    lambda_min = NA,
    lambda_1se = NA,
    scenario_number = NA
  )
  
  for(j in seq_along(tuning_grid$alpha)) {
    
    # fit CV model for each alpha value
    fit <- cv.glmnet(x, y, family = "mgaussian", alpha = tuning_grid$alpha[j], foldid = fold_id)

    # extract MSE and lambda values
    tuning_grid$mse_min[j]    <- fit$cvm[fit$lambda == fit$lambda.min]
    tuning_grid$mse_1se[j]    <- fit$cvm[fit$lambda == fit$lambda.1se]
    tuning_grid$lambda_min[j] <- fit$lambda.min
    tuning_grid$lambda_1se[j] <- fit$lambda.1se
    tuning_grid$scenario_number[j] <- i
  }
  
  aggregate_tuning_grid <- rbind(aggregate_tuning_grid, tuning_grid)
  
  print(tuning_grid %>%
    mutate(se = mse_1se - mse_min) %>%
    ggplot(aes(alpha, mse_min)) +
    geom_line(size = 2) +
    geom_ribbon(aes(ymax = mse_min + se, ymin = mse_min - se), alpha = .25) +
    ggtitle(paste("MSE ± one standard error",i,sep = " ")))
  
}
  
#b. based on results here, I think probably safe to go with lasso in all cases

for (i in 1:length(combns)){
  
  #subset x and y variables
  x <- as.matrix(combns_x[[i]])
  y <- as.matrix(combns_y[[i]])
    
  # fit model
  fit <- cv.glmnet(x, y, family = "mgaussian", alpha = 1)
  
  # The numbers at the top of the plot just refer to the number of variables in the model.
  # The first and second vertical dashed lines represent the λ value with the minimum MSE and the largest λ value within one standard error of the minimum MSE.
  plot(fit)
  
  fit1 <- glmnet(x, y, family = "mgaussian", alpha = 1)
  plot(fit1, xvar = "lambda")
  
  min(fit$cvm)       # minimum MSE
  fit$lambda.min     # lambda for this min MSE

  fit$cvm[fit$lambda == fit$lambda.1se]  # 1 st.error of min MSE
  fit$lambda.1se  # lambda for this MSE
  
  # NEED TO EDIT THIS B/C RESPONSE IS MULTIVARIATE
  # coef(fit, s = "lambda.1se") %>%
  #   tidy() %>%
  #   filter(row != "(Intercept)") %>%
  #   ggplot(aes(value, reorder(row, value), color = value > 0)) +
  #   geom_point(show.legend = FALSE) +
  #   ggtitle("Influential variables") +
  #   xlab("Coefficient") +
  #   ylab(NULL)

}

# NEXT STEPS (28MAR22):

#1. Revisit this tmrw and re-evaluate:
#   does this still seem reasonable after a night's sleep?

# Specifically:
# a. the different combinations of data we are using (a-o)
# b. fitting ecah of these combinations and continuing with them rather than using
# elastic net to fit a single model using all possible data
# c. decision to use lasso for all final model fitting based on elastic net results

#2. Pull info and write visualizations:

# a. log lambda vs. MSE
# b. coef values for each functional group at min and 1se lambda values

#3. Figure out how you would write these models in R Code OR how you would use the
#   built-in "predict" function in either Jake's or Quinn's EnKF templates



