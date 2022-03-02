#Elastic net regularized model selection
#Author: Mary Lofton
#Date: 02MAR22

#Tasks to do in this script:

#1. Subset collated data for various test models; in all cases, target is FP data at
#   the depth of the EXO sonde in 2018-2020 (1 m before 2019-05-20, 1.6 m after);
#   for each subset, need to pare down to complete.cases only as glmnet cannot
#   handle NAs
#a. Input: met, water temp, inflow, CTD @ depth of EXO sonde
#b. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below
#c. Input: met, water temp, inflow, EXO chl-a
#d. Input: met, water temp, inflow, EXO chl-a and EXO phyco
#e. Input: met, water temp, inflow, FP data at 0.5 m above and below the EXO sonde
#f. Input: met, water temp, inflow, CTD @ depth of EXO sonde, 0.5 m above, and 0.5 m
#   below, EXO chl-a and EXO phyco, FP data at 0.5 m above and below the EXO sonde

#2. For each subset of data (a-f), check to see if each target vector and input
#   variable needs to be log-transformed (glmnet does standardization)

#3. Fit models to each subset
#Notes for Task 3:
#
#Use elastic net (0<alpha<1)
#Use cv.glmnet to make sure cross-validation is done automatically
#Use code from https://uc-r.github.io/regularized_regression#predict to tune
# both alpha and lambda parameters

#4. Plot relevant calibration/tuning/validation output
#a. Lambda vs. coefficient values
#b. Lambda vs. MSE
#c. Input vars vs. coefficient value
#d. Alpha vs. MSE

#5. Select best-fit model for test cases (a-f)