#Multivariate multiple regressions for data assimilation in JAGS
#Author: Mary Lofton
#Date: 18APR22

#Tasks to do in this script: ####

#1. Subset collated data for various test models
#   a. In all cases, target is FP data at
#   the depth of the EXO sonde in 2014-2020 (1 m before 2019-05-20, 1.6 m after);
#   
#   b. For each subset of data (EXO, CTD, or EXO + CTD), check to see if each target vector and input
#   variable needs to be log-transformed 

#2. Fit models to each subset
#Notes for Task 2:
#
#swapping to using Bayes b/c can extract
#PG-specific variances for use in uncertainty propagation as well as
#all the other uncertainty

#3. Plot relevant calibration/tuning/validation output
#summary
#variance for each group


##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate, mvtnorm, readxl, rjags, runjags, moments, coda)

##1. SUBSET COLLATED DATA FOR MODELS####

#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") 
obs <- obs[,c(1,2:4,8:9,25,28,31,34)]

#estimate prior for covariance matrix using 2014 data
get_cov <- obs %>%
  filter(year(Date) == 2014) %>%
  select(GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL) %>%
  filter(complete.cases(.))
#remember to log cause that's how model will be operating
for(i in 1:ncol(get_cov)){
  dat <- subset(get_cov[,i],get_cov[,i] > 0)
  my.min <- min(dat)
  get_cov[,i] <- log(get_cov[,i]+my.min)
}
cov_prior <- cov(get_cov)
cov_prior
#convert that to R; assuming k = 4 b/c want prior on correlation to be
#uniform distributions from (-1,1)
#sigma.guess <- sqrt(diag(R)/4)
#rearrange the above to solve for R
R = cov_prior^2*4
R
R = matrix(data = R, nrow = 3, ncol = 3)

# #see if we can return to where we started
# sigma.guess <- sqrt(diag(R)/4)
# sigma.guess
# #great ok so we can return what we started with - phew!
# 
# #how do the Wishart draws look?
# TAU <- rWishart(1, df = 4, Sigma = R)
# TAU <- matrix(data = TAU[,,1], nrow = 3, ncol = 3)
# VCOV <- solve(TAU)
# 
# pred <- rmvnorm(1, mean = unlist(get_cov[1,]), sigma = VCOV)
# ugL <- exp(pred)
# #not unreasonable; green algae and cyanos at 0, brown algae between 1-2

#check for skewness and log-transform if needed
for(i in 2:ncol(obs)){
  hist(unlist(obs[,i]),main = colnames(obs)[i])
  hist(unlist(log(obs[,i])),main = colnames(obs)[i])
}
#based on above think should log all
for(i in 2:ncol(obs)){
  dat <- subset(obs[,i],obs[,i] > 0)
  my.min <- min(dat)
  obs[,i] <- log(obs[,i]+my.min)
}

obs <- obs %>% filter(year(Date) %in% c(2015:2020))
colnames(obs)


#a. Input: EXO + CTD
a <- obs[complete.cases(obs[ , ]), ]

#b. Input: EXO sonde chl-a only
b <- obs[complete.cases(obs[ ,c(1,5,7:10)]), c(1,5,7:10)]

#c. Input: EXO sonde 
c <- obs[complete.cases(obs[ ,c(1,5:6,7:10)]), c(1,5:6,7:10)]

#d. Input: CTD
d <- obs[complete.cases(obs[ ,c(1,2:4,7:10)]), c(1,2:4,7:10)]

combns <- list(a=a, b=b, c=c, d=d)
combns_x <- list()
combns_y <- list()

for (i in 1:length(combns)){
  combns_x[[i]] <- combns[[i]][,2:(length(combns[[i]])-4)]
  combns_y[[i]] <- combns[[i]][,(length(combns[[i]])-3):length(combns[[i]])]
}

##2. FIT MODELS

#set a directory to use as a local file repository for plots if desire to write to file
my_directory <- "C:/Users/Mary Lofton/Documents/RProjects/phyto-data-assimilation/2_Model_selection/plots"
write_plots <- TRUE

##List of steps

####1. Calibrate model####

#a. Source helper functions ---------------------------------------------------------
source('2_Model_selection/JAGS_functions/model_calibration_plots.R')

#b. Model options => pick model -----------------------------------------------------

model_name = "CTD_EXO"
model=paste0("./2_Model_selection/JAGS_models/",model_name, '.R') #Do not edit


#c. Read in data for model ------------------------------------------------------------------------------------------------------------
fp = matrix(data = NA, nrow = nrow(combns_y[[1]]), ncol = 3)
fp[,1] <- unlist(combns_y[[1]][,1])
fp[,2] <- unlist(combns_y[[1]][,2])
fp[,3] <- unlist(combns_y[[1]][,3])

ctd_up <- unlist(combns_x[[1]][,1])
ctd <- unlist(combns_x[[1]][,2])
ctd_down <- unlist(combns_x[[1]][,3])

exo_chla <- unlist(combns_x[[1]][,4])
exo_phy <- unlist(combns_x[[1]][,5])

data = list(fp = fp,
            ctd_up = ctd_up,
            ctd = ctd,
            ctd_down = ctd_down,
            exo_chla = exo_chla,
            exo_phy = exo_phy,
            R = R,
            k = 4,
            beta.mg = as.vector(rep(0,times = 6)), 
            beta.vg = solve(diag(1E-03,6)),
            beta.mc = as.vector(rep(0,times = 6)), 
            beta.vc = solve(diag(1E-03,6)),
            beta.mb = as.vector(rep(0,times = 6)), 
            beta.vb = solve(diag(1E-03,6)),
            n = nrow(fp))

inits <- list()
inits[[1]] <- list(beta_g = c(rep(0.5,times = 6)),beta_c = c(rep(0.5,times = 6)),beta_b = c(rep(0.5,times = 6)))
inits[[2]] <- list(beta_g = c(rep(0,times = 6)),beta_c = c(rep(0,times = 6)),beta_b = c(rep(0,times = 6)))
inits[[3]] <- list(beta_g = c(rep(-0.5,times = 6)),beta_c = c(rep(-0.5,times = 6)),beta_b = c(rep(-0.5,times = 6)))

monitor = c("beta_g","beta_c","beta_b","TAU")


#e. Run model  -------------------------------------------------------------
j.model   <- jags.model (file = model,
                         data = data,
                         inits = inits,
                         n.chains = 3)

jags.out <- run.jags(model = model,
                     data = data,
                     adapt =  5000,
                     burnin =  10000,
                     sample = 50000,
                     n.chains = 3,
                     inits=inits,
                     monitor = monitor)

#convert to an MCMC list to calculate cross-correlation later
jags.out.mcmc <- as.mcmc.list(jags.out)

#f. Save output for calibration assessment -------------------------------

#plot parameters
params <- c("beta_g[1]","beta_g[2]","beta_g[3]","beta_g[4]","beta_g[5]","beta_g[6]","beta_c[1]","beta_c[2]","beta_c[3]","beta_c[4]","beta_c[5]","beta_c[6]","beta_b[1]","beta_b[2]","beta_b[3]","beta_b[4]","beta_b[5]","beta_b[6]","TAU[1,1]" , "TAU[2,1]" , "TAU[3,1]" , "TAU[1,2]" , "TAU[2,2]" , "TAU[3,2]", 
            "TAU[1,3]" , "TAU[2,3]" , "TAU[3,3]" )
plot_parameters(params = params,
                write_plots = write_plots,
                my_directory = my_directory)

#calculate parameter summaries, effective sample size, and cross-correlations
sum <- summary(jags.out, vars = params)
crosscorr <- crosscorr(jags.out.mcmc[,c(params)])

#save results
sink(file = file.path("./2_Model_selection/plots/",paste0(model_name,'_param_summary.txt')))
print("Parameter summary")
print(sum)
print("Parameter cross-correlations")
print(crosscorr)
sink()




