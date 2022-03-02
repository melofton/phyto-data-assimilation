#2.2 Create phytoplankton data product
#Author: Mary Lofton
#Date: 22FEB22
#Adapted from scripts used in:


##################################SET-UP##############################################

#install to load and install other packages as needed
#install.packages('pacman')

#load packages
pacman::p_load(tidyverse, readxl, rjags, runjags, moments, coda)

#set a directory to use as a local file repository for plots if desire to write to file
my_directory <- "C:/Users/Mary Lofton/Dropbox/Ch5/Bayes_model_calibration_output"
write_plots <- FALSE

#assign model name
model_name <- c("phyto-model-data-product")

########################CALIBRATE MODELS##############################################

  #1) Grab data ---------------------------------------------------------

  y <- as.matrix(read_csv("./1_Data_wrangling/collated_phyto_data.csv"))
  y_exo_chla <- as.numeric(y[,3])
  y_ctd <- as.numeric(y[,2])
  y_fp <- as.numeric(y[,9])
  y_count <- as.numeric(y[,14])
  N = length(y[,1])

  #2) Assign specifications for model run -----------------------------------------------------
  
  data <- list(y_exo_chla=y_exo_chla, N=N)
  variable.names <- c("tau_obs")
  variable.namesout <- c("tau_obs","mu")
  init <- list(list(tau_obs = 0.1), list(tau_obs = 1), list(tau_obs = 5))
  params <- c("tau_obs")
  model=paste0("2.1_JAGS_models/",model_name, '.R') 
  
  #3) Run model (no edits, unless you want to change # of iterations) -------------------------------------------------------------
  j.model   <- jags.model (file = model,
                           data = data,
                           inits = init,
                           n.chains = 3)
  
  jags.out <- run.jags(model = model,
                       data = jags_plug_ins$data.model,
                       adapt =  5000,
                       burnin =  10000,
                       sample = 50000,
                       n.chains = 3,
                       inits=jags_plug_ins$init.model,
                       monitor = jags_plug_ins$variable.namesout.model)
  
  #convert to an MCMC list to calculate cross-correlation later
  jags.out.mcmc <- as.mcmc.list(jags.out)
  
  
  #6) Save output for calibration assessment
  
  #save predicted states
  Nmc = 10000
  out <- as.matrix(jags.out.mcmc)
  srow <- sample.int(nrow(out),Nmc,replace=TRUE)
  mus <- out[srow,grep("mu",colnames(out))]
  write.csv(mus,file = file.path("./5_Model_output/5.1_Calibration",paste0(model_name,'_predicted_states.csv')),row.names = FALSE)
  
  
  #plot parameters
  plot_parameters(params = jags_plug_ins$params.model,
                  write_plots = write_plots,
                  my_directory = my_directory)
  
  #calculate parameter summaries, effective sample size, and cross-correlations
  sum <- summary(jags.out, vars = jags_plug_ins$variable.names.model)
  crosscorr <- crosscorr(jags.out.mcmc[,c(jags_plug_ins$params.model)])
  
  #save results
  sink(file = file.path("./5_Model_output/5.1_Calibration",paste0(model_name,'_param_summary.txt')))
  print("Parameter summary")
  print(sum)
  print("Parameter cross-correlations")
  print(crosscorr)
  sink()
  
  
  #7) Save runjags output
  write.jagsfile(jags.out, file=file.path("./5_Model_output/5.1_Calibration",paste0(model_name,'_calibration.txt')),
                 remove.tags = TRUE, write.data = TRUE, write.inits = TRUE)
  
}

#Congratulations! You have run all the models. You may have a cookie.
