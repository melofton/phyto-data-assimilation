#5 Analysis
#Author: Mary Lofton
#Date: 20APR22

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate, cowplot)
scenarios <- c("FP","FP_EXO_1","FP_EXO_2","FP_CTD","all")


##' retreive the model time steps based on start and stop dates and time step ----
#'
#' @param model_start model start date in date class
#' @param model_stop model stop date in date class
#' @param time_step model time step, defaults to daily timestep
get_model_dates = function(model_start, model_stop, time_step = 'days'){
  
  model_dates = seq.Date(from = as.Date(model_start), to = as.Date(model_stop), by = time_step)
  
  return(model_dates)
}

#RMSE 
#m = model values, o = observed values
rmse = function(m, o){
  sqrt(mean((m - o)^2, na.rm = TRUE))
}

##Tasks for this script:
#'1. Calculate RMSE at each forecast horizon for total biomass 
#'averaged from May-Oct for each scenario, and then do the same
#'for each FP group
#' a. find dates in 2021 with FP obs
#' b. retrieve forecasts for the 10 days preceding those dates (date - 10)
#' c. filter each forecast file to get the forecast for the date with obs
#' d. calculate RMSE for total biomass and each PG
#' e. assign forecast horizon to each forecast RMSE
#' f. aggregate all horizons
#' g. aggregate all obs
#' h. repeat for each scenario
#' 
## READ IN OBS ----
#read in collated data frame
obs <- read_csv("./1_Data_wrangling/collated_obs_data.csv") %>%
  filter(Date >= "2021-04-27" & Date <= "2021-11-09") %>%
  select(Date, GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL, MixedAlgae_ugL) %>%
  filter(complete.cases(.)) %>%
  mutate(TotalBiomass_ugL = rowSums(across(where(is.numeric)))) %>%
  mutate(GreenAlgae_ugL = log(GreenAlgae_ugL + 0.0001),
         Bluegreens_ugL = log(Bluegreens_ugL + 0.0001),
         BrownAlgae_ugL = log(BrownAlgae_ugL + 0.0001),
         MixedAlgae_ugL = log(MixedAlgae_ugL + 0.0001),
         TotalBiomass_ugL = log(TotalBiomass_ugL + 0.0001))

forecast_horizon = 10

for(i in 1:length(scenarios)){
  
  hindcast_output_folder = paste0("./4_Data_assimilation/hindcasts/",scenarios[i],"/")
  
  for(j in 1:length(obs$Date)){
    
    fc_dates <- get_model_dates(model_start = (as.Date(unlist(obs[j,"Date"]), origin = "1970-01-01") - forecast_horizon), model_stop = as.Date(unlist(obs[j,"Date"]), origin = "1970-01-01"), time_step = 'days') 
    
    if(j == 1){
      fc_dates <- fc_dates[-1]
    }
    
    for(k in 1:length(fc_dates)){
      
      fc <- read_csv(file.path(paste0("./4_Data_assimilation/hindcasts/",scenarios[i],"/phytos_",fc_dates[k],".csv"))) %>%
        filter(fc_date == as.Date(unlist(obs[j,"Date"]), origin = "1970-01-01")) %>%
        select(-issue_date) %>%
        mutate(GreenAlgae_ugL_exp = exp(GreenAlgae_ugL),
               Bluegreens_ugL_exp = exp(Bluegreens_ugL),
               BrownAlgae_ugL_exp = exp(BrownAlgae_ugL),
               MixedAlgae_ugL_exp = exp(MixedAlgae_ugL)) %>%
        rowwise() %>% 
        mutate(TotalBiomass_ugL = sum(c_across(ends_with("exp")), na.rm = T)) %>%
        ungroup() %>%
        mutate(TotalBiomass_ugL = log(TotalBiomass_ugL + 0.0001)) %>%
        summarize(GreenAlgae_ugL_pred = median(GreenAlgae_ugL),
                  Bluegreens_ugL_pred = median(Bluegreens_ugL),
                  BrownAlgae_ugL_pred = median(BrownAlgae_ugL),
                  MixedAlgae_ugL_pred = median(MixedAlgae_ugL),
                  TotalBiomass_ugL_pred = median(TotalBiomass_ugL)) %>%
        add_column(fc_horizon = as.Date(fc_dates[k]) - as.Date(unlist(obs[j,"Date"]), origin = "1970-01-01"),
                   obs_date = as.Date(unlist(obs[j,"Date"]), origin = "1970-01-01"),
                   Scenario = scenarios[i],
                   fc_date = fc_dates[k],
                   GreenAlgae_ugL_obs = obs$GreenAlgae_ugL[j],
                   Bluegreens_ugL_obs = obs$Bluegreens_ugL[j],
                   BrownAlgae_ugL_obs = obs$BrownAlgae_ugL[j],
                   MixedAlgae_ugL_obs = obs$MixedAlgae_ugL[j],
                   TotalBiomass_ugL_obs = obs$TotalBiomass_ugL[j]) 

      if(i == 1 & j == 1 & k==1){
        appended_fc <- fc
      } else {
        appended_fc <- rbind(appended_fc,fc)
      }
      
      
    }
  }
}

appended_fc <- appended_fc %>%
  mutate(fc_horizon = as.numeric(fc_horizon))
write.csv(appended_fc,"./5_Analysis/pred_v_obs_log.csv",row.names = FALSE)

appended_fc_notlog <- read_csv("./5_Analysis/pred_v_obs.csv")

fc_horizons <- unique(appended_fc$fc_horizon)

for(s in 1:length(scenarios)){
  
  
  
  for(h in 1:length(fc_horizons)){
    
    dat <- appended_fc %>%
      filter(Scenario == scenarios[s] & fc_horizon == fc_horizons[h])
    
    assessment <- data.frame(rmse_ga = rmse(dat$GreenAlgae_ugL_pred,dat$GreenAlgae_ugL_obs),
                             rmse_bg = rmse(dat$Bluegreens_ugL_pred,dat$Bluegreens_ugL_obs),
                             rmse_ba = rmse(dat$BrownAlgae_ugL_pred,dat$BrownAlgae_ugL_obs),
                             rmse_ma = rmse(dat$MixedAlgae_ugL_pred,dat$MixedAlgae_ugL_obs),
                             rmse_tot = rmse(dat$TotalBiomass_ugL_pred,dat$TotalBiomass_ugL_obs),
                             fc_horizon = fc_horizons[h],
                             scenario = scenarios[s])
    
    if(h == 1 & s == 1){
      appended_assessment <- assessment
    } else {
      appended_assessment <- rbind(appended_assessment,assessment)
    }
    
  }
    
    
  
}

ggplot(dat = appended_assessment, aes(x = fc_horizon, y = rmse_tot, group = scenario,color = scenario))+
  geom_point()+
  theme_classic()+
  xlim(-10,0)

ggplot(dat = appended_assessment, aes(x = fc_horizon, y = rmse_ga, group = scenario,color = scenario))+
  geom_point()+
  theme_classic()+
  xlim(-10,0)

ggplot(dat = appended_assessment, aes(x = fc_horizon, y = rmse_bg, group = scenario,color = scenario))+
  geom_point()+
  theme_classic()+
  xlim(-10,0)

ggplot(dat = appended_assessment, aes(x = fc_horizon, y = rmse_ba, group = scenario,color = scenario))+
  geom_point()+
  theme_classic()+
  xlim(-10,0)

ggplot(dat = appended_assessment, aes(x = fc_horizon, y = rmse_ma, group = scenario,color = scenario))+
  geom_point()+
  theme_classic()+
  xlim(-10,0)

#consider plotting these as percent increase or decrease in RMSE compared to null (FP) - HA!

null_model <- appended_assessment %>%
  filter(scenario == "FP")

for(s in 2:length(scenarios)){
  
  scen <- appended_assessment %>%
    filter(scenario == scenarios[s]) 
  
  perc_change <- data.frame(perc_change_rmse_ga = ((scen$rmse_ga - null_model$rmse_ga)/null_model$rmse_ga)*100,
                            perc_change_rmse_bg = ((scen$rmse_bg - null_model$rmse_bg)/null_model$rmse_bg)*100,
                            perc_change_rmse_ba = ((scen$rmse_ba - null_model$rmse_ba)/null_model$rmse_ba)*100,
                            perc_change_rmse_ma = ((scen$rmse_ma - null_model$rmse_ma)/null_model$rmse_ma)*100,
                            perc_change_rmse_tot = ((scen$rmse_tot - null_model$rmse_tot)/null_model$rmse_tot)*100,
                            scenario = scenarios[s],
                            fc_horizon = scen$fc_horizon)
  
  if(s==2){
    appended_perc_change <- perc_change
  } else {
    appended_perc_change <- rbind(appended_perc_change,perc_change)
  }
  
  
  
}

tot <- ggplot(dat = appended_perc_change, aes(x = fc_horizon, y = perc_change_rmse_tot, group = scenario, color = scenario))+
  geom_point(size = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  scale_y_continuous(limits = c(-25,25), breaks = seq(from = -25, to = 25,by = 5))+
  xlab("Forecast horizon (days)")+
  ylab("RMSE percent change from null")+
  scale_x_continuous(limits = c(-10,0), breaks = c(-10:0))+
  labs(color='Assimilation experiment') +
  theme(legend.spacing.y = unit(1, 'cm'), panel.grid.minor = element_blank())+
  guides(color = guide_legend(byrow = TRUE))+
  ggtitle(expression(paste("Total Biomass ","(",mu,"g ",L^-1,")")))
tot
ggsave(tot, filename = "./5_Analysis/perc_change_rmse_tot.tif",height = 4, width = 6,
       units = "in", dpi = 300, dev = "tiff")

ga <- ggplot(dat = appended_perc_change, aes(x = fc_horizon, y = perc_change_rmse_ga, group = scenario, color = scenario))+
  geom_point(size = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  scale_y_continuous(limits = c(-25,25), breaks = seq(from = -25, to = 25,by = 5))+
  xlab("Forecast horizon (days)")+
  ylab("RMSE percent change from null")+
  scale_x_continuous(limits = c(-10,0), breaks = c(-10:0))+
  labs(color='Assimilation experiment') +
  theme(legend.spacing.y = unit(1, 'cm'), panel.grid.minor = element_blank())+
  guides(color = guide_legend(byrow = TRUE))+
  ggtitle(expression(paste("Green Algae ","(",mu,"g ",L^-1,")")))
ga
ggsave(ga, filename = "./5_Analysis/perc_change_rmse_ga.tif",height = 4, width = 6,
       units = "in", dpi = 300, dev = "tiff")

bg <- ggplot(dat = appended_perc_change, aes(x = fc_horizon, y = perc_change_rmse_bg, group = scenario, color = scenario))+
  geom_point(size = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  scale_y_continuous(limits = c(-25,25), breaks = seq(from = -25, to = 25,by = 5))+
  xlab("Forecast horizon (days)")+
  ylab("RMSE percent change from null")+
  scale_x_continuous(limits = c(-10,0), breaks = c(-10:0))+
  labs(color='Assimilation experiment') +
  theme(legend.spacing.y = unit(1, 'cm'), panel.grid.minor = element_blank())+
  guides(color = guide_legend(byrow = TRUE))+
  ggtitle(expression(paste("Cyanobacteria ","(",mu,"g ",L^-1,")")))
bg
ggsave(bg, filename = "./5_Analysis/perc_change_rmse_bg.tif",height = 4, width = 6,
       units = "in", dpi = 300, dev = "tiff")

ba <- ggplot(dat = appended_perc_change, aes(x = fc_horizon, y = perc_change_rmse_ba, group = scenario, color = scenario))+
  geom_point(size = 2)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  scale_y_continuous(limits = c(-25,25), breaks = seq(from = -25, to = 25,by = 5))+
  xlab("Forecast horizon (days)")+
  ylab("RMSE percent change from null")+
  scale_x_continuous(limits = c(-10,0), breaks = c(-10:0))+
  labs(color='Assimilation experiment') +
  theme(legend.spacing.y = unit(1, 'cm'), panel.grid.minor = element_blank())+
  guides(color = guide_legend(byrow = TRUE))+
  ggtitle(expression(paste("Brown Algae ","(",mu,"g ",L^-1,")")))
ba
ggsave(ba, filename = "./5_Analysis/perc_change_rmse_ba.tif",height = 4, width = 6,
       units = "in", dpi = 300, dev = "tiff")

#How often did each scenario correctly predict the dominant group?
groups <- appended_fc %>%
  select(GreenAlgae_ugL_obs, Bluegreens_ugL_obs, BrownAlgae_ugL_obs, MixedAlgae_ugL_obs)
dom_obs <- colnames(groups)[apply(groups,1,which.max)]
dom_obs <- dom_obs %>%
  strsplit( "_" ) %>%
  sapply( "[", 1 )

pred.groups <- appended_fc %>%
  select(GreenAlgae_ugL_pred, Bluegreens_ugL_pred, BrownAlgae_ugL_pred, MixedAlgae_ugL_pred)
dom_pred <- colnames(pred.groups)[apply(pred.groups,1,which.max)]
dom_pred <- dom_pred %>%
  strsplit( "_" ) %>%
  sapply( "[", 1 )

dom_pred_correct <- ifelse(dom_obs == dom_pred, "yes","no")

appended_fc <- appended_fc %>%
  add_column(dom_group_obs = dom_obs,
             dom_group_pred = dom_pred,
             dom_pred_correct = dom_pred_correct)

check <- appended_fc %>%
  filter(dom_group_obs == "MixedAlgae")

dom <- ggplot(data = appended_fc, aes(x = obs_date, y = dom_pred_correct, group = Scenario, color = Scenario))+
  geom_jitter(size = 2)+
  theme_classic()+
  xlab("")+
  ylab("Correct dominant group predicted?")+
  scale_x_date(date_breaks = "1 week", date_labels = "%m/%d")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(legend.spacing.y = unit(0.5, 'cm'))+
  guides(color = guide_legend(byrow = TRUE))+
  scale_color_discrete(type = c("slategray","slategray1","slategray2","slategray3","slategray4"))+
  labs(color='Assimilation experiment') 
dom
ggsave(dom, filename = "./5_Analysis/dom_group_predictions.tif",height = 4, width = 10,
       units = "in", dpi = 300, dev = "tiff")

fp <- appended_fc_notlog %>%
  select(obs_date, GreenAlgae_ugL_obs, Bluegreens_ugL_obs, BrownAlgae_ugL_obs) %>%
  gather(GreenAlgae_ugL_obs:BrownAlgae_ugL_obs, key = "PG",value = "ugL")

fp.plot <- ggplot(data = fp, aes(x = obs_date, y = ugL, group = PG, color = PG))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  theme_classic()+
  xlab("")+
  ylab(expression(paste("Biomass ","(",mu,"g ",L^-1,")")))+
  scale_x_date(date_breaks = "1 week", date_labels = "%m/%d")+
  scale_color_discrete(type = c("cyan4","brown","lightgreen"), labels = c("Cyanobacteria","Brown Algae","Green Algae"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(legend.spacing.y = unit(0.5, 'cm'),axis.text.x = element_blank())+
  guides(color = guide_legend(byrow = TRUE))+
  labs(color='Phytoplankton group') 
fp.plot

plot<-plot_grid(fp.plot,dom, align='v', scale = 1,
                nrow = 2, ncol = 1)

ggsave(plot, filename = "./5_Analysis/dom_group_predictions.tif",height = 5, width = 6.5,
       units = "in", dpi = 300, dev = "tiff")


#Example hindcast
ggplot(data = appended_fc, aes(x = fc_date, y = Bluegreens_ugL_pred))+
  geom_point()+
  theme_classic()

ggplot(data = appended_fc, aes(x = fc_date, y = Bluegreens_ugL_obs))+
  geom_point()+
  theme_classic()

which.min(abs(appended_fc$Bluegreens_ugL_pred - appended_fc$Bluegreens_ugL_obs))
check <- appended_fc %>%
  filter(obs_date == "2021-09-21") %>%
  select(Scenario, fc_date, obs_date, Bluegreens_ugL_obs, Bluegreens_ugL_pred)

for (i in 1:length(scenarios)){
  fc <- read_csv(paste0("./4_Data_assimilation/hindcasts/",scenarios[i],"/phytos_2021-10-28.csv")) %>%
    add_column(scenario = scenarios[i]) %>%
    select(scenario, fc_date, Bluegreens_ugL) %>%
    group_by(scenario, fc_date) %>%
    summarize(Bluegreens_ugL_med = median(Bluegreens_ugL),
              Bluegreens_ugL_upper = quantile(Bluegreens_ugL, probs = c(0.975)),
              Bluegreens_ugL_lower = quantile(Bluegreens_ugL, probs = c(0.025)))
  
  if(i == 1){
    example_fc <- fc
  } else {
    example_fc <- rbind(example_fc, fc)
  }
  
}

# example_fc <- example_fc %>%
#   gather(Bluegreens_ugL_med:Bluegreens_ugL_lower, key = "PI", value = "log_ugL")

my.obs <- appended_fc %>%
  filter(obs_date == "2021-11-03")
my.obs <- my.obs[c(11),]

ggplot(data = example_fc, aes(x = fc_date, y = Bluegreens_ugL_med, color = scenario)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Bluegreens_ugL_lower, ymax = Bluegreens_ugL_upper, group = scenario, color = scenario, fill = scenario), alpha = 0.2) +
  geom_point(data = my.obs, aes(x = fc_date, y = Bluegreens_ugL_obs), color = "black", size = 3) +
  labs(x = "", y = expression(paste("log ","(",mu,"g ",L^-1,")")), title = "") +
  theme_classic()

which.min(abs(appended_fc$GreenAlgae_ugL_pred - appended_fc$GreenAlgae_ugL_obs))
check <- appended_fc %>%
  filter(obs_date == "2021-09-21")

for (i in 1:length(scenarios)){
  fc <- read_csv(paste0("./4_Data_assimilation/hindcasts/",scenarios[i],"/phytos_2021-09-21.csv")) %>%
    add_column(scenario = scenarios[i]) %>%
    select(scenario, fc_date, GreenAlgae_ugL) %>%
    group_by(scenario, fc_date) %>%
    summarize(med = median(GreenAlgae_ugL),
              upper = quantile(GreenAlgae_ugL, probs = c(0.975)),
              lower = quantile(GreenAlgae_ugL, probs = c(0.025)))
  
  
  if(i == 1){
    example_fc2 <- fc
  } else {
    example_fc2 <- rbind(example_fc2, fc)
  }
  
}

my.obs <- appended_fc %>%
  filter(obs_date == "2021-09-21")
my.obs <- my.obs[c(1),]

ggplot(data = example_fc2, aes(x = fc_date, y = med, color = scenario)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, group = scenario, color = scenario, fill = scenario), alpha = 0.2) +
  geom_point(data = my.obs, aes(x = obs_date, y = GreenAlgae_ugL_obs), color = "black", size = 3) +
  labs(x = "", y = expression(paste("log ","(",mu,"g ",L^-1,")")), title = "") +
  theme_classic()
