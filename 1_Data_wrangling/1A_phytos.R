##DATA WRANGLING TO GET PHYTO OBSERVATION INPUT FILE FOR REGULARIZED REGRESSION
##Author: Mary Lofton
##Date: 08FEB22

##To-do in this script:####
#1. Pull in all relevant data products from EDI (EXO, CTD, FP)
#2. Retrieve relevant measurements (0.5, 1, 1.5 m until 2019-05-20, then
#   1, 1.6, 2 m until 2021-12-31); all EXO measurements from 2018-2021
#3. Collate data from different sources and write to file

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

#Write a function that returns the closest value
#xv is vector, sv is specific value
closest<-function(xv, sv){
  xv[which.min(abs(xv-sv))]}
closest_condition<-function(df, colname1, colname2,sv){
  df[[colname1]][df[[colname2]]==closest(df[[colname2]],sv)]
}


##EXO####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.271.6&entityid=23a191c1870a5b18cbc17f2779f719cf"
# 
# download.file(data, destfile = "./00_Data_files/FCR_Catwalk_2018_2021.csv", method='libcurl')

exo <- data.table::fread("./00_Data_files/FCR_Catwalk_2018_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  select(Date,DateTime,EXOChla_ugL_1,EXOBGAPC_ugL_1,Flag_Chla_ugL,Flag_Phyco_ugL) %>%
  filter(year(Date) %in% c(2014:2021))%>%
  group_by(Date) %>%
  summarize(daily_EXOChla_ugL_1 = mean(EXOChla_ugL_1, na.rm = TRUE),
            daily_EXOBGAPC_ugL_1 = mean(EXOBGAPC_ugL_1, na.rm = TRUE))

##CTD####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.200.12&entityid=0a62d1946e8d9a511bc1404e69e59b8c"
# 
# download.file(data, destfile = "./00_Data_files/CTD_final_2013_2021.csv", method='libcurl')

ctd <- data.table::fread("./00_Data_files/CTD_final_2013_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(Date)) %>%
  filter(year(Date) %in% c(2014:2021) & Reservoir == "FCR" & Site == 50)%>%
  select(Date,Chla_ugL,Flag_Chla,Depth_m) 

dates <- unique(ctd$Date)

for (i in 1:length(dates)){
  profile <- ctd %>%
    filter(Date == dates[i])
  
  if(dates[i] < "2019-05-20"){
  answer <- profile %>%
    summarize(CTDChla_ugL_above = mean(closest_condition(.,"Chla_ugL","Depth_m", 0.5),na.rm = TRUE),
           CTDChla_ugL = mean(closest_condition(.,"Chla_ugL","Depth_m", 1),na.rm = TRUE),
           CTDChla_ugL_below = mean(closest_condition(.,"Chla_ugL","Depth_m", 1.5),na.rm = TRUE)) %>%
    select(CTDChla_ugL_above, CTDChla_ugL, CTDChla_ugL_below)
  } else {
  answer <- profile %>%
    summarize(CTDChla_ugL_above = mean(closest_condition(.,"Chla_ugL","Depth_m", 1),na.rm = TRUE),
           CTDChla_ugL = mean(closest_condition(.,"Chla_ugL","Depth_m", 1.6),na.rm = TRUE),
           CTDChla_ugL_below = mean(closest_condition(.,"Chla_ugL","Depth_m", 2),na.rm = TRUE)) %>%
    select(CTDChla_ugL_above, CTDChla_ugL, CTDChla_ugL_below)
  }
  if(i==1){temp <- answer[0,]}
  temp <- rbind(temp,answer)
}

ctd1 <- cbind(dates,temp) %>%
  rename(Date = dates)

##FP####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.272.6&entityid=6b3151c0fdd913e02641363c2b00ae57"
# 
# download.file(data, destfile = "./00_Data_files/FluoroProbe_2014_2021.csv", method='libcurl')

fp <- data.table::fread("./00_Data_files/FluoroProbe_2014_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  filter(year(Date) %in% c(2014:2021) & Reservoir == "FCR" & Site == 50) %>%
  select(Date,Depth_m,GreenAlgae_ugL,Bluegreens_ugL,BrownAlgae_ugL,MixedAlgae_ugL,TotalConc_ugL,Flag_GreenAlgae,Flag_BluegreenAlgae,Flag_BrownAlgae,Flag_MixedAlgae,Flag_TotalConc) 

dates <- unique(fp$Date)

for (i in 1:length(dates)){
  profile <- fp %>%
    filter(Date == dates[i])
  
  if(dates[i] < "2019-05-20"){
    answer <- profile %>%
      summarize(GreenAlgae_ugL_above = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 0.5),na.rm = TRUE),
                GreenAlgae_ugL = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                GreenAlgae_ugL_below = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 1.5),na.rm = TRUE),
                Bluegreens_ugL_above = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 0.5),na.rm = TRUE),
                Bluegreens_ugL = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 1),na.rm = TRUE),
                Bluegreens_ugL_below = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 1.5),na.rm = TRUE),
                BrownAlgae_ugL_above = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 0.5),na.rm = TRUE),
                BrownAlgae_ugL = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                BrownAlgae_ugL_below = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 1.5),na.rm = TRUE),
                MixedAlgae_ugL_above = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 0.5),na.rm = TRUE),
                MixedAlgae_ugL = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                MixedAlgae_ugL_below = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 1.5),na.rm = TRUE)) 
  } else {
    answer <- profile %>%
      summarize(GreenAlgae_ugL_above = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                GreenAlgae_ugL = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 1.6),na.rm = TRUE),
                GreenAlgae_ugL_below = mean(closest_condition(.,"GreenAlgae_ugL","Depth_m", 2),na.rm = TRUE),
                Bluegreens_ugL_above = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 1),na.rm = TRUE),
                Bluegreens_ugL = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 1.6),na.rm = TRUE),
                Bluegreens_ugL_below = mean(closest_condition(.,"Bluegreens_ugL","Depth_m", 2),na.rm = TRUE),
                BrownAlgae_ugL_above = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                BrownAlgae_ugL = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 1.6),na.rm = TRUE),
                BrownAlgae_ugL_below = mean(closest_condition(.,"BrownAlgae_ugL","Depth_m", 2),na.rm = TRUE),
                MixedAlgae_ugL_above = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 1),na.rm = TRUE),
                MixedAlgae_ugL = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 1.6),na.rm = TRUE),
                MixedAlgae_ugL_below = mean(closest_condition(.,"MixedAlgae_ugL","Depth_m", 2),na.rm = TRUE)) 
    
  }
  if(i==1){temp <- answer[0,]}
  temp <- rbind(temp,answer)
}

fp1 <- cbind(dates,temp) %>%
  rename(Date = dates)

##COMBINE DATA STREAMS####

#date range: 2014-05-05 to 2021-12-31

#define function to get data product date series (thx Jake Zwart!)
get_model_dates = function(model_start, model_stop, time_step = 'days'){
  
  model_dates = seq.Date(from = as.Date(model_start), to = as.Date(model_stop), by = time_step)
  
  return(model_dates)
}
mydates <- data.frame(get_model_dates(model_start = "2014-05-05",
                model_stop = "2021-12-31",
                time_step = 'days'))
colnames(mydates) <- c("Date")
f1 <- left_join(mydates, ctd1, by = "Date")
f2 <- left_join(f1, exo, by = "Date")
f3 <- left_join(f2, fp1, by = "Date")
write.csv(f3, "./1_Data_wrangling/collated_phyto_data.csv",row.names = FALSE)
