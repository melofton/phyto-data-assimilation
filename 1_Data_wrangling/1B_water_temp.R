##DATA WRANGLING TO GET WATER TEMP OBSERVATION INPUT FILE FOR REGULARIZED REGRESSION
##Author: Mary Lofton
##Date: 08FEB22

##To-do in this script:####
#1. Pull in all relevant data products from EDI (EXO, CTD, YSI)
#2. Retrieve relevant measurements (0.5, 1, 1.5 m until 2019-05-20, then
#   1, 1.6, 2 m until 2021-12-31); all EXO temp measurements from 2018-2021
#3. Collate various data sources and write to file

#Note: for dates before 2018-08-06 (when EXO was deployed), we will use water 
#temperature readings from the thermistors, the CTD and then the YSI in that 
#order of preference

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
  select(Date,DateTime,EXOTemp_C_1, ThermistorTemp_C_1) %>%
  filter(year(Date) %in% c(2014:2021))%>%
  group_by(Date) %>%
  summarize(daily_EXOTemp_C_1 = mean(EXOTemp_C_1, na.rm = TRUE),
            daily_ThermistorTemp_C_1 = mean(ThermistorTemp_C_1, na.rm = TRUE))

##CTD####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.200.12&entityid=0a62d1946e8d9a511bc1404e69e59b8c"
# 
# download.file(data, destfile = "./00_Data_files/CTD_final_2013_2021.csv", method='libcurl')

ctd <- data.table::fread("./00_Data_files/CTD_final_2013_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(Date)) %>%
  filter(year(Date) %in% c(2014:2021) & Reservoir == "FCR" & Site == 50)%>%
  select(Date,Temp_C,Depth_m) 

dates <- unique(ctd$Date)

for (i in 1:length(dates)){
  profile <- ctd %>%
    filter(Date == dates[i])
  
  if(dates[i] < "2019-05-20"){
  answer <- profile %>%
    summarize(CTDTemp_C = mean(closest_condition(.,"Temp_C","Depth_m", 1),na.rm = TRUE)) 
  } else {
  answer <- profile %>%
    summarize(CTDTemp_C = mean(closest_condition(.,"Temp_C","Depth_m", 1.6),na.rm = TRUE)) 
  }
  if(i==1){temp <- answer[0,]}
  temp <- rbind(temp,answer)
}

ctd1 <- cbind(dates,temp) %>%
  rename(Date = dates)

##YSI####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.198.10&entityid=b3bd353312f9e37ca392e2a5315cc9da"
# 
# download.file(data, destfile = "./00_Data_files/YSI_PAR_profiles_2014_2021.csv", method='libcurl')

ysi <- data.table::fread("./00_Data_files/YSI_PAR_profiles_2014_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  filter(year(Date) %in% c(2014:2021) & Reservoir == "FCR" & Site == 50 & Depth_m >=1 & Depth_m <= 2) %>%
  select(Date,Depth_m,Temp_C) %>%
  filter(complete.cases(.))

dates <- unique(ysi$Date)

for (i in 1:length(dates)){
  profile <- ysi %>%
    filter(Date == dates[i])
  
  if(dates[i] < "2019-05-20"){
    answer <- profile %>%
      summarize(YSITemp_C = mean(closest_condition(.,"Temp_C","Depth_m", 1),na.rm = TRUE)) 
  } else {
    answer <- profile %>%
      summarize(YSITemp_C = mean(closest_condition(.,"Temp_C","Depth_m", 1.6),na.rm = TRUE)) 
  }
  if(i==1){temp <- answer[0,]}
  temp <- rbind(temp,answer)
}

ysi1 <- cbind(dates,temp) %>%
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
f3 <- left_join(f2, ysi1, by = "Date")
f4 <- f3 %>%
  mutate(Temp_C = ifelse(!is.na(daily_EXOTemp_C_1),daily_EXOTemp_C_1,
                         ifelse(!is.na(daily_ThermistorTemp_C_1),daily_ThermistorTemp_C_1,
                                ifelse(!is.na(CTDTemp_C),CTDTemp_C,
                                       ifelse(!is.na(YSITemp_C),YSITemp_C,NA)))))
write.csv(f4, "./1_Data_wrangling/collated_water_temp_data.csv",row.names = FALSE)
