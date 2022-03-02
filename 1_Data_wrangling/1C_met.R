##DATA WRANGLING TO GET MET OBSERVATION INPUT FILE FOR REGULARIZED REGRESSION
##Author: Mary Lofton
##Date: 02MAR22

##To-do in this script:####
#1. Pull in all relevant data products from EDI (met station)
#2. Aggregate measurements (air temp, rel hum, windspeed, shortwave,
#   longwave aggregated to daily)
#3. Write to file

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

##MET####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.389.6&entityid=a5524c686e2154ec0fd0459d46a7d1eb"
# 
# download.file(data, destfile = "./00_Data_files/FCR_Met_final_2015_2021.csv", method='auto')

met <- data.table::fread("./00_Data_files/FCR_Met_final_2015_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  select(Date,AirTemp_Average_C,RH_percent,Rain_Total_mm,WindSpeed_Average_m_s,WindDir_degrees,ShortwaveRadiationDown_Average_W_m2,InfraredRadiationDown_Average_W_m2,PAR_Average_umol_s_m2) %>%
  filter(year(Date) %in% c(2014:2021))%>%
  group_by(Date) %>%
  summarize(daily_AirTemp_C = mean(AirTemp_Average_C, na.rm = TRUE),
            daily_RH_percent = mean(RH_percent, na.rm = TRUE),
            daily_Rain_Total_mm = sum(Rain_Total_mm, na.rm = TRUE),
            daily_WindSpeed_Average_m_s = mean(WindSpeed_Average_m_s, na.rm = TRUE),
            daily_WindDir_degrees = mean(WindDir_degrees, na.rm = TRUE),
            daily_ShortwaveRadiationDown_Average_W_m2 = mean(ShortwaveRadiationDown_Average_W_m2, na.rm = TRUE),
            daily_InfraredRadiationDown_Average_W_m2 = mean(InfraredRadiationDown_Average_W_m2, na.rm = TRUE),
            daily_PAR_Average_umol_s_m2 = mean(PAR_Average_umol_s_m2, na.rm = TRUE))


##WRITE TO FILE####

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
f1 <- left_join(mydates, met, by = "Date")
write.csv(f1, "./1_Data_wrangling/collated_met_data.csv",row.names = FALSE)
