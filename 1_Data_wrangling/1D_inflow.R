##DATA WRANGLING TO GET INFLOW OBSERVATION INPUT FILE FOR REGULARIZED REGRESSION
##Author: Mary Lofton
##Date: 02MAR22

##To-do in this script:####
#1. Pull in all relevant data products from EDI (weir)
#2. Aggregate measurements (WVWA gauge and VT gauge to daily)
#3. Write to file

#Note: using VT gauge when available and then WVWA gauge, in that order

##SET-UP####
rm(list = ls())
pacman::p_load(tidyverse, lubridate)

##INFLOW####
# data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.202.8&entityid=cc045f9fe32501138d5f4e1e7f40d492"
# 
# download.file(data, destfile = "./00_Data_files/Inflow_2013_2021.csv", method='libcurl')

inf <- data.table::fread("./00_Data_files/Inflow_2013_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  select(Date,WVWA_Flow_cms,VT_Flow_cms) %>%
  filter(year(Date) %in% c(2014:2021))%>%
  group_by(Date) %>%
  summarize(daily_WVWA_Flow_cms = mean(WVWA_Flow_cms, na.rm = TRUE),
            daily_VT_Flow_cms = mean(VT_Flow_cms, na.rm = TRUE))


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
f1 <- left_join(mydates, inf, by = "Date")
f2 <- f1 %>%
  mutate(Flow_cms = ifelse(!is.na(daily_VT_Flow_cms),daily_VT_Flow_cms,
                           ifelse(!is.na(daily_WVWA_Flow_cms),daily_WVWA_Flow_cms,NA)))
write.csv(f2, "./1_Data_wrangling/collated_inflow_data.csv",row.names = FALSE)
