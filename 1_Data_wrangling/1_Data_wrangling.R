##DATA WRANGLING TO GET OBSERVATION INPUT FILE FOR PHYTO DATA PRODUCT
##Author: Mary Lofton
##Date: 08FEB22

##To-do in this script:####
#1. Pull in all relevant data products from EDI (EXO, CTD, FP, counts)
#2. Aggregate to daily timestep between 1-2 m for 2016-2019
#3. Write to file

##SET-UP####
pacman::p_load(tidyverse, lubridate, data.table)

##EXO####
data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.271.6&entityid=23a191c1870a5b18cbc17f2779f719cf"

download.file(data, destfile = "./00_Data_files/FCR_Catwalk_2018_2021.csv", method='libcurl')

exo <- fread("./00_Data_files/FCR_Catwalk_2018_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  select(Date,DateTime,EXOChla_ugL_1,EXOBGAPC_ugL_1,Flag_Chla_ugL,Flag_Phyco_ugL) %>%
  filter(year(Date) %in% c(2014:2020))%>%
  group_by(Date) %>%
  summarize(daily_EXOChla_ugL_1 = mean(EXOChla_ugL_1, na.rm = TRUE),
            daily_EXOBGAPC_ugL_1 = mean(EXOBGAPC_ugL_1, na.rm = TRUE))

ggplot(data = exo, aes(x = Date, y = daily_EXOChla_ugL_1))+
  geom_line(size = 1, color = "darkgreen")+
  theme_bw()

ggplot(data = exo, aes(x = Date, y = daily_EXOBGAPC_ugL_1))+
  geom_line(size = 1, color = "cyan4")+
  theme_bw()

ggplot(data = exo, aes(x = daily_EXOChla_ugL_1, y = daily_EXOBGAPC_ugL_1))+
  geom_point(size = 1)+
  theme_bw()

##CTD####
data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.200.11&entityid=d771f5e9956304424c3bc0a39298a5ce"

download.file(data, destfile = "./00_Data_files/CTD_final_2013_2020.csv", method='libcurl')

ctd <- fread("./00_Data_files/CTD_final_2013_2020.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(Date)) %>%
  filter(year(Date) %in% c(2014:2020) & Reservoir == "FCR" & Site == 50 & Depth_m >=1 & Depth_m <=2)%>%
  select(Date,Chla_ugL,Flag_Chla) %>%
  group_by(Date) %>%
  summarize(integrated_Chla_ugL = mean(Chla_ugL, na.rm = TRUE))

ggplot(data = ctd, aes(x = Date, y = integrated_Chla_ugL))+
  geom_line(size = 1, color = "darkgreen")+
  theme_bw()

##FP####
data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.272.6&entityid=6b3151c0fdd913e02641363c2b00ae57"

download.file(data, destfile = "./00_Data_files/FluoroProbe_2014_2021.csv", method='libcurl')

fp <- fread("./00_Data_files/FluoroProbe_2014_2021.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(DateTime)) %>%
  filter(year(Date) %in% c(2014:2020) & Reservoir == "FCR" & Site == 50 & Depth_m >=1 & Depth_m <=2) %>%
  select(Date,GreenAlgae_ugL,Bluegreens_ugL,BrownAlgae_ugL,MixedAlgae_ugL,TotalConc_ugL,Flag_GreenAlgae,Flag_BluegreenAlgae,Flag_BrownAlgae,Flag_MixedAlgae,Flag_TotalConc) %>%
  group_by(Date) %>%
  summarize(integrated_GreenAlgae_ugL = mean(GreenAlgae_ugL, na.rm = TRUE),
            integrated_Bluegreens_ugL = mean(Bluegreens_ugL, na.rm = TRUE),
            integrated_BrownAlgae_ugL = mean(BrownAlgae_ugL, na.rm = TRUE),
            integrated_MixedAlgae_ugL = mean(MixedAlgae_ugL, na.rm = TRUE),
            integrated_TotalConc_ugL = mean(TotalConc_ugL, na.rm = TRUE))

ggplot(data = fp, aes(x = Date, y = integrated_GreenAlgae_ugL))+
  geom_line(size = 1, color = "green")+
  theme_bw()

ggplot(data = fp, aes(x = Date, y = integrated_Bluegreens_ugL))+
  geom_line(size = 1, color = "cyan4")+
  theme_bw()

ggplot(data = fp, aes(x = Date, y = integrated_BrownAlgae_ugL))+
  geom_line(size = 1, color = "darkgoldenrod3")+
  theme_bw()

ggplot(data = fp, aes(x = Date, y = integrated_MixedAlgae_ugL))+
  geom_line(size = 1, color = "brown")+
  theme_bw()

ggplot(data = fp, aes(x = Date, y = integrated_TotalConc_ugL))+
  geom_line(size = 1, color = "darkgreen")+
  theme_bw()

##COUNTS####
data <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.875.1&entityid=00d64d23fc2b75d973f69cc09bb5d083"

download.file(data, destfile = "./00_Data_files/phytoplankton.csv", method='libcurl')

counts <- fread("./00_Data_files/phytoplankton.csv", fill = TRUE, blank.lines.skip = FALSE) %>%
  mutate(Date = date(Date)) %>%
  filter(Depth_m <= 2 & Site == 50)

length(unique(counts$Date))

#calculate total BV for each sample day
total_bv <- counts %>%
  mutate(BV_um3mL = as.double(BV_um3mL)) %>%
  group_by(Date) %>%
  summarize(BV_TOTAL = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Date)) 

#spectral groups
#set up vectors of different divisions
cyano <- c("Pseudanabaena","Aphanocapsa","Dolichospermum","Synechococcus",
           "Microcystis","unicell Cyanobacterium","Woronichinia","Aphanocapsa",
           "Merismopedia","Snowella","Dactylococcopsis","Oscillatoria" )
chloro <- c("Chlorophyte sp. 1","Ankistrodesmus" ,"Chlorophyte sp. 2",
            "Oocystis","Chlorophyte sp. 6","Selenastrum","Chlorophyte sp. 4",
            "Dictyosphaerium" ,"Schroederia","unicell Eukaryote","Tetraselmis",
            "Dysmorphococcus","Chlamydomonas","Nephroselmis","Elakatothrix",
            "Pandorina","Eudorina","Oonephris","Monomastix","Platydorina",
            "Pediastrum","Micractinium","Gloeocystis","Chlorophyte sp. 3",
            "Botryococcus","Kirchneriella","Quadrigula","Actinastrum",
            "Tetraspora","Nephrocytium","Dunaliella","Coccomyxa","Volvulina",
            "Polytomella")
baci <- c("Nitzchia","Synedra","Asterionella","Cyclotella" ,"Fragilaria",
          "Ceratoneis","Tabellaria","Navicula" )
chryso <- c("Synura","Dinobryon","Oochromonas")
dino <- c("Gymnodinium","dinoflagellate cyst","Peridinium","naked dino","Parvodinium","Gloeodinium",
          "Prorocentrum")
desmid <- c("Spondylosium","Staurastrum","Staurodesmus","Closterium",
            "Desmid sp. 1","Bambusina","Actinotaenium","Cosmarium" )
crypto <- c("Cryptomonas","Rhodomonas" )
eugleno <- c("Trachelomonas","Euglena","Euglena?","Lepocinclis" )
raphid <- c("Gonyostomum")

genera <- c(cyano, chloro, baci, chryso, dino, desmid, crypto, eugleno,raphid)
check <- counts %>%
  filter(!Genus %in% genera)
bad_genera <- unique(check$Genus)

#get relative abundances of divisions
counts1 <- counts %>%
  mutate(Phyto_group = ifelse(Genus %in% cyano,"Cyanobacteria",
                              ifelse(Genus %in% chloro,"Chlorophytes",
                                     ifelse(Genus %in% baci, "Bacillaria",
                                            ifelse(Genus %in% chryso, "Chrysophytes",
                                                   ifelse(Genus %in% dino, "Dinoflagellates",
                                                          ifelse(Genus %in% desmid,"Desmids",
                                                                 ifelse(Genus %in% crypto, "Cryptophytes",
                                                                        ifelse(Genus %in% eugleno, "Euglenoids",
                                                                               ifelse(Genus %in% raphid, "Raphids",""))))))))))

counts2 <- counts1 %>%
  group_by(Date, Phyto_group) %>%
  summarize(BV_group = sum(BV_um3mL, na.rm = TRUE)) %>%
  mutate(Date = as.Date(Date))

dates <- unique(counts2$Date)
groups <- unique(counts2$Phyto_group)
combinations <- expand.grid(Date = dates, Phyto_group = groups)

counts3 <- full_join(counts2, combinations, by = c("Date","Phyto_group")) %>%
  mutate(BV_group = ifelse(is.na(BV_group), 0, BV_group)) 

counts4 <- counts3 %>%
  ungroup() %>%
  select(Date, BV_group, Phyto_group) %>%
  spread(BV_group,key = Phyto_group)

counts5 <- left_join(counts4, total_bv, by = "Date") 


colnames(counts5)[2:9] <- paste("BV", colnames(counts5)[2:9], sep = "_")

counts6 <- counts5 %>%
  mutate(BV_green = BV_Chlorophytes + BV_Desmids,
         BV_cyano = BV_Cyanobacteria,
         BV_diatom = BV_Bacillaria,
         BV_flagellates = BV_Dinoflagellates + BV_Cryptophytes) %>%
  select(Date, BV_green,BV_cyano,BV_diatom,BV_flagellates)

##COMBINE DATA STREAMS####