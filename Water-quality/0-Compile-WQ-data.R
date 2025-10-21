rm(list=ls());graphics.off();rm(.SavedPlots);gc();windows(record=T)
library(sp)
library(sf)
library(rgdal)
library(dplyr)
library(ggplot2)
library(ggmap)

##################################################################
##
## READ IN WATERQUALITY DATA SOURCES
## (and handle dates)

########################
## FDACS
fdacs = read.csv("./data/fdacs_wq_CK-SS-HB.csv")
names(fdacs) = tolower(names(fdacs))
fdacs$x = NULL
fdacs$date = as.Date(fdacs$date, format = ("%m/%d/%Y"))
fdacs$year = as.numeric(format(fdacs$date, format = "%Y"))
fdacs$month = as.character(format(fdacs$date, format = "%m"))
fdacs$ym = paste(fdacs$year, fdacs$month, sep = "-")
fdacs$day = as.numeric(format(fdacs$date, format = "%d"))

fdacs = fdacs %>% 
  filter (status == "ACTIVE") %>%  
  filter (year   >= 1997); fdacs$status = NULL

########################
## LCR
lcr = read.csv("./data/lcr_wq_total.csv")
names(lcr) = tolower(names(lcr))
lcr$x = NULL; lcr$id=NULL
lcr = lcr %>% 
  filter(in_service == 1) %>% 
  filter(between(site, 1, 10)) %>% 
  filter(between(salinity, 0, 45))
lcr$in_service = NULL

lcr$obs_date = as.Date(lcr$obs_date, format = ("%m/%d/%Y"))
lcr$date_time = as.POSIXct(lcr$date_time, format= "%m/%d/%Y %H:%M")
lcr$read_time = round(lcr$date_time, "hour")

## Read in and merge sites
sites_lcr = read.csv("./data/lcr_sites_lat_long.csv")
names(sites_lcr) = tolower(names(sites_lcr))
lcr2 = merge(lcr, sites_lcr, by.x = "site")


########################
## Lakewatch LCR data
lakewatch = read.csv("./data/lab.csv")
lakewatch$Sun_Code = NULL
lakewatch$Date = as.POSIXct(lakewatch$Date, format = "%m/%d/%Y %H:%M")
lakewatch$read_time = round(lakewatch$Date, "hour")
lakewatch$obs_date = as.Date(lakewatch$Date, format= "%Y/%m/%d")
lakewatch = lakewatch %>% filter (obs_date <= "2020-12-31")


###############################################
## Try to find missing salinity and temperature 
## measurements for Lakewatch data

## First remove YSI from "lab" dataset and join columns
lakewatch_only = subset(lakewatch, lakewatch$Sensor_Type=="LAKEWATCH"); nrow(lakewatch_only)
ysi_only = subset(lakewatch, lakewatch$Sensor_Type=="YSI"); nrow(ysi_only)
subysi = data.frame(read_time=ysi_only$read_time, Site=ysi_only$Site, 
                    Salinity =ysi_only$Salinity,  Temperature=ysi_only$Temperature)
lakewatch2 = lakewatch_only
lakewatch2$Salinity = NULL; lakewatch2$Temperature = NULL
lakewatch2 = merge(lakewatch2, subysi, by=c("Site", "read_time"), all=TRUE)
lakewatch2$Sensor_Type[is.na(lakewatch2$Sensor_Type)] = "YSI only"
nrow(lakewatch)-nrow(lakewatch2) ## 115 observations merged

## Next Join with LCR data by hour
sublcr = data.frame(read_time = lcr$read_time, Site = lcr$site, Sal_LCR = lcr$salinity, Temp_LCR = lcr$temperature)
lakewatch3 = left_join(lakewatch2, sublcr, by=c("Site", "read_time"))
lakewatch3$Diff_Sal =  lakewatch3$Salinity - lakewatch2$Sal_LCR
lakewatch3$Diff_Temp = lakewatch3$Temperature - lakewatch3$Temp_LCR

## Read in and merge sites
sites_lcr = read.csv("./data/lcr_sites_lat_long.csv")
lakewatch4 = merge(lakewatch3, sites_lcr, by.x = "Site"); head(lakewatch4)

write.csv(lakewatch4, "./data/lakewatch_sal-temp_joined.csv", row.names = F)

########################
## Read in Frazer data
frazer <- read.csv("./data/frazer_suwannee_97-15.csv")

frazer$Date = as.Date(frazer$Date, format="%Y-%m-%d") 
frazer %>% group_by(Station) %>% 
  summarise(Lat=mean(Lat), Long=mean(Long)) %>% 
  as.data.frame()
frazer$Station = as.factor(frazer$Station)

########################
## FIM physical data set
fim = read.csv("./data/FIM CK full physical dataset 2020.csv")

## Get dates
fim$date = as.Date(fim$date)
fim$Year = as.numeric(format(fim$date, format="%Y"))
fim$Month = as.numeric(format(fim$date, format="%m"))
fim$Day = as.numeric(format(fim$date, format="%d"))
#fim %>% group_by(Year, Month) %>% summarise(Sal = mean(salinity, na.rm=T)) %>% as.data.frame()

########################
## USF Optical virtual buoys
vbuoys = read.csv("./data/virtual_buoy_temps_compiled.csv", stringsAsFactors = T)
vbuoys$Date = as.Date(paste(vbuoys$YEAR, vbuoys$MONTH, 15, sep = "-"), format = "%Y-%m-%d")



##################################################################
##
## MAKE COMPILED DATA SET
##
##################################################################

## FDACS
Source = "FDACS"
attach(fdacs)
x.fdacs <- data.frame(Source = Source, Date = as.Date(date), Time = time,
                      Year = NA, Month = NA, YM = NA, 
                      Lat = lat, Long = long,
                      Site_ID = paste(Source, station_num, sep="-"),
                      Salinity = sal, Temperature = temp, 
                      DO = bdo, FC = fc, 
                      Chl = as.numeric(NA), TN = as.numeric(NA), TP = as.numeric(NA))
ncol(x.fdacs); nrow(x.fdacs); head(x.fdacs); str(x.fdacs)
detach(fdacs)

## LCR
Source = "LCR"


attach(lcr2)
x.lcr <- data.frame(Source = Source, Date = obs_date, Time = format(date_time, format = "%H:%M"),
                    Year = NA, Month = NA, YM = NA,
                    Lat = lat, Long = long,
                    Site_ID = paste(Source, site, sep="-"),
                    Salinity = salinity, Temperature = temperature, 
                    DO = as.numeric(NA), FC = as.numeric(NA),
                    Chl = as.numeric(NA), TN = as.numeric(NA), TP = as.numeric(NA))
str(x.lcr)
detach(lcr)

## Lakewatch
attach(lakewatch4)
Source = "Lakewatch"
x.lakewatch <- data.frame(Source = Source, Date = Date, Time = format(read_time, format = "%H:%M"),
                          Year = NA, Month = NA, YM = NA,
                          Lat = Lat, Long = Long,
                          Site_ID = paste(Source, Site, sep="-"),
                          Salinity = Salinity, Temperature = Temperature, 
                          DO = DO, FC = as.numeric(NA),
                          Chl = Chlorophyll, TN = Nitrogen, TP = Phosphorus)
detach(lakewatch)
nrow(x.lakewatch); head(x.lakewatch); str(x.lakewatch)

## Frazer
attach(frazer)
Source = "Frazer"
x.frazer <- data.frame(Source = Source, Date = Date, Time = NA, Year = NA, Month = NA, YM = NA,
                       Lat = Lat, Long = Long,
                       Site_ID = paste(Source, Station, sep="-"),
                       Salinity = Sal, Temperature = Temp_C, 
                       DO = as.numeric(DO_mg_L),FC = as.numeric(NA),
                       Chl = Chl, TN = TN, TP = TP)
nrow(x.frazer); head(x.frazer)
str(x.frazer)
detach(frazer)

## FIM
attach(fim)
Source = "FIM"
x.fim <- data.frame(Source = Source, Date = date, Time = NA, Year = NA, Month = NA, YM = NA,
                    Lat = latitude, Long = longitude,
                    Site_ID = reference,
                    Salinity = salinity, Temperature = temperature, 
                    DO = dissolvedO2, FC = as.numeric(NA),
                    Chl = as.numeric(NA), TN = as.numeric(NA), TP = as.numeric(NA))
detach(fim)
head(x.fim); tail(x.fim); nrow(x.fim)

## Virtual buoys
attach(vbuoys)
Source = "VBuoys"
x.vbuoy <- data.frame(Source = Source, Date = Date, Time = NA, Year = NA, Month = NA, YM = NA,
                      Lat = Lat, Long = Long,
                      Site_ID = Station,
                      Salinity = as.numeric(NA), Temperature = MEAN, 
                      DO = as.numeric(NA), FC = as.numeric(NA),
                      Chl = as.numeric(NA), TN = as.numeric(NA), TP = as.numeric(NA))
detach(vbuoy)
head(x.vbuoy); tail(x.vbuoy); nrow(x.vbuoy)


#################################
## Make compiled data set

## Rowbind dataframes
physcomp = rbind(x.fdacs, x.lcr, x.lakewatch, x.frazer, x.fim, x.vbuoy, x.vbuoy) 

## Classify factors
physcomp$Date = as.Date(physcomp$Date)
#physcomp$Date = as.POSIXct(physcomp$Date, format = "%m/%d/%Y %H:%M")
physcomp$Source = as.factor(physcomp$Source); physcomp$Site_ID = as.factor(physcomp$Site_ID)
physcomp$DO = as.numeric(physcomp$DO); physcomp$Salinity = as.numeric(physcomp$Salinity); physcomp$Temperature = as.numeric(physcomp$Temperature)
physcomp$Chl = as.numeric(physcomp$Chl); physcomp$TN = as.numeric(physcomp$TN); physcomp$TP = as.numeric(physcomp$TP)

## Make Year, Month, YM, Date
physcomp$Year = as.numeric(format(physcomp$Date, format="%Y"))
physcomp$Month = as.numeric(format(physcomp$Date, format="%m"))
physcomp$YM = as.factor(substr(physcomp$Date, 1, 7))
physcomp$TNP = physcomp$TN + physcomp$TP

## Filter and arrange by date
physcomp = physcomp %>% 
  filter (Date >= "1997-01-01" & Date <= "2020-12-31") %>% 
  arrange(Date, Time)

plot(Temperature ~ Date, data = subset(physcomp, Site_ID == "LCR-3")) ## Check a site

##########################################################################################
##
## AGGREGATE WITH FLOW DATA
##
##########################################################################################


##################################################################
## Get daily lagged flow rates
## Daily discharge from Wilcox station

## Get daily water data
library(waterData)
dailyflow = importDVs(staid="02323500",code="00060",sdate="1996-10-01",edate="2020-12-31")
dailyflow$staid = NULL; dailyflow$qualcode = NULL
names(dailyflow) = c("Flow", "Date"); head(dailyflow); tail(dailyflow)
dailyflow$Date = as.Date(dailyflow$Date)

## Calculate moving average
library(zoo)
dflow = dailyflow %>% 
  mutate(ma01=rollapply(Flow, 1,  mean, align='right', fill=NA),
         ma05=rollapply(Flow, 5,  mean, align='right', fill=NA),
         ma15=rollapply(Flow, 15,  mean, align='right', fill=NA),
         ma30=rollapply(Flow, 30, mean, align='right', fill=NA),
         sd05=rollapply(Flow, 5,  sd, align='right', fill=NA),
         sd15=rollapply(Flow, 15,  sd, align='right', fill=NA),
         sd30=rollapply(Flow, 30,  sd, align='right', fill=NA)
  ) 

dflow$Flow = NULL
dflow = dflow %>% filter(Date >= "1997-01-01") 
dflow[,5:10] = round(dflow[,5:10],1)
head(dflow, n = 30L)


#######################################
## Data set with lags
## Create columns for each day
lagflow = dailyflow
nlags = 45
s = ncol(lagflow)+1

for (j in 0:nlags) {
  lagflow$lagflow = lag(lagflow$Flow, n = j)
  names(lagflow)[s+j] = paste0("lag", j)
  lagflow$lagflow = NULL
}
lagflow = lagflow %>% filter(Date >= "1997-01-01"); lagflow$Flow = NULL; head(lagflow, n=30L)

## Write dataset with moving averages and lags
write.csv(merge(dflow, lagflow, by.x="Date"), "./out/Flow-with-avgs-lags.csv", row.names = F)

## Plot visual
plot(ma01 ~ Date, data = dflow, ylab="Daily Flow") 
lines(ma15 ~ Date, data = dflow, col = "green")
lines(ma30 ~ Date, data = dflow, col = "blue")

##########################################################################################
## Merge physical data and flow
subflow = dflow[,1:5]
names(subflow) = c("Date", "ma01", "ma05", "ma15", "ma30"); head(sublow)
subflow[,2:5]=round(subflow[,2:5],1)
head(subflow)

## Classify factors
physcomp$Source = as.factor(physcomp$Source); physcomp$Site_ID = as.factor(physcomp$Site_ID)
physcomp$DO = as.numeric(physcomp$DO); physcomp$Salinity = as.numeric(physcomp$Salinity); physcomp$Temperature = as.numeric(physcomp$Temperature)
physcomp$Chl = as.numeric(physcomp$Chl); physcomp$TN = as.numeric(physcomp$TN); physcomp$TP = as.numeric(physcomp$TP)
physcomp$Date = as.Date(physcomp$Date)

physcomp2 = merge(physcomp, subflow, by="Date", all.x=T)

str(physcomp2)
head(physcomp2); tail(physcomp2)

## Write out data 
write.csv(physcomp2, "./out/Sptatial-temp_Phys_all-measurements-with-Flow.csv")

## Number observations by source
physcomp2 %>% group_by(Source) %>% summarise(n(), n_distinct(Year)) 

##########################################################################################
##
## AGGREGATE BY DAY AND MONTH
##
##########################################################################################


##################################################################
##
## Get YM Averages per Site

## First group/average multiple readings in a day
physcomp3 = physcomp %>% 
  group_by(Source, Site_ID, Date, YM, Year, Month, Lat, Long) %>% 
  summarise(
    n_samp = n(),
    Sal_mean = mean(Salinity, na.rm=TRUE), 
    Sal_min = min(Salinity), 
    Sal_max = max(Salinity),
    Temp_mean = mean(Temperature, na.rm=TRUE), 
    Temp_min = min(Temperature), 
    Temp_max = max(Temperature),
    DO = mean(DO, na.rm=TRUE),
    FC = mean(FC, na.rm=TRUE),
    Chl = mean(Chl, na.rm=TRUE),
    TN = mean(TN, na.rm=TRUE),
    TP = mean(TP, na.rm=TRUE),
    TNP = mean(TNP, na.rm=TRUE),
  )

physcomp3flow = merge(physcomp3, subflow, by="Date") ## Merge with flow data
write.csv(physcomp3flow, "./out/Spatial-temp_Phys-Flow_xDay.csv", row.names = FALSE)

## Now group by YM
physcomp4 = physcomp3 %>% 
  group_by(Source, Site_ID,  YM, Year, Month, Lat, Long) %>% 
  summarise(
    n_days = n(),
    n_samp = sum(n_samp),
    Salinity   = round(mean(Sal_mean, na.rm=TRUE), 1), 
    Temperature  = round(mean(Temp_mean, na.rm=TRUE), 1),
    DO    = round(mean(DO, na.rm=TRUE), 1),
    FC = mean(FC, na.rm=TRUE),
    Chl   = round(mean(Chl, na.rm=TRUE), 0),
    TN    = round(mean(TN, na.rm=TRUE), 0),
    TP    = round(mean(TP, na.rm=TRUE), 0),
    TNP   = round(mean(TNP, na.rm=TRUE), 0),
  )

## Merge with monthly flow
dailyflow$Date = as.Date(dailyflow$Date)
dailyflow = dailyflow %>% filter(Date >= "1997-01-01")
dailyflow$YM = as.factor(paste(format(dailyflow$Date, "%Y"), format(dailyflow$Date, "%m"), sep = "-"))
monthlyflow = dailyflow %>% group_by(YM) %>% summarise(Flow = mean(Flow, na.rm=T)); monthlyflow

physcomp4flow = merge(physcomp4, monthlyflow, by="YM")


## Total number sample numbers
physcomp %>% group_by(Source) %>% 
  summarise(n_samples = n(), Years = n_distinct(Year), Months = n_distinct(YM),
            Days = n_distinct(Date), Sites = n_distinct(Site_ID), 
            Salinity = n_distinct(Salinity)-1, Temperature = n_distinct(Temperature)-1,
            FC = n_distinct(FC)-1, TN = n_distinct(TN)-1)

## Total number samples (individual days)
physcomp3 %>% group_by(Source) %>% 
  summarise(n_samples = n(), Years = n_distinct(Year), Months = n_distinct(YM),
            Days = n_distinct(Date), Sites = n_distinct(Site_ID), 
            Salinity = n_distinct(Sal_mean)-1, Temperature = n_distinct(Temp_mean)-1,
            FC = n_distinct(FC)-1, TN = n_distinct(TN)-1)


## Summarize physical samples by YM
summYM = physcomp4 %>% group_by(YM) %>% 
  summarise(Temp = n_distinct(Temperature),
            Sal  = n_distinct(Salinity), 
            TNP  = n_distinct(TNP))

sum(summYM$Temp); min(summYM$Temp); max(summYM$Temp); mean(summYM$Temp); sd(summYM$Temp)
sum(summYM$Sal); min(summYM$Sal); max(summYM$Sal); mean(summYM$Sal); sd(summYM$Sal)
sum(summYM$TNP); min(summYM$TNP); max(summYM$TNP); mean(summYM$TNP); sd(summYM$TNP)


write.csv(physcomp4flow, "./out/Spatial-temp_Phys-Flow_xMonth.csv", row.names = FALSE)
##--> This will be used to make rster maps for Ecospace





