##########################################################################
#
# LITTLE OWL BIWEEKLY SURVIVAL ANALYSIS FROM TELEMETRY - DATA PREPARATION
#
##########################################################################
# written by Steffen Oppel, July 2023
# based on data by Marco Perrig

## inp data files were all set up to have simultaneous start
## juveniles fledged at different times, so this needs to be accounted for

## extracted data from 15 June onwards, because overwinter data run until 15 Jun

## added simple summary for manuscript on 16 Oct 2023

## UPDATED ON 22 NOV TO ALIGN RECAPTURE MATRIX (which had a mismatch between 0 effort observations and 0-recaps)

library(tidyverse)
library(data.table)
library(lubridate)
library(tidyverse)
library(geosphere)
filter<-dplyr::filter
select<-dplyr::select
library(MCMCvis)
library(RMark)
library(stringr)
library(readxl)
# install.packages("remotes")
# remotes::install_github("nFrechen/RgetDWDdata")
# library(RgetDWDdata)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN RAW DATA FROM EXCEL FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")

dat2009<-read_excel("data/EH Master_juveniles_Moggi.xlsx",sheet="EH 2009 daily") %>%
  filter(Adult==0) %>%
  select(-Adult,-`Ring Number`,-`Remove after last capture`,-`Hatching 2`,-Remarks) %>%
  rename(LastAliveDate=`Last known to be alive`,DeathDate=`Known to be dead`) %>%
  gather(key="Date",value="Status",-Individual, -Sex,-Treatment,-Hatching,-LastAliveDate,-DeathDate) %>%
  mutate(Status=as.numeric(Status), DeathDate=as.numeric(DeathDate), Date=as.numeric(Date)) %>%
  filter(!is.na(Status)) %>%
  filter(!is.na(Date)) %>%
  mutate(Date=as.Date(Date, origin = "1899-12-30"), DeathDate=as.Date(DeathDate, origin = "1899-12-30"), LastAliveDate=ymd(LastAliveDate), Hatching=ymd(Hatching)) %>%
  mutate(Day=yday(Date),Month=month(Date),OCC=week(Date),HatchDay=yday(Hatching),HatchOCC=week(Hatching)) %>%
  filter(Date>ymd("2009-05-15")) %>%
  filter(Date<ymd("2009-08-01")) %>%
  group_by(Individual, Sex,Treatment,HatchOCC,OCC) %>%
  summarise(STATE=max(Status, na.rm=T)) %>%
  spread(key=OCC, value=STATE)
  

dat2010<-read_excel("data/EH Master_juveniles_Moggi.xlsx",sheet="EH 2010 daily") %>%
  filter(Adult==0) %>%
  select(-Adult,-`Ring Number`,-`Remove after last capture`,-`Hatching 2`,-Remarks) %>%
  rename(LastAliveDate=`Last known to be alive`,DeathDate=`Known to be dead`) %>%
  gather(key="Date",value="Status",-Individual, -Sex,-Treatment,-Hatching,-LastAliveDate,-DeathDate) %>%
  mutate(Status=as.numeric(Status), DeathDate=as.numeric(DeathDate), Date=as.numeric(Date)) %>%
  filter(!is.na(Status)) %>%
  filter(!is.na(Date)) %>%
  mutate(Date=as.Date(Date, origin = "1899-12-30"), DeathDate=as.Date(DeathDate, origin = "1899-12-30"), LastAliveDate=ymd(LastAliveDate), Hatching=ymd(Hatching)) %>%
  mutate(Day=yday(Date),Month=month(Date),OCC=week(Date),HatchDay=yday(Hatching),HatchOCC=week(Hatching)) %>%
  filter(Date>ymd("2010-05-15")) %>%
  filter(Date<ymd("2010-08-01")) %>%
  group_by(Individual, Sex,Treatment,HatchOCC,OCC) %>%
  summarise(STATE=max(Status, na.rm=T)) %>%
  spread(key=OCC, value=STATE)

dat2011<-read_excel("data/EH Master_juveniles_Moggi.xlsx",sheet="EH 2011 daily") %>%
  filter(Adult==0) %>%
  select(-Adult,-`Ring Number`,-`Remove after last capture`,-`Hatching 2`,-Remarks) %>%
  rename(LastAliveDate=`Last known to be alive`,DeathDate=`Known to be dead`) %>%
  gather(key="Date",value="Status",-Individual, -Sex,-Treatment,-Hatching,-LastAliveDate,-DeathDate) %>%
  mutate(Status=as.numeric(Status), DeathDate=as.numeric(DeathDate), Date=as.numeric(Date)) %>%
  filter(!is.na(Status)) %>%
  filter(!is.na(Date)) %>%
  mutate(Date=as.Date(Date, origin = "1899-12-30"), DeathDate=as.Date(DeathDate, origin = "1899-12-30"), LastAliveDate=ymd(LastAliveDate), Hatching=ymd(Hatching)) %>%
  mutate(Day=yday(Date),Month=month(Date),OCC=week(Date),HatchDay=yday(Hatching),HatchOCC=week(Hatching)) %>%
  filter(Date>ymd("2011-05-15")) %>%
  filter(Date<ymd("2011-08-01")) %>%
  group_by(Individual, Sex,Treatment,HatchOCC,OCC) %>%
  summarise(STATE=max(Status, na.rm=T)) %>%
  spread(key=OCC, value=STATE)

keepinds<- bind_rows(dat2009,dat2010,dat2011) %>%
  gather(key="OCC",value="STATE",-Individual, -Sex,-Treatment,-HatchOCC) %>%
  group_by(Individual, Sex,Treatment,HatchOCC) %>%
  summarise(N=sum(STATE, na.rm=T)) %>%
  filter(N>0)

ALLDAT<- bind_rows(dat2011,dat2010,dat2009) %>%
  filter(Individual %in% keepinds$Individual)

dim(ALLDAT)


### create fortnightly capture history to match data frame with LIOWpf below

ALLDAT2w<-ALLDAT %>%
  gather(key="OCC",value="STATE",-Individual, -Sex,-Treatment,-HatchOCC) %>%
  filter(!is.na(as.numeric(STATE))) %>%
  #mutate(OCC=round(as.numeric(OCC)/2)*2) %>%
  mutate(OCC=ceiling(as.numeric(OCC)/2)) %>%
  mutate(HatchOCC=ceiling(as.numeric(HatchOCC)/2)) %>%
  group_by(Individual, Sex,Treatment,HatchOCC,OCC) %>%
  summarise(STATE=max(STATE, na.rm=T)) %>%
  spread(key=OCC, value=STATE, fill=0) %>%
  rename(bird_id=Individual)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN WEATHER DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### pre-prepared winter weather data - these are indexed by encounter occasion and year
wincov<-fread("data/LIOW_winter_covariates.csv")
head(wincov)

### PREPARE WINTER COVARIATES
allcov<-wincov %>% gather(key="variable", value="value",-occ,-year) %>%
  group_by(variable) %>%
  mutate(value=scale(value,center=F)[,1]) %>%   ## scale all variables for better estimation
  ungroup() %>%
  spread(key=occ, value=value) %>%
  arrange(variable,year)

### ADD 6 OCCASIONS FOR POST_FLEDGING PERIOD
allcov[,26:31]<-allcov[,3]
names(allcov)[26:31]<-paste("pf",seq(1:6),sep="")

### revisited on 16 Oct 2023 to include weather and rain from post-fledging period
## manually downloaded air temp and precip from: https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/
## manually downloaded snow cover data from: https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/

## extract following variables:
# mean.precip, total.precip -> PRECIP
# mean.temp,mean.temp.min,mean.temp.max -> AirTemp
# mean.temp.ground.min, day.below.zero.ground -> GroundTemp, but only available at 5 cm depth, so less prone to freezing?
# mean.snow.cover, day.snow.cover0,day.snow.cover3,day.snow.cover5	-> ??


### READ IN THE RAW DATA ###

#ground_temp<-read.table("data/weather/produkt_eb_stunde_19951004_20221231_04931.txt", header=T, sep=";") %>%
ground_temp<-read.table("data/weather/produkt_eb_stunde_19881101_20221231_04349.txt", header=T, sep=";") %>%
  select(STATIONS_ID,MESS_DATUM,V_TE005) %>%
  mutate(Date=ymd_h(MESS_DATUM)) %>%
  mutate(year=year(Date),DAY=yday(Date)) %>%
  filter(year %in% c(2009,2010,2011,2012)) %>%
  mutate(OCC=ceiling(week(Date)/2)) %>%
  rename(ground_temp=V_TE005) %>%
  filter(ground_temp>-999) %>%
  group_by(STATIONS_ID,OCC,year,DAY) %>%
  summarise(min.grd.tmp=min(ground_temp)) %>%
  ungroup() %>%
  mutate(frz.grd=ifelse(min.grd.tmp<0,1,0)) %>%
  group_by(STATIONS_ID,OCC,year) %>%
  summarise(start=min(DAY),mean.temp.ground.min=mean(min.grd.tmp),day.below.zero.ground=sum(frz.grd)) 

#air_temp<-read.table("data/weather/produkt_tu_stunde_19880101_20221231_04931.txt", header=T, sep=";") %>%
air_temp<-read.table("data/weather/produkt_tu_stunde_20040601_20221231_04349.txt", header=T, sep=";") %>%
  select(STATIONS_ID,MESS_DATUM,TT_TU) %>%
  mutate(Date=ymd_h(MESS_DATUM)) %>%
  mutate(year=year(Date),DAY=yday(Date)) %>%
  filter(year %in% c(2009,2010,2011,2012)) %>%
  mutate(OCC=ceiling(week(Date)/2)) %>%
  rename(air_temp=TT_TU) %>%
  filter(air_temp>-999) %>%
  group_by(STATIONS_ID,OCC,year,DAY) %>%
  summarise(min.air.tmp=min(air_temp),mean.air.tmp=mean(air_temp)) %>%
  ungroup() %>%
  group_by(STATIONS_ID,OCC,year) %>%
  summarise(start=min(DAY),mean.temp=mean(mean.air.tmp),mean.temp.min=mean(min.air.tmp),mean.temp.max=max(mean.air.tmp)) 

#precip<-read.table("data/weather/produkt_rr_stunde_19951004_20221231_04931.txt", header=T, sep=";") %>%
precip<-read.table("data/weather/produkt_rr_stunde_20040601_20221231_04349.txt", header=T, sep=";") %>%
  select(STATIONS_ID,MESS_DATUM,R1,WRTR) %>%
  mutate(Date=ymd_h(MESS_DATUM),DAY=yday(Date)) %>%
  mutate(year=year(Date)) %>%
  filter(year %in% c(2009,2010,2011,2012)) %>%
  mutate(OCC=ceiling(week(Date)/2)) %>%
  rename(rain=R1) %>%  ### if prec_type ==7 then this is snow
  mutate(snow=ifelse(WRTR==7,rain,ifelse(WRTR==8,rain*0.5,0))) %>%
  filter(rain>-999) %>%
  group_by(STATIONS_ID,OCC,year,DAY) %>%
  summarise(day_precip=sum(rain),snow_sum=sum(snow)) %>%
  ungroup() %>%
  group_by(STATIONS_ID,OCC,year) %>%
  summarise(start=min(DAY),mean.precip=mean(day_precip), total.precip=sum(day_precip), mean.snow=mean(snow_sum),total.snow=sum(snow_sum))

#snow<-read.table("data/weather/produkt_klima_tag_19530101_20221231_04931.txt", header=T, sep=";") %>%
snow<-read.table("data/weather/produkt_klima_tag_19871101_20221231_04349.txt", header=T, sep=";") %>%
  select(STATIONS_ID,MESS_DATUM,SHK_TAG) %>%
  mutate(Date=ymd(MESS_DATUM)) %>%
  mutate(year=year(Date)) %>%
  filter(year %in% c(2009,2010,2011,2012)) %>%
  mutate(OCC=ceiling(week(Date)/2)) %>%
  rename(snow=SHK_TAG) %>%  ### snow depth
  filter(snow>-999) %>%
  mutate(snowday1=ifelse(snow>0,1,0),snowday3=ifelse(snow>3,1,0),snowday5=ifelse(snow>5,1,0)) %>%
  group_by(STATIONS_ID,OCC,year) %>%
  summarise(day.snow.cover0=sum(snowday1),day.snow.cover3=sum(snowday3),day.snow.cover5=sum(snowday5),mean.snow.cover=mean(snow))


### COMBINE ALL WEATHER DATA
dim(ground_temp)
dim(air_temp)
dim(precip)
weather.data<- ground_temp %>% left_join(air_temp, by=c('start','STATIONS_ID','OCC','year')) %>%
                    left_join(precip, by=c('start','STATIONS_ID','OCC','year')) %>%
                    left_join(snow, by=c('STATIONS_ID','OCC','year')) %>%
                    gather(key='variable', value='value', -STATIONS_ID,-OCC,-year,-start) %>%
                    filter(!variable %in% c("mean.snow","total.snow"))


### ADJUST ENCOUNTER OCCASIONS TO START on 15 MAY (to match with CMR data)
yday(ymd("2009-05-15"))
week(ymd("2009-05-15"))
ceiling(week(ymd("2009-05-15"))/2)  
  
weather.cov.matrix<- 
  weather.data %>% ungroup() %>%
  mutate(ch.occ.prel=OCC-9) %>%
  mutate(ch.year=ifelse(ch.occ.prel<1,year-1,year)) %>%
  mutate(ch.occ=ifelse(ch.occ.prel<1,ch.occ.prel+max(ch.occ.prel)+9,ch.occ.prel)) %>%
  filter(ch.year>2008) %>%
  spread(key=variable, value=value) %>%
  arrange(ch.year,ch.occ)
summary(weather.cov.matrix)  

### CREATE A MATRIX SIMILAR TO ALLCOV ABOVE BUT WITH ACTUAL DATA

allcov.new<-weather.cov.matrix %>% 
  select(-STATIONS_ID,-OCC,-year,-start,-ch.occ.prel) %>%
  gather(key="variable", value="value",-ch.occ,-ch.year) %>%
  group_by(variable) %>%
  #mutate(value=scale(value,center=F)[,1]) %>%   ## scale all variables for better estimation - removed to facilitate easier backtransformation
  ungroup() %>%
  spread(key=ch.occ, value=value) %>%
  arrange(variable,ch.year)

dim(allcov.new)

### ADD 6 OCCASIONS FOR SUBSEQUENT BREEDING PERIOD FROM THE NEXT YEAR's START
allcov.new$'28'<-0
allcov.new$'29'<-0
allcov.new$'30'<-0
for (l in unique(allcov.new$variable)) {
  for (y in 2009:2011) {
    allcov.new[allcov.new$variable==l & allcov.new$ch.year==y,27:32]<-allcov.new[allcov.new$variable==l & allcov.new$ch.year==y+1,3:8]
  }
}
allcov.new<-allcov.new %>% filter(ch.year<2012)


### IDENTIFY MISSING VALUES AND IMPUTE THEM
apply(is.na(allcov.new),2,which)
allcov.new[c(6,9,12,18),27:32]<-0 ### missing snow data in June/July imputed with 0



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPARE DATA FROM MOGGI AND DWD TO ASSESS HOW SIMILAR THEY ARE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

testplot<-allcov.new %>% gather(key='OCC', value=value, -ch.year,-variable) %>%
  mutate(OCC=as.numeric(OCC)) %>%
  filter(OCC>6) %>%
  mutate(OCC=OCC-6) %>%
  filter(OCC<24) %>%
  rename(year=ch.year) 

allcov %>% gather(key='OCC', value=ORIGvalue, -year,-variable) %>%
  mutate(OCC=as.numeric(OCC)) %>%
  filter(!is.na(OCC)) %>%
  left_join(testplot, by=c("year","variable","OCC")) %>%
  
  ggplot(aes(x=ORIGvalue,y=value, colour=as.factor(year))) +
  geom_point(size=2) +
  facet_wrap(~variable) +
  theme(legend.position=c(0.85,0.15)) +
  xlab("Moggis Wetterdaten")+
  ylab("Steffens Wetterdaten")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM MARK INP DATA FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")

LIOWch<-convert.inp(inp.filename='data/1st year.inp',
                    group.df=data.frame(year=c("2009","2010","2011")),
                 covariates = c('Male MS group','Age on Aug 1st','oldest residual.weight','oldest residual.tarsus'),
                 use.comments = TRUE)
str(LIOWch)




############ read in post-fledging data #####
LIOWpf<-convert.inp(inp.filename='data/DOB post-fledging.inp',
                    group.df=data.frame(cohort=c('2009 / Unfed / Original' , '2010 / Fed / Exchanged' , '2010 / Fed / Original' , '2010 / Unfed / Exchanged' , '2010 / Unfed / Original' , '2011 / Fed / Exchanged' , '2011 / Fed / Original' , '2011 / Unfed / Exchanged' , '2011 / Unfed / Original' )),
                    covariates = c('Hatching date' , 'Brood size' , 'Rank' , 'Start feeding' , 'Residual weight' , 'Residual wing' , 'Residual tarsus' , 'Residual beak' , 'Relative residual weight' , 'Relative residual wing' , 'Relative residual tarsus' , 'Relative residual beak' , 'Male'),
                    use.comments = TRUE)
str(LIOWpf)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT COVARIATES AND MATCH THE TWO ENCOUNTER HISTORIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## EXTRACT BIRD ID
LIOWpf$bird_id<-str_trim(row.names(LIOWpf))
LIOWch$bird_id<-str_trim(row.names(LIOWch))


## EXTRACT YEAR AND FEEDING REGIME
dim(LIOWpf)
LIOWpf[,18:20] <- str_split_fixed(string=LIOWpf$cohort,pattern=" / ",n=3)
names(LIOWpf)[c(4,8,10,18:20)] <-c("hatch_date","residual.weight","residual.tarsus","year","feeding","origin")
names(LIOWch)[5]<-"age_dept"



### try and merge post-fledging datasets
LIOWpf<-LIOWpf %>% left_join(ALLDAT2w, by="bird_id") %>%
  #mutate(CH2=paste(`20`,`22`,`24`,`26`,`28`,`30`,sep=""))
  mutate(CH2=paste(`10`,`11`,`12`,`13`,`14`,`15`,sep=""))

LIOWpf %>%
  select(bird_id,ch,hatch_date,CH2)



#### MERGE ALL DATA FROM POST-FLEDGING AND REST OF YEAR

LIOW<-LIOWpf %>% select(bird_id,CH2,hatch_date,Male,year,residual.weight,residual.tarsus,feeding,origin) %>%
  left_join(LIOWch[,c(1,3,5,8)],by=c("bird_id","year")) %>%
  mutate(ch=ifelse(is.na(ch),"000000000000000000000000",ch)) %>%
  mutate(ch=paste(CH2,ch,sep="")) %>%
  select(-CH2)

dim(LIOW)
head(LIOW)


## INSPECT HATCH DATE DISTRIBUTION
## the EncHist are all structured to be starting at 1, but that is difficult to reconcile with time-varying covariates later on during the year
## I suspect that the reason for the weird age effects is the issue that age and time are not properly aligned
hist(LIOW$hatch_date)
LIOW %>% select(bird_id,hatch_date,ch) %>% filter(hatch_date < 10 | hatch_date >40)


## CALCULATE AGE AT DEPARTURE FOR IND THAT DID NOT SURVIVE
summary(LIOW$hatch_date)
plot(LIOW$age_dept,LIOW$hatch_date)
summary(lm(LIOW$age_dept~LIOW$hatch_date))
LIOW$age_dept==(92-LIOW$hatch_date)  ### this is the conversion for age

LIOW<-LIOW %>% mutate(age_dept=92-hatch_date)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE INPUT DATA FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N.occ<-nchar(LIOW$ch[1])

### EXTRACT COLUMN VALUES INTO numeric vectors
y <- data.frame(str_split_fixed(LIOW$ch, "", N.occ)) %>% mutate_at(1:N.occ,as.numeric)
CH<-as.matrix(y, dimnames=F)
year <- as.numeric(LIOW$year)-2008
feeding <- ifelse(LIOW$feeding=="Unfed",0,1)
#season<-c(rep(1,6),rep(2,10),rep(3,5),rep(4,2)) ## Dispersal x 6, Winter x 10, Incubation x 5, Brood rearing x 2
#season<-c(rep(1,6),rep(2,10),rep(3,7)) ## Dispersal x 6, Winter x 10, Breeding x 7 - CHANGED ON 14 SEPT BECAUSE MS specifies only 3 stages
season<-c(rep(1,6),rep(2,6),rep(3,10),rep(4,8)) ## Summer x 6, Autumn x 6, Winter x 10, Spring x 7 - CHANGED ON 22 NOV TO MATCH 30 OCCASIONS

winter<-ifelse(season==3,1,0) ## binary variable for winter 

age <- LIOW[,10]   # age in days on 1 Aug
agemat<-CH
agemat[,7]<-age
for(col in 6:1){
 agemat[,col]<-agemat[,col+1]-14
}
for(col in 8:N.occ){
  agemat[,col]<-agemat[,col-1]+14
}

age_scale<-scale(agemat,center=F)  ## to prevent column-specific centering we need to set center=F
simpleage_scale<-scale(age)  ## only use age on 1 Aug as offset rather than temporal progression
weight <- LIOW[,5] # residual weight (seems to be standardized already)
weight_scale <- scale(LIOW[,5]) # scaled residual weight corrected for age at weighing and growth curve
size <- LIOW[,6] # residual tarsus (seems to be standardized already)
size_scale <- scale(LIOW[,6]) # scaled residual tarsus corrected for age at measurement

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


### PREPARE SEX COVARIATE
## 80 individuals do not have an assigned sex, so we use the overall proportion and then randomly allocate birds to a sex
sex <- LIOW[,3]
table(sex)
known.male.ratio<-table(sex)[4]/sum(table(sex)[c(1,4)])
sex[!(sex %in% c(0,1))]<-rbinom(n=sum(table(sex)[2:3]),size=1,prob=known.male.ratio)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE MATRIX FOR RECAPTURE PROBS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 - general recap prob
# 2 - limited effort: 15-18 in 2009, and 11 and 16 in 2010
# 3 - zero in 14, 19, 20 and 21 in 2009 because no field effort at all

## updated on 27 Sept by adding 7 encounter occasions
## changed on 12 Oct to 6 encounter occasions
## changed on 22 Nov to align with 0 captures

recap.mat<-matrix(1, nrow=nrow(CH),ncol=ncol(CH))

recap.mat[year==1,(c(15,16,17,18)+7)] <- 2
recap.mat[year==2,(c(11,16)+7)] <- 2
recap.mat[year==1,(c(14,19,20,21)+7)] <- 3

## check that 0 effort and 0 sightings are in same columns
which(apply(CH[LIOW$year==2009,],2,sum)==0)
which(apply(recap.mat[LIOW$year==2009,],2,max)==3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE WORKSPACE AND R ENVIRONMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image("data/LIOW_SURV_INPUT.RData")
#renv::init()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NUMBERS NEEDED FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dim(LIOW)[1] ## number of individuals
length(unique(gsub("\\..*","",LIOWpf$bird_id))) ### number of broods

### CALCULATE PRIMITIVE SURVIVAL AS % OF INDIVIDUALS RECORDED IN LAST OCCASION - important for first line in manuscript
sum(y[,dim(y)[2]])/dim(y)[1]

