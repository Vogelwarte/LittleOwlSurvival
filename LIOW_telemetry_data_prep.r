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
setwd("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
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
  mutate(OCC=round(as.numeric(OCC)/2)*2) %>%
  group_by(Individual, Sex,Treatment,HatchOCC,OCC) %>%
  summarise(STATE=max(STATE, na.rm=T)) %>%
  spread(key=OCC, value=STATE, fill=0) %>%
  rename(bird_id=Individual)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN WEATHER DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### abandoned because difficult to replicate what Moggi did

# ## https://rdrr.io/github/nFrechen/RgetDWDdata/man/getDWDdata.html
# KlimadatenLB <- getDWDdata(Messstelle=04349, historisch=T)
# # lädt den aktuellen Datensatz der Station Cottbus herunter
# head(KlimadatenLB)
# tail(KlimadatenLB)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM MARK INP DATA FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")

LIOWch<-convert.inp(inp.filename='data/1st year.inp',
                    group.df=data.frame(year=c("2009","2010","2011")),
                 covariates = c('Male MS group','Age on Aug 1st','oldest residual.weight','oldest residual.tarsus'),
                 use.comments = TRUE)
str(LIOWch)


wincov<-fread("data/LIOW_winter_covariates.csv")
head(wincov)



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
  mutate(CH2=paste(`20`,`22`,`24`,`26`,`28`,`30`,sep=""))

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
season<-c(rep(1,6),rep(2,6),rep(3,10),rep(4,7)) ## Summer x 6, Autumn x 6, Winter x 10, Spring x 7 - CHANGED ON 27 SEPT BECAUSE WE NOW INCLUDE THE WHOLE YEAR

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
size <- LIOW[,6] # residual tarsus (seems to be standardized already)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


### PREPARE WINTER COVARIATES
allcov<-wincov %>% gather(key="variable", value="value",-occ,-year) %>%
  group_by(variable) %>%
  mutate(value=scale(value)[,1]) %>%   ## scale all variables for better estimation
  ungroup() %>%
  spread(key=occ, value=value) %>%
  arrange(variable,year)

### ADD 7 OCCASIONS FOR POST_FLEDGING PERIOD
allcov[,26:31]<-allcov[,3]
names(allcov)[26:31]<-paste("pf",seq(1:6),sep="")


### PREPARE SEX COVARIATE
## 80 individuals do not have an assigned sex, so we use the overall proportion and then randomly allocate birds to a sex
sex <- LIOW[,3]
table(sex)
known.male.ratio<-table(sex)[4]/sum(table(sex)[c(1,4)])
sex[!(sex %in% c(0,1))]<-rbinom(n=sum(table(sex)[2:3]),size=1,prob=known.male.ratio)


### CALCULATE PRIMITIVE SURVIVAL AS % OF INDIVIDUALS RECORDED IN LAST OCCASION - important for first line in manuscript
sum(y[,dim(y)[2]])/dim(y)[1]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE MATRIX FOR RECAPTURE PROBS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 - general recap prob
# 2 - limited effort: 15-18 in 2009, and 11 and 16 in 2010
# 3 - zero in 14, 19, 20 and 21 in 2009 because no field effort at all

## updated on 27 Sept by adding 7 encounter occasions
## changed on 12 Oct to 6 encounter occasions

recap.mat<-matrix(1, nrow=nrow(CH),ncol=ncol(CH))

recap.mat[year==1,(c(15,16,17,18)+6)] <- 2
recap.mat[year==2,(c(11,16)+6)] <- 2
recap.mat[year==1,(c(14,19,20,21)+6)] <- 3



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE WORKSPACE AND R ENVIRONMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image("data/LIOW_SURV_INPUT.RData")
#renv::init()


