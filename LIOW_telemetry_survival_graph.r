##########################################################################
#
# LITTLE OWL FIRST YEAR SURVIVAL PROJECTION
#
##########################################################################
# written by Steffen Oppel, July 2023
# based on data by Marco Perrig


## 24 Jan 2024: created plot after brief meeting with Tschumi and Murdy

library(data.table)
library(lubridate)
library(tidyverse)
filter<-dplyr::filter
select<-dplyr::select
library(scales)
library(janitor)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM PREPARED WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### data preparation moved to LIOW_telemetry_data_prep.r
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival"),silent=T)
load("LIOW_survival_output.RData")
head(MCMCpred)
seasons<-data.frame(season=c(rep(1,4),rep(2,6),rep(3,10),rep(4,7)), ## Summer x 4, Autumn x 6, Winter x 10, Spring x 6 - MATCH TABLE 1 in manuscript
                    Fortnight=seq(1,27,1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A BLANK DATA FRAME FOR PROJECTING SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### after initial screening decided to remove sex

AnnualSurvival<-crossing(Fortnight=seq(1,27,1),
                           feeding=c(0,1), 
                           #sex=c(0,1),
                           winter=c("mild","harsh")) %>%
  mutate(season=seasons$season[match(Fortnight,seasons$Fortnight)]) %>%
  mutate(snow=ifelse(winter=="mild",0,
                     ifelse(Fortnight %in% c(13,14,15),4,
                            ifelse(Fortnight %in% c(16,17,18),8,
                                   ifelse(Fortnight %in% c(19,20),12,0))))) %>%
  mutate(med.N=ifelse(Fortnight==1,100,NA),lcl.N=ifelse(Fortnight==1,100,NA),ucl.N=ifelse(Fortnight==1,100,NA))

AnnualSurvival %>% filter(winter=='harsh' & feeding==0 & season==3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOOP OVER EACH SCENARIO TO SAMPLE SURVIVORS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### after initial screening decided to remove sex

#for (s in c(0,1)){  ### different sexes
  for(f in c(0,1)) { ### different feeding regimes
    for(w in c("mild","harsh")){ ### snow conditions
      for(t in 1:27){  ### for every fortnight
        surv.est.fort<-MCMCpred %>%
          #filter(sex==s) %>%
          filter(feeding==f) %>%
          filter(season==seasons$season[t]) %>%
          filter(scaleweight==0) %>%
          filter(snow==ifelse(w=="mild",0,
                              ifelse(t %in% c(13,14,15),4,
                                     ifelse(t %in% c(16,17,18),8,
                                            ifelse(t %in% c(19,20),12,0))))) %>%
          summarise(med=quantile(surv,0.5),lcl=quantile(surv,0.025),ucl=quantile(surv,0.975))
        
        AnnualSurvival$med.N[#AnnualSurvival$sex==s &
                               AnnualSurvival$winter==w &
                               AnnualSurvival$feeding==f &
                               AnnualSurvival$Fortnight==t+1] <- AnnualSurvival$med.N[#AnnualSurvival$sex==s &
                                                                                               AnnualSurvival$winter==w &
                                                                                               AnnualSurvival$feeding==f &
                                                                                               AnnualSurvival$Fortnight==t]*surv.est.fort$med
        
        AnnualSurvival$lcl.N[#AnnualSurvival$sex==s &
                               AnnualSurvival$winter==w &
                               AnnualSurvival$feeding==f &
                               AnnualSurvival$Fortnight==t+1] <- AnnualSurvival$lcl.N[#AnnualSurvival$sex==s &
                                                                                                 AnnualSurvival$winter==w &
                                                                                                 AnnualSurvival$feeding==f &
                                                                                                 AnnualSurvival$Fortnight==t]*surv.est.fort$lcl
        
        AnnualSurvival$ucl.N[#AnnualSurvival$sex==s &
                               AnnualSurvival$winter==w &
                               AnnualSurvival$feeding==f &
                               AnnualSurvival$Fortnight==t+1] <- AnnualSurvival$ucl.N[#AnnualSurvival$sex==s &
                                                                                                 AnnualSurvival$winter==w &
                                                                                                 AnnualSurvival$feeding==f &
                                                                                                 AnnualSurvival$Fortnight==t]*surv.est.fort$ucl
        
      }
    }
  }
#}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE GRAPH OVER TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### CREATE PLOT

AnnualSurvival %>%
  mutate(Condition=ifelse(feeding==0,"natural","food supplemented")) %>%
  mutate(Date=ymd("2010-04-19")+weeks(Fortnight*2)) %>%

ggplot(aes(x=Date, y=med.N,ymin=lcl.N, ymax=ucl.N, colour=Condition, fill=Condition, linetype=Condition))+
  geom_line(linewidth=1.5) +
  geom_ribbon(alpha=0.4, colour=NA)+
  facet_wrap(~winter, ncol=2)+
  
  ## format axis ticks
  #scale_x_continuous(name="Season", limits=c(1,365), breaks=plotdat$age[c(3,5,8,10)], labels=plotdat$Season[c(3,5,8,10)]) +
  scale_x_date(name="Month",labels = date_format("%b"),breaks = "1 month")+
  scale_y_continuous(name="Proportion of fledglings remaining alive", limits=c(0,100), breaks=seq(0,100,20), labels=seq(0,100,20)) +

  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_line(color = "grey87"),
        # panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=16, color="black"),
        legend.background=element_rect(fill="white", colour="white"),
        legend.key = element_rect(fill = NA),
        legend.position=c(0.85,0.91), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))


# ggsave("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/MANUSCRIPTS/LittleOwlSurvival/Fig_2.jpg", height=7, width=11)
# ggsave("C:/STEFFEN/OneDrive - Vogelwarte/General - Little owls/MANUSCRIPTS/LittleOwlSurvival/Fig_2.jpg", height=7, width=11)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMARY TABLE HOW MANY % DIED PER SEASON
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Table2<-
AnnualSurvival %>% group_by(season, feeding, winter) %>%
  summarise(min=min(med.N), start=max(med.N)) %>%
  arrange(feeding, winter, season, start) %>%
  ungroup() %>%
  mutate(end=dplyr::lead(start, n=1)) %>%
  mutate(end=ifelse(is.na(end),min,end)) %>%
  mutate(end=ifelse(end==100,min,end)) %>%
  mutate(loss=start-end) %>%
  group_by(feeding,winter) %>%
  mutate(prop.loss=(loss/sum(loss))*100) %>%
  ungroup() %>%
  mutate(group=paste(feeding,winter,sep="_")) %>%
  select(season, group,prop.loss) %>%
  spread(key=group, value=prop.loss) %>%
  select(season,`0_mild`,`0_harsh`,`1_mild`,`1_harsh`) %>%
  janitor::adorn_totals()

#fwrite(Table2,"C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/MANUSCRIPTS/LittleOwlSurvival/Table2_surv.csv")
#fwrite(Table2,"C:/STEFFEN/OneDrive - Vogelwarte/General - Little owls/MANUSCRIPTS/LittleOwlSurvival/Table2_surv.csv")



##### SUMMARIES FOR TEXT ##########

finsurv<-AnnualSurvival %>%
  mutate(Condition=ifelse(feeding==0,"natural","food supplemented")) %>%
  mutate(Date=ymd("2010-05-01")+weeks(Fortnight*2)) %>% filter(Fortnight==27)

1-(finsurv[1,6]/finsurv[3,6])
1-(finsurv[2,6]/finsurv[4,6])
1-(finsurv[1,6]/finsurv[2,6])
1-(finsurv[3,6]/finsurv[4,6])