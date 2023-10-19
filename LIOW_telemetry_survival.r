##########################################################################
#
# LITTLE OWL BIWEEKLY SURVIVAL ANALYSIS FROM TELEMETRY
#
##########################################################################
# written by Steffen Oppel, July 2023
# based on data by Marco Perrig

## added winter covariates on 10 Aug 2023
## completed winter variable selection on 11 Aug 2023
## added null model on 14 Aug 2023, moved variable selection section to end of script

## added GOF test and detection variables on 14 Aug 2023
## revised model and repeated variable selection

## added different recapture probabilities for certain occasions when no field team was available on 21 Aug 2023
## removed age and size effects as age was confounded with stage and led to very low survival during dispersal with massive uncertainty

## REVISED 15 SEPT 2023: distinguish only 3 stages (not 4) by combining incubation and chick rearing

## INCLUDED INFORMATIVE PRIORS FOR SURVIVAL
## annual survival: 0.258 (Le Gouar 2011: https://link.springer.com/article/10.1007/s00442-010-1868-x)
## monthly survival: 0.83 - 0.98 (Thorup 2010: https://link.springer.com/article/10.1007/s10336-012-0885-4/tables/2


## REVISED 19 SEPT 2023: need to include age as offset (single value) and sex (Tschumi et al. 2019)

## REVISED 27 SEPT 2023: trying to include post-fledging data and re-instate age (data from Perrig et al. 2017)
## added creation of R env file for outsourcing the model selection
## moved model selection into separate script to be run on server

## when including all data the temperature has the most outstanding effect, not snow cover,
## therefore need to build model with weather only affecting winter survival

## REVISED 19 OCT 2023 after exhaustive model selection decided on final model


library(runjags)
library(tidyverse)
library(data.table)
library(lubridate)
library(tidyverse)
# library(geosphere)
filter<-dplyr::filter
select<-dplyr::select
library(MCMCvis)
library(foreach)
library(doParallel)
# library(RMark)
# library(stringr)
# library(R2jags)
# library(renv)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM PREPARED WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### data preparation moved to LIOW_telemetry_data_prep.r
setwd("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
# renv::init()   ### need to re-run this when you add a new library that needs to be run on the server
# renv::snapshot()

load("data/LIOW_SURV_INPUT.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY SURVIVAL MODEL IN JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# based on BPA book Kery and Schaub 2012, Chapter 7
# 7.5. Model with individual variation


# Specify model in JAGS language
sink("models/LIOW_CJS_FINAL.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu[season[t]] +
                        #beta.yr[year[i]] +
                        beta.mass*weight[i] +
                        beta.feed*feeding[i] + 
                        #beta.age*age[i,t] +   ## structure age so as to be only used for post-fledging phase
                        beta.win*env[year[i],t]*winter[t] +
                        beta.male*sex[i] +
                        epsilon[i]    ##  beta.simpleage*simpleage[i] + beta.mass*weight[i] + beta.size*size[i] + 
      logit(p[i,t]) <- mu.p[recap.mat[i,t]] + beta.p.win*env[year[i],t] + epsilon.p[i]  ##  beta.p.yr[year[i]] + 
      } #t
   } #i
for (i in 1:nind){
   epsilon[i] ~ dnorm(0, tau)
   epsilon.p[i] ~ dnorm(0, tau.p)
}
   
  for (s in 1:4){   ### baseline for the 3 seasons dispersal, winter, breeding
   mean.phi[s] ~ dbeta(95, 10)                   # Prior for mean biweekly survival from Thorup et al. 2013, converted to beta
   mu[s] <- log(mean.phi[s] / (1-mean.phi[s]))       # Logit transformation
  }
   
   mean.p[1] ~ dunif(0.9, 1)                     # Prior for mean recapture during full effort periods
   mean.p[2] ~ dunif(0.3, 0.9)                  # Prior for mean recapture during reduced effort periods
   for (y in 1:2) {
    mu.p[y] <- log(mean.p[y] / (1-mean.p[y]))       # Logit transformation 
   }
  mu.p[3] <- -999999999999999999      # recapture probability of zero on logit scale 

sigma ~ dunif(0, 2)                      # Prior for standard deviation for random survival effect
tau <- pow(sigma, -2)
sigma.p ~ dunif(0, 2)                      # Prior for standard deviation for random detection effect
tau.p <- pow(sigma.p, -2)

# for (y in 1:3) {
#  beta.yr[y] ~ dnorm(0, 1)                     # Prior for year effect 
#  #beta.p.yr[y] ~ dnorm(0, 1)                 # Prior for ANNUAL DETECTION effect
# }

#beta.size ~ dnorm(0, 1)                     # Prior for size effect 
#beta.age ~ dnorm(0, 1)                     # Prior for age effect 
beta.mass ~ dnorm(0, 1)                     # Prior for mass effect
#beta.simpleage ~ dnorm(0, 1)                # Prior for age offset (simple value for each bird according to age at 1 Aug) 
beta.male ~ dnorm(0, 1)                     # Prior for sex effect (for males, females are 0)
beta.win ~ dunif(-2, 0)                     # Prior for winter weather effect, which we know is negative
#beta.p.win ~ dnorm(0, 1)                     # Prior for winter weather DETECTION effect
beta.feed ~ dnorm(0, 1)                # Prior for effect of supplementary feeding


# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   z.rep[i,f[i]] <- 1 # replicate z (true state)
   y.rep[i,f[i]] <- 1 # replicate y (data)
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
      z.rep[i,t] ~ dbern(phi[i,t-1] * z.rep[i,t-1]) # replicate z (true state)
      # Observation process
      y[i,t] ~ dbern(p[i,t-1] * z[i,t])
      y.rep[i,t] ~ dbern(p[i,t-1] * z.rep[i,t]) # replicate y (observations)

    
    } #t end
      #Derived parameters
      
        ## GOODNESS OF FIT TEST SECTION
        ## Discrepancy observed data
        E.obs[i] <- pow((sum(y[i,(f[i]+1):n.occasions]) - sum(p[i,f[i]:(n.occasions-1)] * z[i,(f[i]+1):n.occasions])), 2) / (sum(p[i,f[i]:(n.occasions-1)] * z[i,(f[i]+1):n.occasions]) + 0.001)

        ## Discrepancy replicated data
        E.rep[i] <- pow((sum(y.rep[i,(f[i]+1):n.occasions]) - sum(p[i,f[i]:(n.occasions-1)] * z.rep[i,(f[i]+1):n.occasions])), 2) / (sum(p[i,f[i]:(n.occasions-1)] * z.rep[i,(f[i]+1):n.occasions]) + 0.001)
      
   } #i end
      fit <- sum(E.obs[])
      fit.rep <- sum(E.rep[])
}
",fill = TRUE)
sink()

# Bundle data
# Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

INPUT <- list(y = CH, f = f,
              nind = dim(CH)[1],
              n.occasions = dim(CH)[2],
              z = known.state.cjs(CH),
              recap.mat=recap.mat,
              season=season,
              feeding=feeding,
              #winter=ifelse(season==3,1,0),
              #age=age_scale,
              #pf=ifelse(season==1,1,0), # to specify the post-fledging season and facilitate an age effect only for that season
              #simpleage=as.numeric(simpleage_scale),
              sex=sex,
              #size=size,
              year=as.numeric(year),
              weight=weight,
              env=as.matrix((allcov.new %>% dplyr::filter(variable=="day.snow.cover0"))[,c(3:32)]))  ### select any of the winter covariates 
              #env=as.matrix((allcov %>% dplyr::filter(variable=="day.snow.cover5"))[,c(26:31,3:25)]))  ### select any of the winter covariates 
              #rain=as.matrix((allcov %>% dplyr::filter(variable=="total.precip"))[,3:25]))  ### select any of the winter covariates 

# Initial values 
# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

inits <- function(){list(z = cjs.init.z(CH, f),
                         mean.phi = rbeta(4, 95, 10),
                         mean.p = c(runif(1, 0.9, 1),runif(1, 0.3, 0.9)),
                         sigma = runif(1, 0, 2))}  

# Parameters monitored
parameters <- c("mu","mean.phi", "mean.p", "beta.yr","beta.male","beta.win","beta.mass","beta.feed","beta.p.win","deviance","fit","fit.rep")

# MCMC settings
nt <- 6
nb <- 200
nc <- 3
nad<-100
ns<-2000
ni=3500

# Call JAGS from R
full.model <- run.jags(data=INPUT, inits=inits, monitor=parameters,
                    model="C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival/models/LIOW_CJS_FINAL.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 

## fitted in R2jags to retrieve sims.list for GoF test - later abandoned
# full.model <- R2jags::jags(data=INPUT, inits=inits, parameters.to.save=parameters,
#                        model.file="C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival/LIOW_CJS_model_GoF.jags",
#                    n.iter=ni, n.chains = nc, n.thin = nt, n.burnin = nb, DIC=T) 

parameters <- c("mu","mean.phi", "mean.p", "beta.male","beta.mass","beta.feed","beta.p.win","deviance","fit","fit.rep")
null.model <- run.jags(data=INPUT, inits=inits, monitor=parameters,
                    model="C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival/models/LIOW_CJS_FINAL_null.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 


#### MODEL COMPARISON ####
## needs more thought from: https://kevintshoemaker.github.io/NRES-746/LECTURE8.html#Explicit_Bayesian_model_selection
## THIS TEST YIELDS A WEIRD OUTPUT AND NO SIGNIFICANT DIFFERENCE
# dic.full<-extract.runjags(full.model,"DIC")
# dic.null<-extract.runjags(null.model,"DIC")
# 
# str(dic.full)
# hist(dic.full$deviance+dic.full$penalty)
# hist(dic.null$deviance+dic.null$penalty)
# t.test(x=dic.full$deviance+dic.full$penalty,y=dic.null$deviance+dic.null$penalty)

## manual calculation SHOWS CLEAR DIFFERENCE IN deviance
## full model has lower deviance than null model (including penalty of additional parameter) - hence is supported
full.model$summary$quantiles[16,c(3,1,5)] +2 <
null.model$summary$quantiles[15,c(3,1,5)]



#### MODEL ASSESSMENT ####
MCMCplot(full.model$mcmc, params=c("mean.phi","beta.win","beta.male","beta.mass","beta.feed","beta.p.win","mean.p"))
MCMCplot(null.model$mcmc, params=c("mean.phi","beta.male","beta.mass","beta.feed","beta.p.win","mean.p"))
MCMCsummary(full.model$mcmc)
MCMCsummary(null.model$mcmc)
MCMCdiag(full.model$mcmc,
         round = 3,
         file_name = 'LIOW_survival',
         dir = 'C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival/output',
         mkdir = 'LIOW_v6',
         add_field = '6.0',
         add_field_names = 'Data version',
         save_obj = TRUE,
         obj_name = 'LIOW-fit-19Oct2023',
         add_obj = list(INPUT, sessionInfo()),
         add_obj_names = c('surv-data-19Oct2023', 'session-info-19Oct2023'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOODNESS OF FIT ASSESSMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## https://agabrioblog.onrender.com/tutorial/gof-tests-jags/gof-tests-jags/
## https://www.sciencedirect.com/topics/earth-and-planetary-sciences/goodness-of-fit


# Evaluation of fit using R2ucare
# this may not be valid due to individual covariates: https://oliviergimenez.github.io/R2ucare/
# 
# library(R2ucare)
# # Replace the NAs in the simulated data sets by 0
# full.model$BUGSoutput$sims.list$y.rep[is.na(full.model$mcmc$y.rep)] <- 0
# # Perform U-CARE GoF for each replicated data set
# chi.3sr <- chi.3sm <- chi.2cl <- chi.2ct <- numeric()
# for (i in 1:dim(full.model$BUGSoutput$sims.list$y.rep)[1]){
#   chi.3sr[i] <- test3sr(full.model$BUGSoutput$sims.list$y.rep[i,,], rep(1, nrow(CH)))$test3sr[1]
#   chi.3sm[i] <- test3sm(full.model$BUGSoutput$sims.list$y.rep[i,,], rep(1, nrow(CH)))$test3sm[1]
#   chi.2ct[i] <- test2ct(full.model$BUGSoutput$sims.list$y.rep[i,,], rep(1, nrow(CH)))$test2ct[1]
#   chi.2cl[i] <- test2cl(full.model$BUGSoutput$sims.list$y.rep[i,,], rep(1, nrow(CH)))$test2cl[1]
# }
# 
# # Perform U-CARE GoF for each the real (actual) data set
# chi.3sr.actual <- test3sr(CH, rep(1, nrow(CH)))$test3sr[1]
# chi.3sm.actual <- test3sm(CH, rep(1, nrow(CH)))$test3sm[1]
# chi.2ct.actual <- test2ct(CH, rep(1, nrow(CH)))$test2ct[1]
# chi.2cl.actual <- test2cl(CH, rep(1, nrow(CH)))$test2cl[1]
# 
# # Produce plots
# par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
# hist(chi.3sr, nclass=25, col="lightblue", main="3.SR", xlab=NA, xlim=c(0, 50))
# abline(v=chi.3sr.actual, col="red", lwd=2)
# p.3sr <- round(mean((chi.3sr.actual - chi.3sr) < 0, na.rm=TRUE), 3)
# legend("topright", legend=paste("P = ", p.3sr), bty="n")
# hist(chi.3sm, nclass=25, col="lightblue", main="3.SM", xlab=NA)
# abline(v=chi.3sm.actual, col="red", lwd=2)
# p.3sm <- round(mean((chi.3sm.actual - chi.3sm) < 0, na.rm=TRUE), 3)
# legend("topright", legend=paste("P = ", p.3sm), bty="n")
# hist(chi.2ct, nclass=25, col="lightblue", main="2.CT", xlab=NA)
# abline(v=chi.2ct.actual, col="red", lwd=2)
# p.2ct <- round(mean((chi.2ct.actual - chi.2ct) < 0, na.rm=TRUE), 3)
# legend("topright", legend=paste("P = ", p.2ct), bty="n")
# hist(chi.2cl, nclass=25, col="lightblue", main="2.CL", xlab=NA)
# abline(v=chi.2cl.actual, col="red", lwd=2)
# p.2cl <- round(mean((chi.2cl.actual - chi.2cl) < 0, na.rm=TRUE), 3)
# legend("topright", legend=paste("P = ", p.2cl), bty="n")

## USING THE CALCULATED FIT VALUES FROM THE JAGS MODEL
OBS <- MCMCpstr(full.model$mcmc, params=c("fit"), type="chains")
REP <- MCMCpstr(full.model$mcmc, params=c("fit.rep"), type="chains")
GOF<-tibble(Rep=as.numeric(REP[[1]]),Obs=as.numeric(OBS[[1]])) %>%
  mutate(P=ifelse(Obs>Rep,1,0))

ggplot(GOF,aes(x=Rep,y=Obs, fill=P)) + geom_point(position=position_jitterdodge()) +
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position="none") +
  annotate("text",label=as.character(round(mean(GOF$P),2)),x=20,y=20,size=8, colour="firebrick")

mean(GOF$P)


## CHECK WHETHER NULL MODEL FITS THE DATA
OBS.null <- MCMCpstr(null.model$mcmc, params=c("fit"), type="chains")
REP.null <- MCMCpstr(null.model$mcmc, params=c("fit.rep"), type="chains")
GOF.null<-tibble(Rep=as.numeric(REP.null[[1]]),Obs=as.numeric(OBS.null[[1]])) %>%
  mutate(P=ifelse(Obs>Rep,1,0))

ggplot(GOF.null,aes(x=Rep,y=Obs, fill=P)) + geom_point(position=position_jitterdodge()) +
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position="none") +
  annotate("text",label=as.character(round(mean(GOF.null$P),2)),x=20,y=20)

mean(GOF.null$P)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT SURVIVAL ESTIMATES FOR THE VARIOUS PHASES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### PREPARE RAW MCMC OUTPUT
parmcols<-dimnames(full.model$mcmc[[1]])[[2]]


# ### COMBINE SAMPLES ACROSS CHAINS
MCMCout<-rbind(full.model$mcmc[[1]],full.model$mcmc[[2]],full.model$mcmc[[3]])
# str(MCMCout)

### SET UP TABLE FOR PLOTTING THE SEASONAL SURVIVAL GRAPH

AnnTab<-data.frame(season=c(1,2,3,3,3,3,4),
                   age=c(45,98,180,190,200,210,300),
                   feeding=0,
                   weight=0,
                   sex=1,
                   snow=scale(c(0,0,0,4,8,12,0))[,1])  %>% 
  #mutate(pf=ifelse(season==1,1,0)) %>%
  mutate(scaleage=(age-attr(age_scale, 'scaled:scale')[10])/attr(age_scale, 'scaled:scale')[10]) 

Xin<-AnnTab

### CALCULATE PREDICTED VALUE FOR EACH SAMPLE

MCMCpred<-data.frame()
for(s in 1:nrow(MCMCout)) {
  
  X<-  Xin %>%
    
    ##CALCULATE MONTHLY SURVIVAL
    mutate(logit.surv=as.numeric(MCMCout[s,grepl("mu",parmcols)])[season]+
             as.numeric(MCMCout[s,match("beta.mass",parmcols)])*weight +
             #as.numeric(MCMCout[s,match("beta.yr[3]",parmcols)])+   #*year + ### categorical year effect - pick the most average year
             as.numeric(MCMCout[s,match("beta.male",parmcols)])*sex +
             as.numeric(MCMCout[s,match("beta.feed",parmcols)])*feeding +
             as.numeric(MCMCout[s,match("beta.win",parmcols)])*snow) %>%
    
    ## BACKTRANSFORM TO NORMAL SCALE
    mutate(surv=plogis(logit.surv)) %>%
    
    ## RENAME THE SEASONS
    mutate(Season=ifelse(season==2,"Autumn",
                         ifelse(season==3,"Winter",
                                ifelse(season==4,"Spring","Summer")))) %>%
    mutate(simul=s)              
  
  MCMCpred<-rbind(MCMCpred,as.data.frame(X)) 
  
}

### CREATE PLOT

plotdat<-  MCMCpred %>% rename(raw.surv=surv) %>%
  mutate(age=rep(c(45,98,180,190,200,210,300), ns*nc)) %>%
  group_by(Season,age,snow) %>%
  summarise(surv=quantile(raw.surv,0.5),surv.lcl=quantile(raw.surv,0.025),surv.ucl=quantile(raw.surv,0.975)) %>%
  ungroup() %>%
  mutate(snow=c(0,0,0,0,4,8,12)) %>%
  arrange(age)

 
ggplot(plotdat)+
  geom_errorbar(aes(x=age, ymin=surv.lcl, ymax=surv.ucl, colour=factor(snow)), width=0.2) +   ##, type=Origin
  geom_point(aes(x=age, y=surv,colour=factor(snow)),size=2)+     ## , linetype=Origin
  
  ## format axis ticks
  #scale_x_continuous(name="Season", limits=c(1,365), breaks=plotdat$age[c(3,5,8,10)], labels=plotdat$Season[c(3,5,8,10)]) +
  scale_x_continuous(name="Season", limits=c(1,365), breaks=plotdat$age[c(1,2,4,7)], labels=plotdat$Season[c(1,2,4,7)]) +
  scale_y_continuous(name="Biweekly survival probability", limits=c(0.7,1), breaks=seq(0.7,1,0.05), labels=seq(0.7,1,0.05)) +
  #scale_y_continuous(name="Monthly survival probability", limits=c(0.8,1), breaks=seq(0.,1,0.05)) +
  labs(y="Biweekly survival probability") +
  scale_colour_manual(name="Days of ≥ 1 cm\nsnow cover", values=c("black", "goldenrod", "darkorange", "firebrick"),
                      breaks=c(0,4,8,12),labels=c(0,4,8,12)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=16, color="black"),
        legend.background=element_blank(),
        legend.key = element_rect(fill = NA),
        legend.position=c(0.12,0.18), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))


ggsave("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival/output/Seasonal_survival_LIOW.jpg", height=7, width=11)
ggsave("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/LittleOwlSurvival/Fig_1.jpg", height=7, width=11)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CUMULATIVE ANNUAL SURVIVAL PREDICTION FOR FIRST YEAR LITTLE OWLS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stage.surv<-  plotdat %>%
  mutate(dur=c(5,6,10,10,10,10,5)) %>%
  mutate(surv=surv^dur,surv.lcl=surv.lcl^dur,surv.ucl=surv.ucl^dur)
stage.surv


Table1<- stage.surv[c(1:3,7),] %>%
  mutate(mild.survival=sprintf("%s (%s - %s)",round(surv,3),round(surv.lcl,3),round(surv.ucl,3))) %>%
  mutate(Duration=dur*2) %>%
  select(Season,Duration,mild.survival) %>%
  bind_rows(data.frame(Season="Annual",Duration=52,
                       mild.survival=sprintf("%s (%s - %s)",
                                        round(prod(stage.surv[c(1:3,7),4]),3),
                                        round(prod(stage.surv[c(1:3,7),5]),3),
                                        round(prod(stage.surv[c(1:3,7),6]),3))))

Table1<- stage.surv[c(1:2,6:7),] %>%
  mutate(harsh.survival=sprintf("%s (%s - %s)",round(surv,3),round(surv.lcl,3),round(surv.ucl,3))) %>%
  select(Season,harsh.survival) %>%
  bind_rows(data.frame(Season="Annual",
                       harsh.survival=sprintf("%s (%s - %s)",
                                        round(prod(stage.surv[c(1:2,6:7),4]),3),
                                        round(prod(stage.surv[c(1:2,6:7),5]),3),
                                        round(prod(stage.surv[c(1:2,6:7),6]),3))))  %>%
  left_join(Table1, by="Season") %>%
  select(Season,Duration,mild.survival,harsh.survival)

Table1
fwrite(Table1,"C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/LittleOwlSurvival/Table1_surv.csv")


  
# save.image("LIOW_survival_output.RData")
# load("LIOW_survival_output.RData")


### CALCULATE REDUCTION IN % ###

(0.125-0.102)/0.125



##### CREATE TABLE S1 FROM MODEL SELECTION RESULTS  ##########

winter.vars<-fread("output/LIOW_model_selection_DIC_table_winenv.csv")%>%
  arrange(DIC)

winter.vars


## create TABLE S1

TABLES1 <- winter.vars %>%
  mutate(winter.effect=sprintf("%s (%s – %s)",round(beta_med,3),round(beta_lcl,3),round(beta_ucl,3))) %>%
  mutate(size.effect=sprintf("%s (%s – %s)",round(beta_size,3),round(beta_size_lcl,3),round(beta_size_ucl,3))) %>%
  mutate(age.effect=sprintf("%s (%s – %s)",round(beta_age_med,3),round(beta_age_lcl,3),round(beta_age_ucl,3))) %>%
  mutate(feed.effect=sprintf("%s (%s – %s)",round(beta_feed,3),round(beta_feed_lcl,3),round(beta_feed_ucl,3))) %>%
  mutate(age.effect=if_else(age=="no","",age.effect)) %>%
  mutate(feed.effect=if_else(feeding=="no","",age.effect)) %>%
  mutate(size.effect=if_else(size=="none","",size.effect)) %>%
  arrange(DIC) %>%
  select(var,size,age,feeding,DIC,winter.effect,size.effect,age.effect,feed.effect) %>%
  rename(Winter.variable=var)
TABLES1
fwrite(TABLES1,"C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/LittleOwlSurvival/TableS1_DIC.csv")







#### MODEL SELECTION VIA DIC ####
## follow post by Bob O'Hara
## http://deepthoughtsandsilliness.blogspot.com/2007/12/focus-on-dic.html


### DID NOT MANAGE TO IMPLEMENT THIS !!!
# #######################################################################
# DevCalc1=function(mcmc, data) {
# NamesGrp=paste("muGrp[", data$Group, "]", sep="")
# -2*sum(dnorm(data$Y, mcmc[NamesGrp], mcmc["sigma.y"], log=TRUE))
# }
# DevCalc2=function(mcmc, data) {
#   Mean=unlist(tapply(data$y, list(data$Group), mean))
#   Mu=mcmc["mu0"]
#   N=length(data$Group)/length(unique(data$Group))
#   -2*sum(dnorm(Mean, Mu, sqrt((mcmc["sigma.y"]^2)/N + mcmc["sigma.Grp"]^2), log=TRUE)) }
# DevCalc2beta=function(mcmc, data) {
#   Mean=unlist(tapply(data$y, 1,mean))
#   Mu=MCMCout[s,grepl("mu",parmcols)] + mcmc["betaGrp"]*(as.numeric(names(Mean)) - 5.5)
#   N=length(data$Group)/length(unique(data$Group))
#   -2*sum(dnorm(Mean, Mu, sqrt((mcmc["sigma.y"]^2)/N + mcmc["sigma.Grp"]^2), log=TRUE))
# }
# DIC.calc=function(dat, mcmc, func) {
#   Dev=apply(mcmc$sims.array, c(1,2), func, data=dat)
#   mean=apply(mcmc$sims.array, 3, mean)
#   pD=mean(Dev) - func(mean, dat)
#   DIC=mean(Dev) + pD
#   return(list(DIC=DIC, pD=pD, Dbar=mean(Dev)))
# }
# PrintDIC=function(DIC, Label) { cat(Label, " DIC:", DIC$DIC, " pD:", DIC$pD, "\n")}
# 
# 
# # Focus on phi
# DIC2.D1M1=DIC.calc(INPUT, modelfit$mcmc, DevCalc2beta)
# PrintDIC(DIC2.D1M1, "Data 1, Model 1")
