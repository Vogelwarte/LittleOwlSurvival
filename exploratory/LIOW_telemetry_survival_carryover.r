##########################################################################
#
# LITTLE OWL BIWEEKLY SURVIVAL ANALYSIS FROM TELEMETRY - CARRYOVER EFFECTS
#
##########################################################################
# written by Steffen Oppel, July 2023
# based on data by Marco Perrig

## additional analysis requested by Beat Naef in January 2024
## explore whether harsh winter also affects spring survival (carry-over effect)

library(runjags)
library(tidyverse)
library(data.table)
library(lubridate)
library(tidyverse)
filter<-dplyr::filter
select<-dplyr::select
library(MCMCvis)
library(foreach)
library(doParallel)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM PREPARED WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### data preparation moved to LIOW_telemetry_data_prep.r
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival"),silent=T)
# renv::init()   ### need to re-run this when you add a new library that needs to be run on the server
# renv::snapshot()

load("data/LIOW_SURV_INPUT.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY SURVIVAL MODEL IN JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# based on BPA book Kery and Schaub 2012, Chapter 7
# 7.5. Model with individual variation


# Specify model in JAGS language
sink("models/LIOW_CJS_carryover.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions)){
      logit(phi[i,t]) <- mu[season[t]] +
                        beta.mass*weight[i]*pf[t] +
                        beta.feed*feeding[i]*pf[t] + 
                        beta.win*env[year[i],t] +
                        beta.spring*winter[year[i]]*sp[t] +
                        beta.male*sex[i] #+
                        ##epsilon[i]    ##  beta.simpleage*simpleage[i] + beta.mass*weight[i] + beta.size*size[i] + 
      logit(p[i,t]) <- mu.p[recap.mat[i,t]] + beta.p.win*env[year[i],t] + epsilon.p[i]  ##  beta.p.yr[year[i]] + 
      } #t
   } #i
for (i in 1:nind){
   epsilon.p[i] ~ dnorm(0, tau.p)
}
   
  for (s in 1:4){   ### baseline for the 3 seasons dispersal, winter, breeding
   mean.phi[s] ~ dbeta(95, 10)                   # Prior for mean biweekly survival from Thorup et al. 2013, converted to beta
   mu[s] <- log(mean.phi[s] / (1-mean.phi[s]))       # Logit transformation
  }
   
   mean.p[1] ~ dunif(0.7, 1)                     # Prior for mean recapture during full effort periods
   mean.p[2] ~ dunif(0.3, 0.9)                  # Prior for mean recapture during reduced effort periods
   for (y in 1:2) {
    mu.p[y] <- log(mean.p[y] / (1-mean.p[y]))       # Logit transformation 
   }
  mu.p[3] <- -999999999999999999      # recapture probability of zero on logit scale 

sigma.p ~ dunif(0, 2)                      # Prior for standard deviation for random detection effect
tau.p <- pow(sigma.p, -2)

beta.mass ~ dnorm(0, 1)                     # Prior for mass effect
beta.male ~ dnorm(0, 1)                     # Prior for sex effect (for males, females are 0)
beta.win ~ dunif(-2, 2)                     # Prior for winter weather effect, which we know is negative
beta.spring ~ dunif(-2, 2)                     # Prior for spring carryover effect, which could be positive or negative
beta.p.win ~ dnorm(0, 1)                     # Prior for winter weather DETECTION effect
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
      y[i,t] ~ dbern(p[i,t] * z[i,t])
      y.rep[i,t] ~ dbern(p[i,t] * z.rep[i,t]) # replicate y (observations)

    
    } #t end
      #Derived parameters
      
        ## GOODNESS OF FIT TEST SECTION
        ## Discrepancy observed data
        E.obs[i] <- pow((sum(y[i,(f[i]+1):n.occasions]) - sum(p[i,(f[i]+1):(n.occasions)] * z[i,(f[i]+1):n.occasions])), 2) / (sum(p[i,(f[i]+1):n.occasions] * z[i,(f[i]+1):n.occasions]) + 0.001)

        ## Discrepancy replicated data
        E.rep[i] <- pow((sum(y.rep[i,(f[i]+1):n.occasions]) - sum(p[i,(f[i]+1):(n.occasions)] * z.rep[i,(f[i]+1):n.occasions])), 2) / (sum(p[i,(f[i]+1):(n.occasions)] * z.rep[i,(f[i]+1):n.occasions]) + 0.001)
      
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


### ENSURE REPEATABLE SCALING OF SNOWMAT ##
snowmean<-mean(as.matrix((allcov.new %>% dplyr::filter(variable=="day.snow.cover0"))[,c(3:32)]))
snowsd<-sd(as.matrix((allcov.new %>% dplyr::filter(variable=="day.snow.cover0"))[,c(3:32)]))
snowmat<-(as.matrix((allcov.new %>% dplyr::filter(variable=="day.snow.cover0"))[,c(3:32)])-snowmean)/snowsd

### QUANTIFY HARSHNESS OF WINTER 
winter<-allcov.new %>% dplyr::filter(variable=="day.snow.cover0") %>%
  select(-variable) %>%
  gather(key='period', value='snow',-ch.year) %>%
  group_by(ch.year) %>%
  summarise(N=sum(snow))

### ENSURE REPEATABLE SCALING OF SNOWMAT ##
INPUT <- list(y = CH, f = f,
              nind = dim(CH)[1],
              n.occasions = dim(CH)[2],
              z = known.state.cjs(CH),
              recap.mat=recap.mat,
              season=season,
              feeding=feeding,
              winter=scale(winter$N)[1:3,1],
              pf=ifelse(season==1,1,0), # to specify the post-fledging season and facilitate an age effect only for that season
              sp=ifelse(season==4,1,0), # to specify the spring season and facilitate a carryover effect only for that season
              sex=sex,
              year=as.numeric(year),
              weight=as.numeric(weight_scale),
              env=snowmat)  ### select any of the winter covariates 
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
                         mean.phi = rbeta(4, 94, 5),
                         mean.p = c(runif(1, 0.71, 0.98),runif(1, 0.31, 0.89)),
                         sigma.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mu","mean.phi", "mean.p", "beta.yr","beta.male","beta.win","beta.spring","beta.mass","beta.feed","beta.p.win","deviance","fit","fit.rep")

# MCMC settings
nt <- 6
nb <- 200
nc <- 3
nad<-100
ns<-1000
ni<-1500

# Call JAGS from R
carryover.model <- run.jags(data=INPUT, inits=inits, monitor=parameters,
                    model="C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival/models/LIOW_CJS_carryover.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 


#### MODEL ASSESSMENT ####
MCMCplot(carryover.model$mcmc, params=c("mean.phi","beta.win","beta.spring","beta.male","beta.mass","beta.feed","beta.p.win","mean.p"))
# ggsave("C:/Users/sop/OneDrive - Vogelwarte/General - Little owls/ANALYSES/LittleOwlSurvival/output/Carryover_parameter_estimates.jpg", height=11, width=8)
MCMCsummary(carryover.model$mcmc)


