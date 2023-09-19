##########################################################################
#
# LITTLE OWL JUVENILE SURVIVAL ANALYSIS FROM TELEMETRY - PRIOR SENSIBILITY TEST
#
##########################################################################
# written by Steffen Oppel, September 2023
# motivated by paper Bannan et al 2020 about appropriate use of priors in ecology
# simulate survival probabilities from specified priors to ensure they cover plausible survival probabilities

### CONVERSION OF PRECISION SPECIFICATION
# in JAGS, precision in the normal distribution is specified by 1/variance
# in R, precision in normal distribution is specified by sqrt(variance)

precconv<-function(x){sqrt(1/x)}



### CREATE DATA FRAME WITH 10000 RANDOM VALUES DRAWN FROM PRIORS

mean.phi <- rbeta(10000,95, 10)   # uninformative prior for all biweekly survival probabilities
lp.mean <- log(mean.phi/(1 - mean.phi))    # logit transformed survival intercept
beta.yr <- rnorm(10000,0, precconv(1))         # Prior for annual variation on survival probability on logit scale
beta.win <- runif(10000,-2, 0)          # Prior for winter weather effect on survival probability on logit scale
beta.simpleage <- rnorm(10000,0, precconv(1))          # Prior for simple age offset on survival probability on logit scale
beta.male <- rnorm(10000,0, precconv(1))            # Prior for male effect on survival probability on logit scale
sigma.surv <- runif(10000,0, 2)                     # Prior for standard deviation of survival for random individual effect
tau.surv <- sigma.surv^-2
surv.raneff <- rnorm(10000,0, precconv(tau.surv))

INPUT<-data.frame(lp.mean,beta.yr,beta.win,beta.simpleage,beta.male,surv.raneff, simul=1)
FAKEDATA<-expand.grid(sex=c(0,1),age=c(0),win=scale(c(0,3,6,9))) %>% mutate(simul=1)


#### SPECIFY EQUATION TO CALCULATE SURVIVAL PROBABILITY FROM PRIORS
FAKEDATA %>% full_join(INPUT, by='simul',relationship = "many-to-many") %>%
  mutate(logit_phi=
           beta.yr +       ### survival varies by year
           beta.win*(win) +       ### survival dependent on winter weather (scaled)
           beta.male*(sex) +     ### survival dependent on sex
           beta.simpleage*(age) +     ### survival dependent on age on 1 Aug (scaled)
            surv.raneff +
           lp.mean) %>%        ### intercept for mean survival 
  mutate(phi=plogis(logit_phi)) %>%
  mutate(win=as.factor(win)) %>%
  group_by(age,sex,win) %>%
  
  ggplot() + geom_histogram(aes(x=phi, fill=win, col=win)) +
   facet_grid(win~sex)


