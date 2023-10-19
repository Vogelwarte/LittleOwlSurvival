##########################################################################
#
# LITTLE OWL BIWEEKLY SURVIVAL ANALYSIS FROM TELEMETRY - MODEL SELECTION
#
##########################################################################
# written by Steffen Oppel, July 2023
# based on data by Marco Perrig
# this script only performs model selection and requires input data to be generated from "LIOW_telemetry_survival.r"

# updated on 17 OCT 2023 to include parallel processing to work through >100 models in less time
# parallel computation taken from: https://stackoverflow.com/questions/53597282/running-multiple-parallel-processes-in-parallel-r

# updated on 18 OCT 2023 after chat with Murdy - snow is biologically most sensible variable
# skip model selection of lots of weather variables, just take the snow and use that

library(runjags)
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
filter<-dplyr::filter
select<-dplyr::select


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM PREPARED ENVIRONMENT FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/General/ANALYSES/LittleOwlSurvival")
load("data/LIOW_SURV_INPUT.RData")

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

# MCMC settings
nt <- 6
nb <- 200
nc <- 3
nad<-10
ns<-2000
ni<-3500

inits <- function(){list(z = cjs.init.z(CH, f),
                         mean.phi = rbeta(4, 95, 10),
                         mean.p = c(runif(1, 0.9, 1),runif(1, 0.3, 0.9)),
                         sigma = runif(1, 0, 2))} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPLORE BEST WINTER SURVIVAL PREDICTOR AND INCLUSION OF AGE AND SIZE - fully revised to final model on 20 Sept
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

winter.vars<-expand.grid(var=unique(allcov.new$variable)[c(2,6)],
                         feeding=c("yes","no"),
                         size=c("size","mass","none"),
                         age=c("yes","no"),
                         dic_med=0, dic_lcl=0,dic_ucl=0,beta_med=0, beta_lcl=0,beta_ucl=0,
                         beta_size=0, beta_size_lcl=0,beta_size_ucl=0,
                         beta_age_med=0, beta_age_lcl=0,beta_age_ucl=0,
                         beta_feed=0, beta_feed_lcl=0,beta_feed_ucl=0,
                         beta_winP=0, beta_winP_lcl=0,beta_winP_ucl=0)
param2 <- c("deviance","beta.win","beta.size","beta.mass","beta.age","beta.feed","beta.p.win")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# START PARALLEL PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n.cores <- 24  ## this is to be run on the server
registerDoParallel(n.cores)

#Fit models in parallel, with chains running in parallel.----


mod.sel.results <- 
  foreach(s = 1:dim(winter.vars)[1],.combine=rbind, .packages='runjags',.inorder=FALSE,.errorhandling="remove",.verbose=FALSE) %dopar% {
    
    
  ##for(s in 1:dim(winter.vars)[1]){   ## the old sequential loop
  output<-winter.vars[s,]
  INPUT <- list(y = CH, f = f,
                nind = dim(CH)[1],
                n.occasions = dim(CH)[2],
                z = known.state.cjs(CH),
                recap.mat=recap.mat,
                season=season,
		    winter=ifelse(season==3,1,0),
                feeding=feeding,
		    age=age_scale,
                sex=sex,
                size=size,
                weight=weight,
                year=as.numeric(year),
		            env=as.matrix((allcov.new %>% dplyr::filter(variable==output$var))[,c(3:32)]))  ### select any of the winter covariates
                #env=as.matrix((allcov %>% dplyr::filter(variable==winter.vars$var[s]))[,3:25]))  ### select any of the winter covariates
         	    #env=as.matrix((allcov %>% dplyr::filter(variable==winter.vars$var[s]))[,c(26:32,3:25)]))  ### select any of the winter covariates  

# Call JAGS from R
  if(output$size=="none" & output$age== "no" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }

  if(output$size=="size" & output$age== "no" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_noage_size.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="mass" & output$age== "no" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_noage_mass.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="size" & output$age== "yes" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_size.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="mass" & output$age== "yes" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_mass.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="none" & output$age== "yes" & output$feeding== "no"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }

### add models with feeding ###
  if(output$size=="none" & output$age== "no" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed_noage.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }

  if(output$size=="size" & output$age== "no" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed_noage_size.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="mass" & output$age== "no" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed_noage_mass.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="size" & output$age== "yes" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed_size.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="mass" & output$age== "yes" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed_mass.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }
  
  if(output$size=="none" & output$age== "yes" & output$feeding== "yes"){
    modelfit <- run.jags(data=INPUT, inits=inits, monitor=param2,
                         model="models/LIOW_CJS_fullyear_feed.jags",
                         n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                         method="parallel")
  }




##### WRITING OUTPUT OF PARAMETERS IN THE MODEL #########output/


output[,5:7]<-modelfit$summary$quantiles[1,c(3,1,5)]
output[,8:10]<-modelfit$summary$quantiles[2,c(3,1,5)]
output[,20:22]<-modelfit$summary$quantiles[dim(modelfit$summary$quantiles)[1],c(3,1,5)]

if(output$size %in% c("mass","size")){
  output[,11:13]<-modelfit$summary$quantiles[3,c(3,1,5)]
}
if(output$age== "yes"){
	if(output$size %in% c("mass","size")){	
  	output[,14:16]<-modelfit$summary$quantiles[4,c(3,1,5)]
	}else{output[,14:16]<-modelfit$summary$quantiles[3,c(3,1,5)]}
}

if(output$feeding== "yes"){
 if(output$age== "yes"){
	if(output$size %in% c("mass","size")){	
  	output[,17:19]<-modelfit$summary$quantiles[5,c(3,1,5)]
	}else{output[,17:19]<-modelfit$summary$quantiles[4,c(3,1,5)]}
 }else{
	if(output$size %in% c("mass","size")){	
  	output[s,17:19]<-modelfit$summary$quantiles[4,c(3,1,5)]
	}else{output[s,17:19]<-modelfit$summary$quantiles[3,c(3,1,5)]}
}}

#write.csv(output,"LIOW_win_var_selection_DIC_table.csv", append=TRUE)

return(output)
#save.image("output/LIOW_survival_output.RData")
} ### end loop over all models


############################-----------------------------------------------------------------------################################################################

mod.sel.results<-mod.sel.results %>% 
  mutate(DIC=if_else(age=="yes",
                     if_else(size=="none",
                             if_else(feeding=="yes",dic_med+4,dic_med+2),
                             if_else(feeding=="yes",dic_med+6,dic_med+4)),
                     if_else(size=="none",
                             if_else(feeding=="yes",dic_med+2,dic_med),
                             if_else(feeding=="yes",dic_med+4,dic_med+2)))) %>%
  arrange(abs(beta_med))

mod.sel.results
fwrite(mod.sel.results,"output/LIOW_model_selection_DIC_table_winenv.csv")

