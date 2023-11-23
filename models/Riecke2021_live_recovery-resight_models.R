############################################################################################
# This script simulates capture-mark-recapture-resight-recovery data and subsequently analyzes
# different subsets of those data using models from Barker (1997), Lindberg et al. (2001),
# and Kendall et al. (2013)
############################################################################################
#
#
#
#
# lines 031-176: Barker (1997) - Barker1997
# lines 188-346: Lindberg et al. (2001) - Lindberg2001
# lines 354-512: Kendall et al. (2013) - Kendall2013
# lines 516-542: specify number of simulations, occasions, releases, and number of iterations/chains
# lines 546-573: begin simulation, simulate and store true parameter values
# lines 578-599: specify and simulate 'true' latent states
# lines 602-645: simulate observation process
# lines 651-677: summarize robust design detection data following Riecke et al. 2018
# lines 680-700: remove individuals that don't contribute to the likelihood (i.e., released in last occasion)
# lines 709-738: create Barker input and initial values
# lines 748-780: create Lindberg input and initial values
# lines 787-816: create Kendall input and initial values
# lines 820-840: JAGS data/inits set-up
# lines 844-900: run models and store output
# lines 901-end: create tables and figures...
#
#
#
#
############################################################################################

############################################################################################
# Specify Barker (1997) model in BUGS language
############################################################################################
sink("Barker1997.jags")
cat("
    model {
    
    ####################################################################################  
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # F: site fidelity
    # p: capture probability
    # R: non-breeding resight probability given survival from t to t+1
    # Rp (R'): non-breeding resight probability given mortality during interval t to t+1
    # f: band recovery probability
    #--------------------------------------
    ####################################################################################  
    
    phi ~ dbeta(1,1)
    F ~ dbeta(1,1)
    r ~ dbeta(1,1)
    R ~ dbeta(1,1)
    Rp ~ dbeta(1,1)
    gammapp <- 1
    gammap <- gammapp
    pstar ~ dbeta(1,1)
    
    
    #####################################################################################
    #--------------------------------------
    # Latent states (Psi):
    # 1) alive, available
    # 2) permanent emigraht
    # 3) recently dead, shot
    # 4) recently dead, resighted
    # 5) dead in previous interval
    #--------------------------------------
    #####################################################################################
    
    Psi[1,1] <- phi * F * gammapp
    Psi[1,2] <- phi * F * (1 - gammapp)
    Psi[1,3] <- phi * (1-F)
    Psi[1,4] <- (1 - phi)*r
    Psi[1,5] <- (1 - phi)*Rp*(1-r)
    Psi[1,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[2,1] <- phi * F * gammap
    Psi[2,2] <- phi * F * (1 - gammap)
    Psi[2,3] <- phi * (1-F)
    Psi[2,4] <- (1 - phi)*r
    Psi[2,5] <- (1 - phi)*Rp*(1-r)
    Psi[2,6] <- (1 - phi)*(1 - Rp)*(1-r)

    Psi[3,1] <- 0
    Psi[3,2] <- 0
    Psi[3,3] <- phi
    Psi[3,4] <- (1 - phi)*r
    Psi[3,5] <- (1 - phi)*Rp*(1-r)
    Psi[3,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[4,1] <- 0
    Psi[4,2] <- 0
    Psi[4,3] <- 0
    Psi[4,4] <- 0
    Psi[4,5] <- 0
    Psi[4,6] <- 1
    
    Psi[5,1] <- 0
    Psi[5,2] <- 0
    Psi[5,3] <- 0
    Psi[5,4] <- 0
    Psi[5,5] <- 0
    Psi[5,6] <- 1

    Psi[6,1] <- 0
    Psi[6,2] <- 0
    Psi[6,3] <- 0
    Psi[6,4] <- 0
    Psi[6,5] <- 0
    Psi[6,6] <- 1    
    
    #####################################################################################
    #--------------------------------------
    # Breeding observations (o ~ Omega):
    # 1) seen breeding, seen non-breeding
    # 2) seen breeding, not seen non-breeding
    # 3) not seen breeding, seen non-breeding
    # 4) not seen breeding, not seen non-breeding
    # 5) shot
    #--------------------------------------
    #####################################################################################
    
    Omega[1,1] <- pstar*R                       
    Omega[1,2] <- pstar*(1 - R)  
    Omega[1,3] <- (1 - pstar)*R 
    Omega[1,4] <- (1 - pstar)*(1 - R)
    Omega[1,5] <- 0 
    
    Omega[2,1] <- 0                       
    Omega[2,2] <- 0   
    Omega[2,3] <- R 
    Omega[2,4] <- (1 - R) 
    Omega[2,5] <- 0 

    Omega[3,1] <- 0                       
    Omega[3,2] <- 0   
    Omega[3,3] <- R 
    Omega[3,4] <- (1 - R) 
    Omega[3,5] <- 0 
    
    Omega[4,1] <- 0    
    Omega[4,2] <- 0 
    Omega[4,3] <- 0 
    Omega[4,4] <- 0 
    Omega[4,5] <- 1
    
    Omega[5,1] <- 0    
    Omega[5,2] <- 0 
    Omega[5,3] <- 1 
    Omega[5,4] <- 0 
    Omega[5,5] <- 0 

    Omega[6,1] <- 0    
    Omega[6,2] <- 0 
    Omega[6,3] <- 0 
    Omega[6,4] <- 1 
    Omega[6,5] <- 0 
    
    
    # Likelihood 
    for (i in 1:n.ind){
    
    z[i, first[i]] <- y[i,first[i]]
    
    for (t in (first[i] + 1):n.occasions){
    
    z[i,t] ~ dcat(Psi[z[i,t-1], ])
    
    y[i,t] ~ dcat(Omega[z[i,t], ])
    
    }}}
    
    
    ",fill = TRUE)
sink()











############################################################################################
# Specify Lindberg et al. (2001) model in BUGS language
############################################################################################
sink("Lindberg2001.jags")
cat("
    model {
    
    ####################################################################################  
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # F: site fidelity
    # p: capture probability
    # R: non-breeding resight probability given survival from t to t+1
    # Rp (R'): non-breeding resight probability given mortality during interval t to t+1
    # f: band recovery probability
    #--------------------------------------
    ####################################################################################  
    
    phi ~ dbeta(1,1)
    F ~ dbeta(1,1)
    r ~ dbeta(1,1)
    R <- 0
    Rp <- 0
    gammapp ~ dbeta(1,1)
    gammap <- gammapp

    for (j in 1:n.sec){
      p[j] ~ dbeta(1,1)
    }

    #####################################################################################
    # model secondary detection probabilities
    #####################################################################################
    for (t in 1:n.occasions){
      for (j in 1:n.sec){
        yes[t,j] ~ dbin(p[j], total[t,j])
      }
    }
    pstar <- 1 - prod(1 - p[])
    
    
    #####################################################################################
    #--------------------------------------
    # Latent states (Psi):
    # 1) alive, available
    # 2) permanent emigraht
    # 3) recently dead, shot
    # 4) recently dead, resighted
    # 5) dead in previous interval
    #--------------------------------------
    #####################################################################################
    
    Psi[1,1] <- phi * F * gammapp
    Psi[1,2] <- phi * F * (1 - gammapp)
    Psi[1,3] <- phi * (1-F)
    Psi[1,4] <- (1 - phi)*r
    Psi[1,5] <- (1 - phi)*Rp*(1-r)
    Psi[1,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[2,1] <- phi * F * gammap
    Psi[2,2] <- phi * F * (1 - gammap)
    Psi[2,3] <- phi * (1-F)
    Psi[2,4] <- (1 - phi)*r
    Psi[2,5] <- (1 - phi)*Rp*(1-r)
    Psi[2,6] <- (1 - phi)*(1 - Rp)*(1-r)

    Psi[3,1] <- 0
    Psi[3,2] <- 0
    Psi[3,3] <- phi
    Psi[3,4] <- (1 - phi)*r
    Psi[3,5] <- (1 - phi)*Rp*(1-r)
    Psi[3,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[4,1] <- 0
    Psi[4,2] <- 0
    Psi[4,3] <- 0
    Psi[4,4] <- 0
    Psi[4,5] <- 0
    Psi[4,6] <- 1
    
    Psi[5,1] <- 0
    Psi[5,2] <- 0
    Psi[5,3] <- 0
    Psi[5,4] <- 0
    Psi[5,5] <- 0
    Psi[5,6] <- 1

    Psi[6,1] <- 0
    Psi[6,2] <- 0
    Psi[6,3] <- 0
    Psi[6,4] <- 0
    Psi[6,5] <- 0
    Psi[6,6] <- 1    
    
    #####################################################################################
    #--------------------------------------
    # Breeding observations (o ~ Omega):
    # 1) seen breeding, seen non-breeding
    # 2) seen breeding, not seen non-breeding
    # 3) not seen breeding, seen non-breeding
    # 4) not seen breeding, not seen non-breeding
    # 5) shot
    #--------------------------------------
    #####################################################################################
    
    Omega[1,1] <- pstar*R                       
    Omega[1,2] <- pstar*(1 - R)  
    Omega[1,3] <- (1 - pstar)*R 
    Omega[1,4] <- (1 - pstar)*(1 - R)
    Omega[1,5] <- 0 
    
    Omega[2,1] <- 0                       
    Omega[2,2] <- 0   
    Omega[2,3] <- R 
    Omega[2,4] <- (1 - R) 
    Omega[2,5] <- 0 

    Omega[3,1] <- 0                       
    Omega[3,2] <- 0   
    Omega[3,3] <- R 
    Omega[3,4] <- (1 - R) 
    Omega[3,5] <- 0 
    
    Omega[4,1] <- 0    
    Omega[4,2] <- 0 
    Omega[4,3] <- 0 
    Omega[4,4] <- 0 
    Omega[4,5] <- 1
    
    Omega[5,1] <- 0    
    Omega[5,2] <- 0 
    Omega[5,3] <- 1 
    Omega[5,4] <- 0 
    Omega[5,5] <- 0 

    Omega[6,1] <- 0    
    Omega[6,2] <- 0 
    Omega[6,3] <- 0 
    Omega[6,4] <- 1 
    Omega[6,5] <- 0 
    
    
    # Likelihood 
    for (i in 1:n.ind){
    
    z[i, first[i]] <- y[i,first[i]]
    
    for (t in (first[i] + 1):n.occasions){
    
    z[i,t] ~ dcat(Psi[z[i,t-1], ])
    
    y[i,t] ~ dcat(Omega[z[i,t], ])
    
    }}}
    
    
    ",fill = TRUE)
sink()







############################################################################################
# Specify Kendall et al. (2013) model in BUGS language
############################################################################################
sink("Kendall2013.jags")
cat("
    model {
    
    ####################################################################################  
    #--------------------------------------
    # Parameters:
    # phi: survival probability
    # F: site fidelity
    # p: capture probability
    # R: non-breeding resight probability given survival from t to t+1
    # Rp (R'): non-breeding resight probability given mortality during interval t to t+1
    # f: band recovery probability
    #--------------------------------------
    ####################################################################################  
    
    phi ~ dbeta(1,1)
    F ~ dbeta(1,1)
    r ~ dbeta(1,1)
    R ~ dbeta(1,1)
    Rp ~ dbeta(1,1)
    gammapp ~ dbeta(1,1)
    gammap <- gammapp

    for (j in 1:n.sec){
      p[j] ~ dbeta(1,1)
    }

    #####################################################################################
    # model secondary detection probabilities
    #####################################################################################
    for (t in 1:n.occasions){
      for (j in 1:n.sec){
        yes[t,j] ~ dbin(p[j], total[t,j])
      }
    }
    pstar <- 1 - prod(1 - p[])
    
    
    #####################################################################################
    #--------------------------------------
    # Latent states (Psi):
    # 1) alive, available
    # 2) permanent emigraht
    # 3) recently dead, shot
    # 4) recently dead, resighted
    # 5) dead in previous interval
    #--------------------------------------
    #####################################################################################
    
    Psi[1,1] <- phi * F * gammapp
    Psi[1,2] <- phi * F * (1 - gammapp)
    Psi[1,3] <- phi * (1-F)
    Psi[1,4] <- (1 - phi)*r
    Psi[1,5] <- (1 - phi)*Rp*(1-r)
    Psi[1,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[2,1] <- phi * F * gammap
    Psi[2,2] <- phi * F * (1 - gammap)
    Psi[2,3] <- phi * (1-F)
    Psi[2,4] <- (1 - phi)*r
    Psi[2,5] <- (1 - phi)*Rp*(1-r)
    Psi[2,6] <- (1 - phi)*(1 - Rp)*(1-r)

    Psi[3,1] <- 0
    Psi[3,2] <- 0
    Psi[3,3] <- phi
    Psi[3,4] <- (1 - phi)*r
    Psi[3,5] <- (1 - phi)*Rp*(1-r)
    Psi[3,6] <- (1 - phi)*(1 - Rp)*(1-r)
    
    Psi[4,1] <- 0
    Psi[4,2] <- 0
    Psi[4,3] <- 0
    Psi[4,4] <- 0
    Psi[4,5] <- 0
    Psi[4,6] <- 1
    
    Psi[5,1] <- 0
    Psi[5,2] <- 0
    Psi[5,3] <- 0
    Psi[5,4] <- 0
    Psi[5,5] <- 0
    Psi[5,6] <- 1

    Psi[6,1] <- 0
    Psi[6,2] <- 0
    Psi[6,3] <- 0
    Psi[6,4] <- 0
    Psi[6,5] <- 0
    Psi[6,6] <- 1    
    
    #####################################################################################
    #--------------------------------------
    # Breeding observations (o ~ Omega):
    # 1) seen breeding, seen non-breeding
    # 2) seen breeding, not seen non-breeding
    # 3) not seen breeding, seen non-breeding
    # 4) not seen breeding, not seen non-breeding
    # 5) shot
    #--------------------------------------
    #####################################################################################
    
    Omega[1,1] <- pstar*R                       
    Omega[1,2] <- pstar*(1 - R)  
    Omega[1,3] <- (1 - pstar)*R 
    Omega[1,4] <- (1 - pstar)*(1 - R)
    Omega[1,5] <- 0 
    
    Omega[2,1] <- 0                       
    Omega[2,2] <- 0   
    Omega[2,3] <- R 
    Omega[2,4] <- (1 - R) 
    Omega[2,5] <- 0 

    Omega[3,1] <- 0                       
    Omega[3,2] <- 0   
    Omega[3,3] <- R 
    Omega[3,4] <- (1 - R) 
    Omega[3,5] <- 0 
    
    Omega[4,1] <- 0    
    Omega[4,2] <- 0 
    Omega[4,3] <- 0 
    Omega[4,4] <- 0 
    Omega[4,5] <- 1
    
    Omega[5,1] <- 0    
    Omega[5,2] <- 0 
    Omega[5,3] <- 1 
    Omega[5,4] <- 0 
    Omega[5,5] <- 0 

    Omega[6,1] <- 0    
    Omega[6,2] <- 0 
    Omega[6,3] <- 0 
    Omega[6,4] <- 1 
    Omega[6,5] <- 0 
    
    
    # Likelihood 
    for (i in 1:n.ind){
    
    z[i, first[i]] <- y[i,first[i]]
    
    for (t in (first[i] + 1):n.occasions){
    
    z[i,t] ~ dcat(Psi[z[i,t-1], ])
    
    y[i,t] ~ dcat(Omega[z[i,t], ])
    
    }}}
    
    
    ",fill = TRUE)
sink()



#############################################################################################
# specify parameters for the simulation run
#############################################################################################
#######################################################
# simulate the data
#######################################################
n.occasions <- 15
n <- c(100, rep(50, n.occasions-1))
n.sec <- 3
n.sims <- 100

######################################################################################
# MCMC settings for models
######################################################################################
n.chains <- 6
n.thin <- 5
n.adapt <- 1000
n.iter <- 25000
n.burnin <- 15000

###
### arrays to store results
###
sim.params <- matrix(NA, 8, n.sims)
results.barker <- array(NA, dim = c(8, n.sims, 3))
results.lindberg <- array(NA, dim = c(8, n.sims, 3))
results.kendall <- array(NA, dim = c(8, n.sims, 3))


#######################################################
# Barker Robust Design Parameters
#######################################################
for (ii in 1:n.sims){
print(ii)
print(Sys.time())

# simulate parameter values
phi <- runif(1, 0.7, 0.9)
F <- runif(1, 0.8, 0.95)
gammapp <- runif(1, 0.6, 0.95)
gammap <- gammapp
p <- runif(n.sec, 0.15, 0.35)
pstar <- 1 - prod(1 - p)
r <- runif(1, 0.15, 0.35)
R <- runif(1, 0.15, 0.35)
Rp <- R/2



# store simulated parameter values
sim.params[1,ii] <- phi
sim.params[2,ii] <- F
sim.params[3,ii] <- gammapp
sim.params[4,ii] <- gammap
sim.params[5,ii] <- pstar
sim.params[6,ii] <- r
sim.params[7,ii] <- R
sim.params[8,ii] <- Rp


# create latent state matrix, temporal variation can be added with a loop, but not dealing with
# that for simple simulations
S <- matrix(c(phi*gammapp*F, phi*(1-gammapp)*F, phi*(1-F), 1-phi, 0,
              phi*gammap*F, phi*(1-gammap)*F, phi*(1-F), 1-phi, 0,
              0, 0, phi, 1-phi, 0,
              0, 0, 0, 0, 1,
              0, 0, 0, 0, 1), 5, 5, byrow = TRUE)

# ensure state-transition probabilities sum to 1
# rowSums(S)

# create a true latent state matrix, fill first available for capture with 1's
true <- matrix(NA, sum(n), n.occasions+1)
first <- NULL
for (j in 1:(n.occasions)){
  true[(1 + sum( n[0:(j-1)] )) : sum(n[1:j]), j] <- 1
  first[(1 + sum( n[0:(j-1)] )) : sum(n[1:j])] <- j
}

# simulate true process
for (i in 1:nrow(true)){
  for (j in (first[i]+1):(n.occasions+1)){
    true[i,j] <- which(rmultinom(1, 1, S[true[i,j-1],]) == 1)
  }}


# secondary occasion array
# this needs to get fixed to account for multiple different numbers of secondary occasions
sec <- array(0, dim = c(nrow(true), n.occasions, n.sec))

# observe secondary occasion data
for (i in 1:nrow(true)){
  for (j in (first[i]):n.occasions){
    for (k in 1:n.sec){
      if (true[i,j] == 1){sec[i,j,k] <- rbinom(1,1,p[k])}
    }}}

# primary occasion, binary
prim <- apply(sec, c(1,2), sum)
prim[prim > 1] <- 1


#######################################################
# get.first, write over first!
#######################################################
true.first <- first
get.first <- function(x){min(which(x != 0))}
first <- apply(prim, 1, get.first)

#######################################################
# clean up...
#######################################################
sec <- sec[which(first != Inf),,]
prim <- subset(prim, first != Inf)
true.first <- subset(true.first, first != Inf)
true <- subset(true, first != Inf)
first <- subset(first, first != Inf)

# non-breeding season array
non <- matrix(0, nrow(true), n.occasions)
# for recently dead, 1 = seen R', 2 = not seen or recovered, 3 = shot...
pr.rec.dead <- c(Rp*(1-r), (1-Rp)*(1-r), r)
for (i in 1:nrow(true)){
  for (j in (first[i]):(n.occasions)){
    if(true[i,j+1] == 1 | true[i,j+1] == 2 | true[i,j+1] == 3){non[i,j] <- rbinom(1, 1, R)}
    if(true[i,j+1] == 4){non[i,j] <- which(rmultinom(1, 1, pr.rec.dead) == 1)}
    if(true[i,j+1] == 5){non[i,j] <- 2}
  }}
non[non == 0] <- 2
table(non)





#######################################################
# robust design data pulled from sec...
#######################################################
test <- apply(sec, c(1,2), sum)
seen <- array(0, dim = c(dim(sec)[1], dim(sec)[2], dim(sec)[3]))
missed <- array(0, dim = c(dim(sec)[1], dim(sec)[2], dim(sec)[3]))

for (i in 1:dim(seen)[1]){
  for (j in 1:dim(seen)[2]){
    for (k in 1:dim(seen)[3]){
      if (test[i,j] > 1 & sec[i,j,k] == 1){seen[i,j,k] <- 1}
      if (test[i,j] >= 1 & sec[i,j,k] == 0){missed[i,j,k] <- 1}
    }
  }
}

yes <- matrix(NA, dim(seen)[2], dim(seen)[3])
no <- matrix(NA, dim(seen)[2], dim(seen)[3])

for (i in 1:nrow(yes)){
  for (j in 1:ncol(yes)){
    yes[i,j] <- sum(seen[,i,j])
    no[i,j] <- sum(missed[,i,j])
  }
}

total <- yes + no


#######################################################
# clean up breeding observations
#######################################################
prim[prim > 1] <- 1
prim[prim == 0] <- 2

sec <- sec[which(first < n.occasions),,]
prim <- subset(prim, first < n.occasions)
true.first <- subset(true.first, first < n.occasions)
non <- subset(non, first < n.occasions)
true <- subset(true, first < n.occasions)
first <- subset(first, first < n.occasions)

summary(as.factor(non[,n.occasions]))

############################################################################################
# combine breeding and non-breeding resights...
# five states:
# 1) seen breeding, seen non-breeding
# 2) seen breeding, not seen non-breeding
# 3) not seen breeding, seen non-breeding
# 4) not seen breeding, not seen non-breeding
# 5) seen, then shot
# 6) not seen, then shot
############################################################################################




############################################################################################
# create barker capture history and initial values
############################################################################################
barker.ch <- matrix(4, nrow(prim), ncol(prim))

for (i in 1:nrow(barker.ch)){
  barker.ch[i,first[i]] <- 1
  for (j in (first[i]+1):ncol(prim)){
    if(prim[i,j] == 1 & non[i,j-1] == 1){barker.ch[i,j] <- 1} # pR
    if(prim[i,j] == 1 & non[i,j-1] == 2){barker.ch[i,j] <- 2} # p(1-R)
    if(prim[i,j] == 2 & non[i,j-1] == 1){barker.ch[i,j] <- 3} # (1-p)R
    if(non[i,j-1] == 3){barker.ch[i,j] <- 5}                  # f
    # default is (1-p)(1-R)
  }
  barker.ch[i,1:first[i]] <- NA
  barker.ch[i,first[i]] <- 1
}


########################################################################
# Barker initial values...
########################################################################
z.init.barker <- matrix(1, nrow(barker.ch), ncol(barker.ch))
for (i in 1:nrow(barker.ch)){
  for (j in first[i]:ncol(barker.ch)){
    if(barker.ch[i,j] == 5){z.init.barker[i,j] <- 4}
    if(barker.ch[i,j] == 5 & j < ncol(barker.ch)){z.init.barker[i,(j+1):(ncol(barker.ch))] <- 6}
  }
  z.init.barker[i,(1:first[i])] <- NA
}









############################################################################################
# create lindberg capture history and initial values
############################################################################################
lindberg.ch <- matrix(4, nrow(prim), ncol(prim))

for (i in 1:nrow(barker.ch)){
  lindberg.ch[i,first[i]] <- 1
  for (j in (first[i]+1):ncol(prim)){
    if(prim[i,j] == 1){lindberg.ch[i,j] <- 2} # pstar
    if(non[i,j-1] == 3){lindberg.ch[i,j] <- 5}                  # f
    # default is 2
  }
  lindberg.ch[i,1:first[i]] <- NA
  lindberg.ch[i,first[i]] <- 1
}


########################################################################
# lindberg initial values...
########################################################################
z.init.lindberg <- matrix(1, nrow(lindberg.ch), ncol(lindberg.ch))
for (i in 1:nrow(lindberg.ch)){
  for (j in first[i]:ncol(lindberg.ch)){
    if(lindberg.ch[i,j] == 5){z.init.lindberg[i,j] <- 4}
    if(lindberg.ch[i,j] == 5 & j < ncol(lindberg.ch)){z.init.lindberg[i,(j+1):(ncol(lindberg.ch))] <- 6}
  }
  z.init.lindberg[i,(1:first[i])] <- NA
}

########################################################################
# lindberg z data
########################################################################
get.last <- function(x){max(which(x == 1))}






############################################################################################
# create kendall capture history and initial values
############################################################################################
kendall.ch <- matrix(4, nrow(prim), ncol(prim))

for (i in 1:nrow(barker.ch)){
  kendall.ch[i,first[i]] <- 1
  for (j in (first[i]+1):ncol(prim)){
    if(prim[i,j] == 1 & non[i,j-1] == 1){kendall.ch[i,j] <- 1} # p*R
    if(prim[i,j] == 1 & non[i,j-1] == 2){kendall.ch[i,j] <- 2} # p*(1-R)
    if(prim[i,j] == 2 & non[i,j-1] == 1){kendall.ch[i,j] <- 3} # (1-p*)R
    if(non[i,j-1] == 3){kendall.ch[i,j] <- 5}                  # f
    # default is (1-p)(1-R)
  }
  kendall.ch[i,1:first[i]] <- NA
  kendall.ch[i,first[i]] <- 1
}


########################################################################
# kendall initial values...
########################################################################
z.init.kendall <- matrix(1, nrow(kendall.ch), ncol(kendall.ch))
for (i in 1:nrow(kendall.ch)){
  for (j in first[i]:ncol(kendall.ch)){
    if(kendall.ch[i,j] == 5){z.init.kendall[i,j] <- 4}
    if(kendall.ch[i,j] == 5 & j < ncol(kendall.ch)){z.init.kendall[i,(j+1):(ncol(kendall.ch))] <- 6}
  }
  z.init.kendall[i,(1:first[i])] <- NA
}



######################################################################################
# Bundle data
######################################################################################
dat.barker <- list(y = barker.ch, n.ind = dim(barker.ch)[1], n.occasions = dim(barker.ch)[2], first = first)
dat.lindberg <- list(y = lindberg.ch, n.ind = dim(lindberg.ch)[1], n.occasions = dim(lindberg.ch)[2], 
                     first = first, yes = yes, n.sec = n.sec, total = total)
dat.kendall <- list(y = kendall.ch, n.ind = dim(kendall.ch)[1], n.occasions = dim(kendall.ch)[2], 
                    first = first, yes = yes, n.sec = n.sec, total = total)

######################################################################################
# Initial values
######################################################################################
init.barker <- function(){list(phi = phi, R = R, Rp = Rp, f = first, p = pstar, z = z.init.barker)}    
init.lindberg <- function(){list(phi = phi, f = first, p = p, z = z.init.lindberg)} 
init.kendall <- function(){list(phi = phi, R = R, Rp = Rp, f = first, p = p, z = z.init.kendall)}


######################################################################################
# Parameters monitored
######################################################################################
pars <- c('phi', 'F', 'R', 'Rp', 'r','gammapp','gammap','pstar')



######################################################################################
# compile model
######################################################################################
require(jagsUI)
# ptm <- proc.time()
Sys.time()
m1 <- jags(dat.barker, init.barker, pars, "barker1997.jags", 
               n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, 
               parallel = T)
# proc.time() - ptm
results.barker[1,ii,1] <- m1$mean$phi; results.barker[1,ii,2] <- m1$q2.5$phi; results.barker[1,ii,3] <- m1$q97.5$phi
results.barker[2,ii,1] <- m1$mean$F; results.barker[2,ii,2] <- m1$q2.5$F; results.barker[2,ii,3] <- m1$q97.5$F
results.barker[3,ii,1] <- m1$mean$gammapp; results.barker[3,ii,2] <- m1$q2.5$gammapp; results.barker[3,ii,3] <- m1$q97.5$gammapp
results.barker[4,ii,1] <- m1$mean$gammap; results.barker[4,ii,2] <- m1$q2.5$gammap; results.barker[4,ii,3] <- m1$q97.5$gammap
results.barker[5,ii,1] <- m1$mean$pstar; results.barker[5,ii,2] <- m1$q2.5$pstar; results.barker[5,ii,3] <- m1$q97.5$pstar
results.barker[6,ii,1] <- m1$mean$r; results.barker[6,ii,2] <- m1$q2.5$r; results.barker[6,ii,3] <- m1$q97.5$r
results.barker[7,ii,1] <- m1$mean$R; results.barker[7,ii,2] <- m1$q2.5$R; results.barker[7,ii,3] <- m1$q97.5$R
results.barker[8,ii,1] <- m1$mean$Rp; results.barker[8,ii,2] <- m1$q2.5$Rp; results.barker[8,ii,3] <- m1$q97.5$Rp


Sys.time()
m2 <- jags(dat.lindberg, init.lindberg, pars, "lindberg2001.jags", 
           n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, 
           parallel = T)

results.lindberg[1,ii,1] <- m2$mean$phi; results.lindberg[1,ii,2] <- m2$q2.5$phi; results.lindberg[1,ii,3] <- m2$q97.5$phi
results.lindberg[2,ii,1] <- m2$mean$F; results.lindberg[2,ii,2] <- m2$q2.5$F; results.lindberg[2,ii,3] <- m2$q97.5$F
results.lindberg[3,ii,1] <- m2$mean$gammapp; results.lindberg[3,ii,2] <- m2$q2.5$gammapp; results.lindberg[3,ii,3] <- m2$q97.5$gammapp
results.lindberg[4,ii,1] <- m2$mean$gammap; results.lindberg[4,ii,2] <- m2$q2.5$gammap; results.lindberg[4,ii,3] <- m2$q97.5$gammap
results.lindberg[5,ii,1] <- m2$mean$pstar; results.lindberg[5,ii,2] <- m2$q2.5$pstar; results.lindberg[5,ii,3] <- m2$q97.5$pstar
results.lindberg[6,ii,1] <- m2$mean$r; results.lindberg[6,ii,2] <- m2$q2.5$r; results.lindberg[6,ii,3] <- m2$q97.5$r
results.lindberg[7,ii,1] <- m2$mean$R; results.lindberg[7,ii,2] <- m2$q2.5$R; results.lindberg[7,ii,3] <- m2$q97.5$R
results.lindberg[8,ii,1] <- m2$mean$Rp; results.lindberg[8,ii,2] <- m2$q2.5$Rp; results.lindberg[8,ii,3] <- m2$q97.5$Rp

Sys.time()
m3 <- jags(dat.kendall, init.kendall, pars, "kendall2013.jags", 
           n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, 
           parallel = T)
results.kendall[1,ii,1] <- m3$mean$phi; results.kendall[1,ii,2] <- m3$q2.5$phi; results.kendall[1,ii,3] <- m3$q97.5$phi
results.kendall[2,ii,1] <- m3$mean$F; results.kendall[2,ii,2] <- m3$q2.5$F; results.kendall[2,ii,3] <- m3$q97.5$F
results.kendall[3,ii,1] <- m3$mean$gammapp; results.kendall[3,ii,2] <- m3$q2.5$gammapp; results.kendall[3,ii,3] <- m3$q97.5$gammapp
results.kendall[4,ii,1] <- m3$mean$gammap; results.kendall[4,ii,2] <- m3$q2.5$gammap; results.kendall[4,ii,3] <- m3$q97.5$gammap
results.kendall[5,ii,1] <- m3$mean$pstar; results.kendall[5,ii,2] <- m3$q2.5$pstar; results.kendall[5,ii,3] <- m3$q97.5$pstar
results.kendall[6,ii,1] <- m3$mean$r; results.kendall[6,ii,2] <- m3$q2.5$r; results.kendall[6,ii,3] <- m3$q97.5$r
results.kendall[7,ii,1] <- m3$mean$R; results.kendall[7,ii,2] <- m3$q2.5$R; results.kendall[7,ii,3] <- m3$q97.5$R
results.kendall[8,ii,1] <- m3$mean$Rp; results.kendall[8,ii,2] <- m3$q2.5$Rp; results.kendall[8,ii,3] <- m3$q97.5$Rp

######################################################################################
# end loop
######################################################################################

}
Sys.time()

###
### save results
###
save.image("O:/Methods/barker_robust_design/latex/simulation_and_results/results_final.RData")


#############################################################################
###
### calculate values for Table 4... MSD and coverage
###
#############################################################################
###
### survival
###
mean(results.barker[1,,1] - sim.params[1,])
length(which(results.barker[1,,2] < sim.params[1,] | sim.params[1,] < results.barker[1,,3]))/n.sims
mean(results.lindberg[1,,1] - sim.params[1,])
length(which(results.lindberg[1,,2] < sim.params[1,] | sim.params[1,] < results.lindberg[1,,3]))/n.sims
mean(results.kendall[1,,1] - sim.params[1,])
length(which(results.kendall[1,,2] < sim.params[1,] & sim.params[1,] < results.kendall[1,,3]))/n.sims

###
### fidelity
###
mean(results.barker[2,,1] - sim.params[2,])
length(which(results.barker[2,,2] < sim.params[2,] & sim.params[2,] < results.barker[2,,3]))/n.sims
mean(results.lindberg[2,,1] - sim.params[2,])
length(which(results.lindberg[2,,2] < sim.params[2,] & sim.params[2,] < results.lindberg[2,,3]))/n.sims
mean(results.kendall[2,,1] - sim.params[2,])
length(which(results.kendall[2,,2] < sim.params[2,] & sim.params[2,] < results.kendall[2,,3]))/n.sims

###
### breeding probability
###
mean(results.lindberg[3,,1] - sim.params[3,])
length(which(results.lindberg[3,,2] < sim.params[3,] & sim.params[3,] < results.lindberg[3,,3]))/n.sims
mean(results.kendall[3,,1] - sim.params[3,])
length(which(results.kendall[3,,2] < sim.params[3,] & sim.params[3,] < results.kendall[3,,3]))/n.sims

###
### primary occasion detection probability
###
mean(results.barker[5,,1] - sim.params[5,])
length(which(results.barker[5,,2] < sim.params[5,] & sim.params[5,] < results.barker[5,,3]))/n.sims
mean(results.lindberg[5,,1] - sim.params[5,])
length(which(results.lindberg[5,,2] < sim.params[5,] & sim.params[5,] < results.lindberg[5,,3]))/n.sims
mean(results.kendall[5,,1] - sim.params[5,])
length(which(results.kendall[5,,2] < sim.params[5,] & sim.params[5,] < results.kendall[5,,3]))/n.sims


###
### recovery probability
###
mean(results.barker[6,,1] - sim.params[6,])
length(which(results.barker[6,,2] < sim.params[6,] & sim.params[6,] < results.barker[6,,3]))/n.sims
mean(results.lindberg[6,,1] - sim.params[6,])
length(which(results.lindberg[6,,2] < sim.params[6,] & sim.params[6,] < results.lindberg[6,,3]))/n.sims
mean(results.kendall[6,,1] - sim.params[6,])
length(which(results.kendall[6,,2] < sim.params[6,] & sim.params[6,] < results.kendall[6,,3]))/n.sims


###
### resighting probability (R)
###
mean(results.barker[7,,1] - sim.params[7,])
length(which(results.barker[7,,2] < sim.params[7,] & sim.params[7,] < results.barker[7,,3]))/n.sims
mean(results.kendall[7,,1] - sim.params[7,])
length(which(results.kendall[7,,2] < sim.params[7,] & sim.params[7,] < results.kendall[7,,3]))/n.sims

###
### resighting probability (R')
###
mean(results.barker[8,,1] - sim.params[8,])
length(which(results.barker[8,,2] < sim.params[8,] & sim.params[8,] < results.barker[8,,3]))/n.sims
mean(results.kendall[8,,1] - sim.params[8,])
length(which(results.kendall[8,,2] < sim.params[8,] & sim.params[8,] < results.kendall[8,,3]))/n.sims

#############################################################################
# Figure 1
#############################################################################
pdf("O:/Methods/barker_robust_design/latex/Ecosphere_final/Figure1.pdf", width=8, height=14)
par(mfrow = c(5,3), oma = c(2.1,4.1,3.1,2.1), mar = c(5.1,2.1,2.1,2.1), family = 'serif')

#survival, top row
plot(results.barker[1,,1] ~ sim.params[1,], ylab = '', xlab = expression(phi['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[1,], results.barker[1,,2], sim.params[1,], results.barker[1,,3], length = 0, lty = 2)
mtext(expression(phi['estimate']), side = 2, line = 3, cex = 1.5)
mtext(expression(Barker~1997), side = 3, line = 2, cex = 1.25)

plot(results.lindberg[1,,1] ~ sim.params[1,], ylab = '', xlab = expression(phi['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[1,], results.lindberg[1,,2], sim.params[1,], results.lindberg[1,,3], length = 0, lty = 2)
mtext(expression(Lindberg~italic(et~al.)~2001), side = 3, line = 1.75, cex = 1.25)

plot(results.kendall[1,,1] ~ sim.params[1,], ylab = '', xlab = expression(phi['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[1,], results.kendall[1,,2], sim.params[1,], results.kendall[1,,3], length = 0, lty = 2)
mtext(expression(Kendall~italic(et~al.)~2013), side = 3, line = 2, cex = 1.25)


#Fidelity, 2nd row
plot(results.barker[2,,1] ~ sim.params[2,], ylab = '', xlab = expression(F['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[2,], results.barker[2,,2], sim.params[2,], results.barker[2,,3], length = 0, lty = 2)
mtext(expression(F['estimate']), side = 2, line = 3, cex = 1.5)

plot(results.lindberg[2,,1] ~ sim.params[2,], ylab = '', xlab = expression(F['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[2,], results.lindberg[2,,2], sim.params[2,], results.lindberg[2,,3], length = 0, lty = 2)

plot(results.kendall[2,,1] ~ sim.params[2,], ylab = '', xlab = expression(F['truth']), 
     pch = 19, xlim = c(0.6,1), ylim = c(0.6,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[2,], results.kendall[2,,2], sim.params[2,], results.kendall[2,,3], length = 0, lty = 2)

#gamma, 3rd row
plot(rep(0,n.sims) ~ sim.params[3,], ylab = '', xlab = expression(gamma['truth']), 
     pch = 19, xlim = c(0.4,1), ylim = c(0.4,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
mtext(expression(gamma['estimate']), side = 2, line = 3, cex = 1.5)

plot(results.lindberg[3,,1] ~ sim.params[3,], ylab = '', xlab = expression(gamma['truth']), 
     pch = 19, xlim = c(0.4,1), ylim = c(0.4,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[3,], results.lindberg[3,,2], sim.params[3,], results.lindberg[3,,3], length = 0, lty = 2)

plot(results.kendall[3,,1] ~ sim.params[3,], ylab = '', xlab = expression(gamma['truth']), 
     pch = 19, xlim = c(0.4,1), ylim = c(0.4,1), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[3,], results.kendall[3,,2], sim.params[3,], results.kendall[3,,3], length = 0, lty = 2)


#detection, 4th row
plot(results.barker[5,,1] ~ sim.params[5,], ylab = '', xlab = expression('p*'['truth']), 
     pch = 19, xlim = c(0.2,0.8), ylim = c(0.2,0.8), cex.lab = 1.5, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[5,], results.barker[5,,2], sim.params[5,], results.barker[5,,3], length = 0, lty = 2)
mtext(expression('p*'['estimate']), side = 2, line = 3, cex = 1.5)

plot(results.lindberg[5,,1] ~ sim.params[5,], ylab = '', xlab = expression('p*'['truth']), 
     pch = 19, xlim = c(0.2,0.8), ylim = c(0.2,0.8), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[5,], results.lindberg[5,,2], sim.params[5,], results.lindberg[5,,3], length = 0, lty = 2)

plot(results.kendall[5,,1] ~ sim.params[5,], ylab = '', xlab = expression('p*'['truth']), 
     pch = 19, xlim = c(0.2,0.8), ylim = c(0.2,0.8), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[5,], results.kendall[5,,2], sim.params[5,], results.kendall[5,,3], length = 0, lty = 2)

#band recovery, 5th row
plot(results.barker[6,,1] ~ sim.params[6,], ylab = '', xlab = expression(r['truth']), 
     pch = 19, xlim = c(0,0.6), ylim = c(0,0.6), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[6,], results.barker[6,,2], sim.params[6,], results.barker[6,,3], length = 0, lty = 2)
mtext(expression(r['estimate']), side = 2, line = 3, cex = 1.5)

plot(results.lindberg[6,,1] ~ sim.params[6,], ylab = '', xlab = expression(r['truth']), 
     pch = 19, xlim = c(0,0.6), ylim = c(0,0.6), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[6,], results.lindberg[6,,2], sim.params[6,], results.lindberg[6,,3], length = 0, lty = 2)

plot(results.kendall[6,,1] ~ sim.params[6,], ylab = '', xlab = expression(r['truth']), 
     pch = 19, xlim = c(0,0.6), ylim = c(0,0.6), cex.lab = 2, cex.axis = 1.5, las = 1); abline(0,1, lty = 2)
# arrows(sim.params[6,], results.kendall[6,,2], sim.params[6,], results.kendall[6,,3], length = 0, lty = 2)

dev.off()





#############################################################################
# Figure 2
#############################################################################
pdf("O:/Methods/barker_robust_design/latex/Ecosphere_final/Figure2.pdf", width=7, height=8)

par(mfrow = c(2,2), oma = c(2.1,5.1,5.1,2.1), mar = c(5.1,2.1,2.1,2.1), family = 'serif')

plot(results.barker[7,,1] ~ sim.params[7,], ylab = '', xlab = expression(R['truth']), 
     pch = 19, xlim = c(0.1,0.4), ylim = c(0.1,0.4), cex.lab = 1.5, las = 1); abline(0,1, lty = 2)
#arrows(sim.params[7,], results.barker[7,,2], sim.params[7,], results.barker[7,,3], length = 0, lty = 2)
mtext(expression(R['estimate']), side = 2, line = 3, cex = 1.5)
mtext(expression(Barker~1997), side = 3, line = 2, cex = 1.25)

plot(results.kendall[7,,1] ~ sim.params[7,], ylab = '', xlab = expression(R['truth']), 
     pch = 19, xlim = c(0.1,0.4), ylim = c(0.1,0.4), cex.lab = 1.5, las = 1); abline(0,1, lty = 2)
#arrows(sim.params[7,], results.kendall[7,,2], sim.params[7,], results.kendall[7,,3], length = 0, lty = 2)
mtext(expression(Kendall~italic(et~al.)~2013), side = 3, line = 2, cex = 1.25)

plot(results.barker[8,,1] ~ sim.params[8,], ylab = '', xlab = expression(R['truth']^"'"), 
     pch = 19, xlim = c(0,0.4), ylim = c(0,0.4), cex.lab = 1.5, las = 1); abline(0,1, lty = 2)
#arrows(sim.params[8,], results.barker[8,,2], sim.params[8,], results.barker[7,,3], length = 0, lty = 2)
mtext(expression(R['estimate']^"'"), side = 2, line = 3, cex = 1.5)

plot(results.kendall[8,,1] ~ sim.params[8,], ylab = '', xlab = expression(R['truth']^"'"), 
     pch = 19, xlim = c(0,0.4), ylim = c(0,0.4), cex.lab = 1.5, las = 1); abline(0,1, lty = 2)
#arrows(sim.params[8,], results.kendall[8,,2], sim.params[8,], results.kendall[7,,3], length = 0, lty = 2)

dev.off()


