library(rjags);library(MCMCpack)
hist(rnorm(1000,0,6), col="grey", main="~ N(0,6)")
### check week vs. tank covariation

# minimum effective sample size for these models
library(mcmcse)
minESS(8)

##### Model Code #####
model_string<- "model {
  # likelihood
  for (i in 1:length(res)) { # loop over mesocosms / observations
    # varying intercept/slope for every mesocosm
    mu[i] <- beta[treat[i], week[i]] + alpha[tank[i]]
    res[i] ~ dnorm(mu[i], tau) #try a gamma dist
  }
  # priors
  tau <- pow(sigma,-2) #precision on res
  sigma ~ dgamma(0.001,.001) #variance on res

  for (u in 1:U){
      alpha[u] ~ dnorm(mu.a[u], tau.a[u]) # potential slope, excluding for now
      mu.a[u] ~ dnorm(0, 6)
      tau.a[u] <- pow(sigma.a[u], -2)
      sigma.a[u] ~ dgamma(1,1)
  } 

  # Control mean response variables
  B0w1 <- mean(beta[1,1]) #control week 1
  B0w2 <- mean(beta[1,2]) #control week 2
  B0w3 <- mean(beta[1,3]) #control week 3
  B0w4 <- mean(beta[1,4]) #control week 4
  B0w5 <- mean(beta[1,5]) #control week 5
  B0w6 <- mean(beta[1,6]) #control week 6
  B0w7 <- mean(beta[1,7]) #control week 7
  # Live mean response variables
  BLw1 <- mean(beta[2,1]) #live week 1
  BLw2 <- mean(beta[2,2]) #live week 2
  BLw3 <- mean(beta[2,3]) #live week 3
  BLw4 <- mean(beta[2,4]) #live week 4
  BLw5 <- mean(beta[2,5]) #live week 5
  BLw6 <- mean(beta[2,6]) #live week 6
  BLw7 <- mean(beta[2,7]) #live week 7

  for (j in 1:J) { # loop over treatments
    for (s in 1:S) { # loop over weeks
      # so beta[j,s] will be the slope for each treat at each week
      beta[j, s] ~ dnorm(mu_b[j], tau_b[j])
    }
    # response ratio for control treatments
    rCw1[j] <- (mean(beta[j,1]))/ifelse(B0w1==0, B0w1+.01,B0w1) 
    rCw2[j] <- (mean(beta[j,2]))/ifelse(B0w2==0, B0w2+.01,B0w2) 
    rCw3[j] <- (mean(beta[j,3]))/ifelse(B0w3==0, B0w3+.01,B0w3) 
    rCw5[j] <- (mean(beta[j,5]))/ifelse(B0w5==0, B0w5+.01,B0w5) 
    rCw6[j] <- (mean(beta[j,6]))/ifelse(B0w6==0, B0w6+.01,B0w6) 
    rCw7[j] <- (mean(beta[j,7]))/ifelse(B0w7==0, B0w7+.01,B0w7) 
   
    # response ratio for live treatments
    rLw1[j] <- (mean(beta[j,1]))/ifelse(BLw1==0, BLw1+.01,BLw1) 
    rLw2[j] <- (mean(beta[j,2]))/ifelse(BLw2==0, BLw2+.01,BLw2) 
    rLw3[j] <- (mean(beta[j,3]))/ifelse(BLw3==0, BLw3+.01,BLw3) 
    rLw5[j] <- (mean(beta[j,5]))/ifelse(BLw5==0, BLw5+.01,BLw5) 
    rLw6[j] <- (mean(beta[j,6]))/ifelse(BLw6==0, BLw6+.01,BLw6) 
    rLw7[j] <- (mean(beta[j,7]))/ifelse(BLw7==0, BLw7+.01,BLw7) 
    
    #hyper parameter priors
    mu_b[j] ~ dnorm (0, 6)
    tau_b[j] <- pow(sigma_b[j], -2)
    sigma_b[j] ~ dgamma(1,1)
  }
    rCw4[2] <- (mean(beta[2,4]))/ifelse(B0w4==0, B0w4+.01,B0w4)
    rCw4[3] <- (mean(beta[3,4]))/ifelse(B0w4==0, B0w4+.01,B0w4)   
    rLw4[1] <- (mean(beta[1,4]))/ifelse(BLw4==0, BLw4+.01,BLw4) 
    rLw4[3] <- (mean(beta[3,4]))/ifelse(BLw4==0, BLw4+.01,BLw4) 

    # BACI response ratio for the control treatment
    meanBACIaC[2] <- (rCw4[2]+rCw5[2]+rCw6[2]+rCw7[2])/4
    meanBACIbC[2] <- (rCw1[2]+rCw2[2]+rCw3[2])/3
    BACIc[2] <- meanBACIaC[2]/ifelse(meanBACIbC[2]==0,meanBACIbC[2]+.001,meanBACIbC[2])
    meanBACIaC[3]<- (rCw4[3]+rCw5[3]+rCw6[3]+rCw7[3])/4
    meanBACIbC[3]<- (rCw1[3]+rCw2[3]+rCw3[3])/3
    BACIc[3]<-meanBACIaC[3]/ifelse(meanBACIbC[3]==0,meanBACIbC[3]+.001,meanBACIbC[3])

  # BACI response ratio for the live treatment
    meanBACIaL[1] <- (rLw4[1]+rLw5[1]+rLw6[1]+rLw7[1])/4
    meanBACIbL[1] <- (rLw1[1]+rLw2[1]+rLw3[1])/3
    BACIl[1] <- meanBACIaL[1]/ifelse(meanBACIbL[1]==0,meanBACIbL[1]+.001,meanBACIbL[1])
    meanBACIaL[3]<- (rLw4[3]+rLw5[3]+rLw6[3]+rLw7[3])/4
    meanBACIbL[3]<- (rLw1[3]+rLw2[3]+rLw3[3])/3
    BACIl[3]<-meanBACIaL[3]/ifelse(meanBACIbL[3]==0,meanBACIbL[3]+.001,meanBACIbL[3])
}"

##### NH3 Bayesian Models #####
source('waterchem.R')
head(WaterNutrients);nrow(WaterNutrients)
WNmodel<-WaterNutrients %>% dplyr::select(FilterdNH3ugL,FilterdSRPugL,
                                          Day, NewTreat,Tank) %>%
  mutate(logFNH4=log10(FilterdNH3ugL),
         logFSRP=log10(FilterdSRPugL),
         TreatF=factor(NewTreat, levels=c("Control","Live","Dead")),
         DayF=factor(Day, levels=c(-20,-15,-3,4,18,25,39)),
         TankF=factor(Tank))
res<- WNmodel %>% pull(logFNH4)
week<- WNmodel %>% pull(DayF)
treat<-WNmodel %>% pull(TreatF)
tank<-WNmodel %>% pull(TankF)
mean(res) # should be beta + alpha
inits<-list(alpha=rep(0.5,18),
            beta=matrix(rep(1,21),nrow=3))
#run the model
NHmodel<-jags.model(textConnection(model_string),
                    data=list(res=res, treat=treat, week=week, tank=tank,
                              J=3, S=7, U=18), 
                    inits=inits,
                    n.chains=3, n.adapt=30000)
update(NHmodel, 10000) # burn in for 2000 samples
NHmcmc<-coda.samples(NHmodel,
                     variable.names=c("BACIc","BACIl", "rCw4","rLw4"), 
                     n.iter=500000, thin=40)
pdf('mcmcDiagnostic/nh4mcmcdiag.pdf')
plot(NHmcmc)
gelman.plot(NHmcmc)
dev.off()

summary(NHmcmc)
NH4mcmc.data<-as.matrix(NHmcmc); colnames(NH4mcmc.data)
write.csv(NH4mcmc.data, "mcmcDiagnostic/Nh4mcmc.csv")

# did we get enough samples to properly model the posterior distribution?
# min ESS gave us 8804
mcse.multi(NH4mcmc.data) #covariance matrix
multiESS(NH4mcmc.data) #multivariate effective sample size
ess(NH4mcmc.data) #univariate effective sample size

### NEED TO CHECK IF I CONVERT BEFORE OR AFTER TO GET PROB CURVE
#probability ammonium increased compared to control treatments
1-ecdf(NH4mcmc.data[,2])( 1 ) #96.7%
1-ecdf(NH4mcmc.data[,2])( 1.3 ) # 51.9% probability of 30% increase
qplot(NH4mcmc.data[,2], geom="histogram", xlim=c(0,5), bins=60,
      ylab="frequency",xlab="BACI d/c", fill= NH4mcmc.data[,2] >= 1)+
  scale_fill_manual(values=c("grey","black"), guide=F)

#### SRP Bayesan Models #####
ressrp<- WNmodel[complete.cases(WNmodel),] %>% pull(logFSRP)
week<- WNmodel[complete.cases(WNmodel),] %>% pull(DayF)
treat<-WNmodel[complete.cases(WNmodel),] %>% pull(TreatF)
tank<-WNmodel[complete.cases(WNmodel),] %>% pull(TankF)
mean(ressrp) # should be beta + alpha
inits<-list(alpha=rep(0.8,18),
            beta=matrix(rep(1.5,21),nrow=3))

#run the model
SRPmodel<-jags.model(textConnection(model_string),
                    data=list(res=ressrp, treat=treat, week=week, tank=tank,
                              J=3, S=7, U=18), 
                    inits=inits,
                    n.chains=3, n.adapt=30000)
update(SRPmodel, 10000) # burn in for 2000 samples
SRPmcmc<-coda.samples(SRPmodel,
                     variable.names=c("BACIc","BACIl", "rCw4","rLw4"), 
                     n.iter=500000, thin=50)
pdf('mcmcDiagnostic/srpmcmcdiag.pdf')
plot(SRPmcmc)
gelman.plot(SRPmcmc)
dev.off()

summary(SRPmcmc)
SRPmcmc.data<-as.matrix(SRPmcmc); colnames(SRPmcmc.data)

# did we get enough samples to properly model the posterior distribution?
# min ESS gave us 8804
mcse.multi(SRPmcmc.data) #covariance matrix
multiESS(SRPmcmc.data) #multivariate effective sample size
ess(SRPmcmc.data) #univariate effective sample size
# higher cov -> lower effective sample size for the SRP data

#probability srp increased compared to control treatments
1-ecdf(SRPmcmc.data[,2])( 1 ) #98.0%
1-ecdf(SRPmcmc.data[,2])( 1.3 ) #40%
qplot(SRPmcmc.data[,2], geom="histogram", xlim=c(0,5), bins=60,
      ylab="frequency",xlab="SRP BACI d/c", fill= SRPmcmc.data[,2] >= 1)+
  scale_fill_manual(values=c("grey","black"), guide=F)
write.csv(SRPmcmc.data, "mcmcDiagnostic/srpmcmc.csv")

##### Benthic Chlorophyll Bayesian Models #####
source("Producers.R")
head(ChlSummary)
# get the data in the appropriate format
# don't want to use log(BenChl) as it gave me a -log distribution 
# (lots near 0, tapers after 1)
Chlmodel<-ChlSummary %>% dplyr::select(BenthicChlA.ug.cm,WaterColChlA.ug.L,
                                          Day, NewTreat,Tank) %>%
  mutate(logBenChl=log10(BenthicChlA.ug.cm),
         logWCChl=log10(WaterColChlA.ug.L),
         TreatF=factor(NewTreat, levels=c("Control","Live","Dead")),
         DayF=factor(Day, levels=c(-20,-15,-3,4,11,25,39)),
         TankF=factor(Tank))
res<- Chlmodel[complete.cases(Chlmodel),] %>% pull(BenthicChlA.ug.cm)
week<- Chlmodel[complete.cases(Chlmodel),] %>% pull(DayF)
treat<-Chlmodel[complete.cases(Chlmodel),] %>% pull(TreatF)
tank<-Chlmodel[complete.cases(Chlmodel),] %>% pull(TankF)
mean(res) # should be beta + alpha
inits<-list(alpha=rep(1,18),
            beta=matrix(rep(2,21),nrow=3))
#run the model
BChlmodel<-jags.model(textConnection(model_string),
                    data=list(res=res, treat=treat, week=week, tank=tank,
                              J=3, S=7, U=18), 
                    #inits=inits,
                    n.chains=3, n.adapt=30000)
update(NHmodel, 10000) # burn in for 2000 samples
BChlmcmc<-coda.samples(BChlmodel,
                     variable.names=c("BACIc","BACIl", "rCw4","rLw4"), 
                     n.iter=500000, thin=40)
pdf('mcmcDiagnostic/bchlmcmcdiag.pdf')
plot(BChlmcmc)
gelman.plot(BChlmcmc)
dev.off()

summary(BChlmcmc)
BChlmcmc.data<-as.matrix(BChlmcmc); colnames(BChlmcmc.data)
write.csv(BChlmcmc.data, "mcmcDiagnostic/BChlmcmc.csv")

# did we get enough samples to properly model the posterior distribution?
# min ESS gave us 8804
mcse.multi(BChlmcmc.data) #covariance matrix
multiESS(BChlmcmc.data) #multivariate effective sample size
ess(BChlmcmc.data) #univariate effective sample size

### NEED TO CHECK IF I CONVERT BEFORE OR AFTER TO GET PROB CURVE
#probability ammonium increased compared to control treatments
1-ecdf(BChlmcmc.data[,2])( 1 ) #66.7%
1-ecdf(BChlmcmc.data[,2])( 1.3 ) #59.6% probability of 30% increase
qplot(BChlmcmc.data[,2], geom="histogram", xlim=c(0,5), bins=60,
      ylab="frequency",xlab="BACI d/c", fill= BChlmcmc.data[,2] >= 1)+
  scale_fill_manual(values=c("grey","black"), guide=F)
#very wide prediction, high cov, likely should add more mcmc iterations


##### Water Column Chlorophyll Bayesian Models #####
##### Model String for three time point samples #####
##### Metabolism Bayesian Models #####
##### Decomposition Bayesian Models #####