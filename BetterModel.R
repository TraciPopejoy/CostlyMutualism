source('waterchem.R')
library(rjags);library(MCMCpack)
head(WaterNutrients)
# TO DO: 
# check likelihood with Dr. Patten (do I need intercept by tank?)
# read more about priors and make sure they are appropriate
# make sure I'm doing the mcmc sampling correctly
# 
yy<-rgamma(100,1,1)
hist(yy)

xx<-runif(1000,.1, 2);hist(res);hist(xx);hist(10^xx)
mean(res);sd(res)
model_string<- "model {
  # likelihood
  for (i in 1:length(res)) { # loop over mesocosms / observations
    # varying intercept/slope for every mesocosm
    mu[i] <- beta[treat[i], week[i]] + alpha[tank[i]]
    res[i] ~ dnorm(mu[i], tau)
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
      # so beta[s,j] will be the slope for each treat at each week
      beta[j, s] ~ dnorm(mu_b[j], tau_b[j])
    }
    # response ratio for control treatments at week 1
    rCw1[j] <- (mean(beta[j,1]))/ifelse(B0w1==0, B0w1+.01,B0w1) 
    rCw2[j] <- (mean(beta[j,2]))/ifelse(B0w2==0, B0w2+.01,B0w2) 
    rCw3[j] <- (mean(beta[j,3]))/ifelse(B0w3==0, B0w3+.01,B0w3) 
    rCw4[j] <- (mean(beta[j,4]))/ifelse(B0w4==0, B0w4+.01,B0w4)
    rCw5[j] <- (mean(beta[j,5]))/ifelse(B0w5==0, B0w5+.01,B0w5) 
    rCw6[j] <- (mean(beta[j,6]))/ifelse(B0w6==0, B0w6+.01,B0w6) 
    rCw7[j] <- (mean(beta[j,7]))/ifelse(B0w7==0, B0w7+.01,B0w7) 
   
    #hyper parameter priors
    mu_b[j] ~ dnorm (0, 6)
    tau_b[j] <- pow(sigma_b[j], -2)
    sigma_b[j] ~ dgamma(1,1)
  }
    # BACI response ratio
    meanBACIa[2] <- (rCw4[2]+rCw5[2]+rCw6[2]+rCw7[2])/4
    meanBACIb[2] <- (rCw1[2]+rCw2[2]+rCw3[2])/3
    BACI[2] <- meanBACIa[2]/ifelse(meanBACIb[2]==0,meanBACIb[2]+.001,meanBACIb[2])
    meanBACIa[3]<- (rCw4[3]+rCw5[3]+rCw6[3]+rCw7[3])/4
    meanBACIb[3]<- (rCw1[3]+rCw2[3]+rCw3[3])/3
    BACI[3]<-meanBACIa[3]/ifelse(meanBACIb[3]==0,meanBACIb[3]+.001,meanBACIb[3])
}"
NHmodel<-WaterNutrients %>% dplyr::select(FilterdNH3ugL,Day, NewTreat,Tank) %>%
  mutate(logFNH4=log10(FilterdNH3ugL),
         TreatF=factor(NewTreat, levels=c("Control","Live","Dead")),
         DayF=factor(Day, levels=c(-20,-15,-3,4,18,25,39)),
         TankF=factor(Tank))
res<- NHmodel %>% pull(logFNH4)
week<- NHmodel %>% pull(DayF)
treat<-NHmodel %>% pull(TreatF)
tank<-NHmodel %>% pull(TankF)
mean(res) # should be beta - alpha
inits<-list(alpha=rep(0.5,18),
  beta=matrix(rep(1,21),nrow=3))

#run the model
model<-jags.model(textConnection(model_string),
                  data=list(res=res, treat=treat, week=week, tank=tank,
                            J=3, S=7, U=18), 
                  inits=inits,
                  n.chains=3, n.adapt=50000)
update(model, 10000) # burn in for 2000 samples
mcmc_samples<-coda.samples(model,
                           variable.names=c("beta", "BACI"), 
                           n.iter=500000, thin=100)
acfplot(mcmc_samples[,1:2])
pdf("trying.pdf")
plot(mcmc_samples[,1:2])
gelman.plot(mcmc_samples[,1:2])
dev.off()
gelman.diag(mcmc_samples[,1:2]) #it fucking worked!
pnorm(abs(geweke.diag(mcmc_samples[3][,1:2])[[1]]$z), lower.tail=F)
#third converged but other two chains did not
heidel.diag(mcmc_samples)
#run was long enough for BACI ratios; failed for some of the betas
summary(mcmc_samples)
test<-as.matrix(mcmc_samples)

#ggplot(NHmodel[NHmodel$Day==-20,], aes(x=TreatF, y=FilterdNH3ugL))+geom_boxplot()
#ggplot(NHmodel[NHmodel$Day==-15,], aes(x=TreatF, y=FilterdNH3ugL))+geom_boxplot()
#ggplot(NHmodel[NHmodel$Day==4,], aes(x=TreatF, y=FilterdNH3ugL))+geom_boxplot()
str(mcmc_samples)
DdC<-cbind(as.matrix(mcmc_samples[1][,2]),
           as.matrix(mcmc_samples[2][,2]),
           as.matrix(mcmc_samples[3][,2]))

1-ecdf(DdC)( 1 ) #probability ammonium increased within the tanks
qplot(DdC, geom="histogram", xlim=c(0,5),ylab="frequency",xlab="BACI d/c")

LdC<-cbind(as.matrix(mcmc_samples[1][,1]),
           as.matrix(mcmc_samples[2][,1]),
           as.matrix(mcmc_samples[3][,1]))

1-ecdf(LdC)( 1 ) #probability ammonium increased within the tanks
qplot(LdC, geom="histogram", xlim=c(0,5),ylab="frequency",xlab="BACI l/c")
