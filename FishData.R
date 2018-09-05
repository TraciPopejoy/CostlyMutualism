library(readxl)
fish<-read_excel("./data/CostMutData.xlsx", sheet="FishData")

treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
library(tidyverse)
library(lubridate)
FishData<-left_join(fish,treat, by="Tank") %>% 
  mutate(LengthChange=DeadLength.mm - StandLength.mm,
         WeightChange=DeadWeight.g - Weight.g,
         InfectionDensALT=InfectionDensity.gloch,
         DaysSurvived=interval(ymd("2018-06-18"),ymd(FishData$DiedALT)) %/% days())
FishData$InfectionRound<-as.factor(FishData$InfectionRound)
#FishData is individual fish in rows with subsequent variables

ggplot(FishData[is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()
ggplot(FishData[!is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()
ggplot(FishData, aes(x=DaysSurvived, y=InfectionDensALT, color=Infected))+
  geom_point()+geom_smooth(method="lm")
ggplot(FishData, aes(x=Died))+
  geom_histogram()+facet_grid(~Treatment+Infected)+
  fungraph+theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1))
ggplot(FishData[FishData$Infected=="Y",], 
       aes(x=InfectionRound, y=InfectionDensity.gloch))+geom_boxplot()+
  ylab("Glochidia per fish")+xlab("Infection round")+fungraph

#Infectiondata provides summary statistics of infection density
Infectiondata<-FishData %>% filter(Infected=="Y") %>% group_by(InfectionRound) %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T),
            meanSurvival=mean(DaysSurvived, na.rm = T),
            sdSurvival=sd(DaysSurvived, na.rm=T))

#Infectiondata provides summary statistics of infection density
FishData %>% filter(Infected=="Y") %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T))

#this adds the mean density to fish that weren't dissected for their gills
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==1,22]<-10
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==2,22]<-3
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==3,22]<-8
#plots survival against infection density
ggplot(FishData, aes(x=DaysSurvived, y=InfectionDensALT, color=Infected))+
  geom_point()+geom_smooth(method="lm")+ylab("Infection Density (gloch per fish)")+
  xlab("Days Survived")+scale_color_manual(values=c("Y"="red","N"="black"))+fungraph
#plots histogram (?) that shows when fish died and their infection density
ggplot(FishData, aes(x=Died,y=Infected, fill=InfectionDensALT))+
  geom_bar(stat="identity")+facet_grid(~Treatment+Infected)+
  scale_fill_continuous("Infection Density",high="red")+fungraph+
  theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
        title=element_text(size=11))
#plots infection density with post mortem weight gain
ggplot(FishData[FishData$Infected!="N",], aes(x=InfectionDensALT, y=WeightChange))+
  geom_point()+geom_smooth(method="lm")+facet_grid(~Treatment+Infected, space="free")+fungraph+
  ylab("Weight Change (g)")+xlab("Infection Density (gloch per fish)")

#DeadPerDayTreat shows the number of fish dead in each treatment/tank per day
DeadPerDayTreat<-FishData %>% group_by(Treatment, Infected, Died) %>% 
  summarize(n=n()) %>% filter(Died >"2018-06-18")

#point graph with infection density and treatment type with weight change
ggplot(FishData, aes(x=Infected, y=WeightChange, color=InfectionDensALT))+
  geom_jitter(size=3, width=.2)+ stat_summary(color="black")+
  theme_bw()+ylab("Post-mortem weight change (grams)")+facet_wrap(~Treatment)
#graph of how many of each type of fish died per day (histogram)
ggplot(DeadPerDayTreat, aes(x=Died, y=n, fill=Infected))+
  geom_bar(stat="identity", position="dodge")+
  ylab("Number of Dead Fish per day")+
  xlab("Day Fish Died, placed Jun 18")+
  facet_wrap(~Treatment)+theme_bw()

### CumDeath is the tank, type of fish, and treatment type and 
### the number of cumualtive dead fish per day
CumDeath<-FishData  %>% 
  filter(is.na(alive)) %>%
  group_by(Treatment,Tank, Infected, DiedALT) %>% 
  filter(DiedALT > "2018-06-19") %>% summarise(nDied=n()) %>% 
  mutate(cumDeath=cumsum(nDied),
         Alive=5-cumDeath,
         TreatType=paste(Treatment,Infected))
#line graph of number of dead fish in a treatment and infection type
ggplot(CumDeathTEST, aes(x=DiedALT, y=cumDeath))+
  geom_smooth(aes(color=Infected), method="lm")+geom_point(position="jitter")+
  ylab("Number of Dead Fish")+
  xlab("Day Fish Died, placed Jun 18")+theme_bw()

#when the fish died based on BOTH treatment type
ggplot(CumDeathTEST, aes(x=DiedALT, y=Alive, color=Treatment))+
  #geom_jitter(alpha=.4, size=3, width=.2)+
  geom_smooth(method="lm", level=.5)+
  ylab("Number of Alive Fish")+fungraph

#when the fish died based on TANK treatment type
ggplot(CumDeathTEST, aes(x=DiedALT, y=Alive, color=Treatment))+
  #geom_jitter(alpha=.4, size=3, width=.2)+
  geom_smooth(method="lm", level=.5)+
  ylab("Number of Alive Fish")+fungraph


FishData$TreatmentType<-paste(FishData$Treatment, FishData$Infected)
ggplot(FishData, aes(x=Died, y=Weight.g, color=TreatmentType))+
  geom_point()+geom_smooth(method="lm", level=.2)
ggplot(FishData, aes(x=Died, y=StandLength.mm, color=TreatmentType))+
  geom_point()+
  geom_smooth(method="lm", level=.2)

###### Profile Analysis #####
#need to reference David & Davenport 2002
#check assumptions
TankMean<-FishData %>% 
  group_by(Tank, Infected, Treatment) %>% 
  summarize(meanDaysSurv=mean(DaysSurvived, na.rm=T),
            meanWeightChange=mean(WeightChange, na.rm=T),
            meanInfect=mean(InfectionDensALT, na.rm=T))

graphing<-TankMean %>% select(-meanInfect, -meanWeightChange) %>% 
  spread(Infected, meanDaysSurv)%>% gather(variable, value, -c(Tank, Treatment))
ggplot(graphing, 
       aes(x=variable, y=value, color=Treatment))+
  geom_point(aes(color=Treatment), size=2,alpha=.4, position=position_dodge(.3))+
  stat_summary(aes(color=Treatment),size=2, alpha=1, position=position_dodge(.3))+
  theme_bw()+
  scale_x_discrete("",labels=c("N"="Not Infected","Y"="Infected"))+
  scale_color_manual(values=c("Mussel"="goldenrod3","Control"="steelblue"))+
  ylab("Tank Mean Days Survived")+theme_bw()
  #geom_line(aes(group=Tank), alpha=.2)

##### Covariate data
ChlSummary
ChlFish <- ChlSummary %>% filter(Date < ymd("2018-07-02"))

ggplot(ChlFish, aes(x=Date, y=mChlA.ug.cm, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(~Compartment+Date, scales="free_x", nrow=2)+
  theme_bw()+xlab(NULL)


##### Sept4 #####
head(TankMean)
profDATA<-TankMean %>% select(-meanWeightChange, -meanInfect) %>% spread(Infected,meanDaysSurv)
#paralellism   H0: u1=u2=u3   manova
fit <- manova(cbind(N,Y) ~ Treatment, data=profDATA)
summary(fit)
#means match

#flatness   H0: (parameter.matrix)(difference.matrix)=0   hotelling T
library(Hotelling)
#comparing means of two samples (t test for multivariate data)
#null hypothesis - vectors of means are equal
hotTANK<-hotelling.test(meanDaysSurv~Infected, data=TankMean)
hotTANK
#no difference between infected & noninfected

#levels   H0: slopes are equal anova on groups
summary(aov(meanDaysSurv~Treatment, data=TankMean))
#difference between groups