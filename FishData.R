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
       aes(x=InfectionRound, y=InfectionDensity.gloch))+geom_boxplot()+fungraph

#Infectiondata provides summary statistics of infection density
Infectiondata<-FishData %>% filter(Infected=="Y") %>% group_by(InfectionRound) %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T),
            meanSurvival=mean(DaysSurvived, na.rm = T),
            sdSurvival=sd(DaysSurvived, na.rm=T))
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
#figure out how to only do two variables?
TankMean<-FishData %>% 
  group_by(Tank, Infected, Treatment) %>% 
  summarize(meanDaysSurv=mean(DaysSurvived, na.rm=T),
            meanWeightChange=mean(WeightChange, na.rm=T),
            meanInfect=mean(InfectionDensALT, na.rm=T))

library(Hotelling)
#comparing means of two samples (t test for multivariate data)
#null hypothesis - vectors of means are equal
hotTANK<-hotelling.test(meanDaysSurv~Treatment, data=TankMean)
hotTANK
summary(aov(lm(meanDaysSurv~Treatment,data=TankMean)))

library(profileR)
### paralellism - interaction between treatment & block
### flatness -  Hotelling T test treatment effects 
### levels - anova on the groupsz
ggplot(TankMean, aes(x=Infected, y=meanWeightChange, color=meanInfect, group=Tank))+
  geom_jitter(size=3, width=.2)+ stat_summary(color="green")+
  ylab("Post-mortem weight change (grams)")+facet_wrap(~Treatment)

summary(lm(meanWeightChange~Treatment+Infected, TankMean))

ggplot(TankMean, aes(Infected, y=meanDaysSurv, color=meanInfect))+
  geom_jitter(size=3)+stat_summary(color="green") + facet_wrap(~Treatment)

summary(lm(meanDaysSurv~Treatment+Infected, TankMean))

#getting data into formats??
DaysSurvprof<-TankMean %>% select(-meanInfect, -meanWeightChange) %>% 
  spread(Infected, meanDaysSurv)%>% as.data.frame()
Weightprof<-TankMean %>% select(-meanInfect, -meanDaysSurv) %>% 
  spread(Infected, meanWeightChange)%>% as.data.frame()
pData<-cbind(DaysSurvprof, Weightprof[,3:4])
names(pData)[3:6]<-c("N.DaysSurv","Y.DaysSurv","N.WeightChange.g","Y.WeightChange.g")

## data frame for good graphs
graphing<-pData %>% gather(variable, value, -c(Tank, Treatment))
# pretty graph comparing days survived
ggplot(graphing[graphing$variable=="Y.DaysSurv" | graphing$variable=="N.DaysSurv",], 
       aes(x=variable, y=value, color=Treatment))+
  stat_summary(aes(color=Treatment),size=2, alpha=.7, position=position_dodge(.3))+
  theme_bw()+
  scale_x_discrete("",labels=c("N.DaysSurv"="Not Infected","Y.DaysSurv"="Infected"))+
  scale_color_manual(values=c("Mussel"="goldenrod3","Control"="steelblue"))+
  ylab("Tank Mean Days Survived")+fungraph
  #geom_line(aes(group=Tank), alpha=.2)
#pretty graph comparing weights
ggplot(graphing[graphing$variable=="Y.WeightChange.g" | graphing$variable=="N.WeightChange.g",], 
       aes(x=variable, y=value, color=Treatment))+
  theme_bw()+stat_summary(aes(color=Treatment),size=2, alpha=.7, position=position_dodge(.3))+
  scale_color_manual(values=c("Mussel"="goldenrod3","Control"="steelblue"))+
  scale_x_discrete("",labels=c("N.WeightChange.g"="Not Infected","Y.WeightChange.g"="Infected"))+
  ylab("Tank Mean Weight Change (g)")+fungraph
  #geom_line(aes(group=Tank), alpha=.2)

### convert to z scores- appropriate for profiles with different variables!
## uses profile analysis to look at all variables that are important (weight and survival)
library(profileR)
Daysz<-TankMean %>%as.data.frame()%>% mutate(meanDaysSurvz=scale(meanDaysSurv)) %>%
  select(-meanInfect, -meanWeightChange, -meanDaysSurv) %>% 
  spread(Infected, meanDaysSurvz)
Weightz<-TankMean %>%as.data.frame()%>% mutate(meanWeightChangez=scale(meanWeightChange)) %>%
  select(-meanInfect, -meanWeightChange, -meanDaysSurv) %>% 
  spread(Infected, meanWeightChangez)

profZZ<-cbind(Daysz, Weightz[,3:4])
names(profZZ)[3:6]<-c("N.DaysSurvZ","Y.DaysSurvZ","N.WeightChange.gZ","Y.WeightChange.gZ")

profZan<-pbg(profZZ[3:6], group=profZZ$Treatment, original.names=T, profile.plot=T)
print(profZan)
summary(profZan)
