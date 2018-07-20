
library(readxl)
fish<-read_excel("./data/CostMutData.xlsx", sheet="FishData")

treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
library(tidyverse)
FishData<-left_join(fish,treat, by="Tank")

DeadPerDayTreat<-FishData %>% group_by(Treatment, Infected, Died) %>% 
  summarize(n=n()) %>% filter(Died >"2018-06-18")

library(ggplot2)
ggplot(FishData, aes(x=Treatment, y=WeightChange, color=Infected))+
  geom_jitter(size=3, width=.2)+ stat_summary(color="black")+
  theme_bw()+ylab("Post-mortem weight change (grams)")

ggplot(DeadPerDayTreat, aes(x=Died, y=n, fill=Infected))+
  geom_bar(stat="identity", position="dodge")+
  ylab("Number of Dead Fish per day")+
  xlab("Day Fish Died, placed Jun 18")+
  facet_wrap(~Treatment)+theme_bw()

CumDeathTEST<-FishData %>% group_by(Treatment,Tank, Infected, Died) %>% 
  filter(Died >"2018-06-18") %>% summarise(nDied=n()) %>% 
  mutate(cumDeath=cumsum(nDied),
         Alive=10-cumDeath,
         TreatType=paste(Treatment,Infected))

ggplot(CumDeathTEST, aes(x=Died, y=cumDeath))+
  geom_smooth(aes(color=Infected), method="lm")+geom_point(position="jitter")+
  ylab("Number of Dead Fish")+
  xlab("Day Fish Died, placed Jun 18")+theme_bw()

ggplot(CumDeathTEST, aes(x=Died, y=Alive, color=TreatType))+
  geom_jitter(alpha=.4, size=3, width=.2)+
  geom_smooth(method="lm", level=.5)+
  ylab("Number of Alive Fish")+
  xlab("Day, placed Jun 18")+
  scale_y_continuous(breaks=seq(1:10))+theme_bw()

ggplot(CumDeathTEST, aes(x=Died, y=Alive, group=Tank))+
  geom_line(aes(group=Tank))+
  ylab("Number of Alive Fish")+
  xlab("Day, placed Jun 18")+
  scale_y_continuous(breaks=seq(1:10))+theme_bw()+facet_wrap(~TreatType)

FishData$TreatmentType<-paste(FishData$Treatment, FishData$Infected)
ggplot(FishData, aes(x=Died, y=Weight.g, color=TreatmentType))+
  geom_point()+geom_smooth(method="lm", level=.2)
ggplot(FishData, aes(x=Died, y=StandLength.mm, color=TreatmentType))+
  geom_point()+
  geom_smooth(method="lm", level=.2)
