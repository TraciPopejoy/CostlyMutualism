#### Mussel Communities ####
library(readxl); library(tidyverse)

MusselRaw<-read_excel("./data/CostMutData.xlsx",sheet = "MusselLengths")
MusselRaw %>% group_by(Tank) %>% summarise(n=n())

MusselCom <- MusselRaw %>% group_by(Tank, Spp) %>% 
  summarise(n=n(), AvgLength.mm=mean(Length.mm))



#### Mussel Tissue Decay ####
DeathWeight<-read_excel("./data/CostMutData.xlsx",sheet = "DeathWeights")

mDW<-DeathWeight %>% 
  gather(Date, Weight, -Tank, -Mussel)
mDW$Date<-as.POSIXct(as.numeric(mDW$Date)*(60*60*24), origin="1899-12-30")

library(ggplot2)
ggplot(mDW, aes(x=Date, y=Weight, group=Mussel))+
  geom_line(aes(color=Tank),size=2)+
  geom_point()+
  ylab("Weight (g)")+theme_bw()

OrigWeights<-mDW %>% group_by(Mussel) %>% filter(Date==min(Date))
names(OrigWeights)[4]="origW"

WCdata <- mDW %>% group_by(Mussel) %>% full_join(select(OrigWeights, -Date)) %>% 
  mutate(lastW=lag(Weight, 1),
         perOrig=round((Weight/origW)*100,1),
         perLost=round(((lastW-Weight)/origW)*100,1),
         daysSince=(Date-ymd_hms("2018-07-02 15:00:00")-5)/24,
         lostWeight=origW-Weight,
         dailLost=lastW-Weight)

ggplot(WCdata[WCdata$Date< ymd("2018-07-28"),], 
       aes(x=Date, y=perOrig))+
  ylab("percent original weight")+
  geom_line(aes(group=Mussel), alpha=.2)+
  geom_smooth(size=2)+theme_bw()
  
ggplot(WCdata, aes(x=Date, y=perLost, group=Mussel))+
    ylab("percent lost daily")+geom_point()

### regressions
summary(lm(dailLost~daysSince, 
           data=WCdata[WCdata$Date<ymd_hms("2018-08-01 12:00:00"),]))
summary(lm(perOrig~daysSince, 
           data=WCdata[WCdata$Date<ymd_hms("2018-08-01 12:00:00"),]))

ggplot(WCdata, aes(x=daysSince, y=perOrig))+geom_point()+
  geom_smooth(method="lm")
