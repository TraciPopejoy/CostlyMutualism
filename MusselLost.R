library(readxl)
DeathWeight<-read_excel("./data/CostMutData.xlsx",sheet = "DeathWeights")
library(tidyverse)


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
         perLost=round(((lastW-Weight)/origW)*100,1))

ggplot(WCdata, aes(x=Date, y=perOrig, group=Mussel))+
  ylab("percent original weight")+
  geom_line()+theme_bw()
  
ggplot(WCdata, aes(x=Date, y=perLost, group=Mussel))+
    ylab("percent lost daily")+geom_point()
