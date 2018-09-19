library(readxl)
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
nitrogen<-read_xlsx("./data/WaterNutrients.xlsx", sheet="Ammonia")

library(tidyverse)
nitstand<-nitrogen %>% filter(Type=="standard")

library(ggplot2)
ggplot(nitstand, aes(x=ActualNH3.ug.L, y=x640nm))+
  geom_point()+geom_smooth(method="lm", formula=y~x)

niteq<-summary(lm(x640nm~ActualNH3.ug.L,data=nitstand))
niteq
Na=niteq$coefficients[2,1]
Nb=niteq$coefficients[1,1]

NH3raw<-nitrogen %>% mutate(PredNH3=((x640nm-Nb)/Na)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))

NH3Filtered<-NH3raw %>% filter(WaterType=="WCNF") %>% 
  select(-Type, -Location, -entry.position, -ActualNH3.ug.L) %>%
  mutate(Week=substring(Water.Sample, 10))
NH3UnFilt<-NH3raw %>% filter(WaterType=="WCN ") %>% 
  select(-Type, -Location, -entry.position, -ActualNH3.ug.L) %>%
  mutate(Week=substring(Water.Sample, 9,9))
NH3WaterNuts<- left_join(NH3Filtered, NH3UnFilt, by=c("Tank","Week")) %>%
  select(Tank,Week, PredNH3.x, PredNH3.y)
names(NH3WaterNuts)<-c("Tank","Week","FilterdNH3","UnFiltNH3")
  
phosphorus<-read_xlsx("./data/WaterNutrients.xlsx", sheet="SRP")
library(tidyverse)
phosstand<-phosphorus %>% filter(Type=="standard")%>% filter(ActualSRP.ug.L!=640)

library(ggplot2)
ggplot(phosstand, aes(x=ActualSRP.ug.L, y=x885nm))+
  geom_point()+geom_smooth(method="lm")

phoeq<-summary(lm(x885nm~ActualSRP.ug.L,data=phosstand))
phoeq
Pa=phoeq$coefficients[2,1]
Pb=phoeq$coefficients[1,1]

SRPraw<-phosphorus %>% mutate(PredSRP=((x885nm-Pb)/Pa)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))

SRPFiltered<-SRPraw %>% filter(WaterType=="WCNF") %>% 
  select(-Type, -Location, -entry.position, -ActualSRP.ug.L) %>%
  mutate(Week=substring(Water.Sample, 10))
SRPUnFilt<-SRPraw %>% filter(WaterType=="WCN ") %>% 
  select(-Type, -Location, -entry.position, -ActualSRP.ug.L) %>%
  mutate(Week=substring(Water.Sample, 9,9))
SRPWaterNuts<- left_join(SRPFiltered, SRPUnFilt, by=c("Tank","Week")) %>%
  select(Tank,Week, PredSRP.x, PredSRP.y)
names(SRPWaterNuts)<-c("Tank","Week","FilterdSRP","UnFiltSRP")

WaterNutrients<-left_join(NH3WaterNuts,SRPWaterNuts) %>% 
  mutate(Filt.element.ratio=FilterdNH3/FilterdSRP*(30.97/14.01),
         Unfilt.element.ratio=UnFiltNH3/UnFiltSRP*(30.97/14.01)) %>%
  left_join(treat[,c(1,7)])

WNmelted<-WaterNutrients %>% gather(variable, value, -c(Tank,Week, NewTreat)) %>%
  filter(Week!=4) 
WNmelted[WNmelted$variable=="Filt.element.ratio", "VarType"]="ratio"
WNmelted[WNmelted$variable=="Unfilt.element.ratio", "VarType"]="ratio"

ggplot(na.omit(WNmelted[WNmelted$VarType=="ratio" &
                        WNmelted$Tank!="T",]), 
       aes(x=Week, y=value, color=NewTreat)) + 
  geom_boxplot() + facet_grid(~variable) +theme_bw()

ggplot(WNmelted[is.na(WNmelted$VarType),], 
       aes(x=Week, y=value, color=NewTreat)) + 
  geom_boxplot() + facet_grid(~variable)+theme_bw()

lm(Unfilt.element.ratio~NewTreat, data=WaterNutrients)
