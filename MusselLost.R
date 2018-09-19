#### Mussel Communities ####
library(readxl); library(tidyverse)

MusselRaw<-read_excel("./data/CostMutData.xlsx",sheet = "MusselLengths") %>%
  add_column(BME=rep(NA,279))
MusselRaw %>% group_by(Tank, Spp) %>% summarise(n=n()) %>% spread(Spp, n)

mreg<-read_excel("./LENGTH-MASS_CLA_CCV_20161208-reg coefficients.xlsx", sheet = 2)
mreg$Spp<-c("boot","all","ACT","AMB","OREF","POCC","QVER","FFLA")

for(j in 1:nrow(MusselRaw)){
  if(is.na(match(MusselRaw$Spp, mreg$Spp)[j])) {
  MusselRaw$BME[j]<-as.numeric(mreg[2,2]*(MusselRaw$Length.mm[j]^mreg[2,5]))
}else{
  MusselRaw$BME[j]<-as.numeric(mreg[match(MusselRaw$Spp, mreg$Spp)[j],2]*
                                 (MusselRaw$Length.mm[j]^mreg[match(MusselRaw$Spp, mreg$Spp)[j],5])) 
}}

MusselComSpp <- MusselRaw %>% group_by(Tank, Spp) %>% 
  summarise(n=n(), AvgLength.mm=mean(Length.mm), AvgBM=mean(BME), sumBM=sum(BME))
MusselCom <- MusselRaw %>% group_by(Tank) %>% 
  summarise(n=n(), AvgLength.mm=round(mean(Length.mm, na.rm=T),0), 
            AvgBM=round(mean(BME, na.rm=T),1), 
            sumBM=round(sum(BME, na.rm=T),1),
            musselNDen=round(n/pi*(1.8288/2)^2,2),
            musselBMDen=round(sumBM/pi*(1.8288/2)^2,1)) %>% select(-n,-musselNDen)
FishCom<- FishData %>% group_by(Tank) %>% 
  summarise(n=n(), FAvgLength.mm=round(mean(StandLength.mm, na.rm=T),0), 
            FAvgBM=round(mean(Weight.g, na.rm=T),1), 
            FsumBM=round(sum(Weight.g, na.rm=T),1),
            FishND=round(n/pi*(1.8288/2)^2,2),
            FishBMDen=round(FsumBM/pi*(1.8288/2)^2,1)) %>% select(-n, -FishND)
bigComdata<-full_join(MusselCom, FishCom)
write.csv(bigComdata, "musselcomdata.csv")
### Mussel density = 8.25/m2, Fish density = 2.66/m2



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
