#### Mussel Communities ####
library(readxl); library(tidyverse)

MusselRaw<-read_excel("./data/Mussel.xlsx",sheet = "MusselLengths") %>%
  add_column(BME=rep(NA,279))
#who is in this community by tank
MusselRaw %>% group_by(Tank, Spp) %>% summarise(n=n()) %>% spread(Spp, n)

#regression from Hopper et al. 2018
mreg<-read_excel("./LENGTH-MASS_CLA_CCV_20161208-reg coefficients.xlsx", sheet = 2)
mreg$Spp<-c("boot","all","ACT","AMB","OREF","POCC","QVER","FFLA")

#apply the correct regression to each species
for(j in 1:nrow(MusselRaw)){
  if(is.na(match(MusselRaw$Spp, mreg$Spp)[j])) {
  MusselRaw$BME[j]<-as.numeric(mreg[2,2]*(MusselRaw$Length.mm[j]^mreg[2,5]))
}else{
  MusselRaw$BME[j]<-as.numeric(mreg[match(MusselRaw$Spp, mreg$Spp)[j],2]*
                                 (MusselRaw$Length.mm[j]^mreg[match(MusselRaw$Spp, mreg$Spp)[j],5])) 
}}
#get tank x species specific lengths & weights
MusselComSpp <- MusselRaw %>% group_by(Tank, Spp) %>% 
  summarise(n=n(), AvgLength.mm=mean(Length.mm, na.rm=T), 
            AvgBM=mean(BME,na.rm=T), sumBM=sum(BME,na.rm=T))
#tank averages of mass, density, length
MusselCom <- MusselRaw %>% group_by(Tank) %>% 
  summarise(n=n(), AvgLength.mm=round(mean(Length.mm, na.rm=T),0), 
            AvgBM=round(mean(BME, na.rm=T),1), 
            sumBM=round(sum(BME, na.rm=T),1),
            musselNDen=round(n/pi*(1.8288/2)^2,2),
            musselBMDen=round(sumBM/pi*(1.8288/2)^2,1))
#FishData found in FishData.R
#source('FishData.R')
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
DeathWeight<-read_excel("./data/Mussel.xlsx",sheet = "DeathWeights")

mDW<-DeathWeight %>% 
  gather(Date, Weight, -Tank, -Mussel)
mDW$Date<-as.POSIXct(as.numeric(mDW$Date)*(60*60*24), origin="1899-12-30")

library(ggplot2)
ggplot(mDW, aes(x=Date, y=Weight, group=Mussel))+
  geom_line(aes(color=Tank),size=2)+
  geom_point()+
  ylab("Weight (g)")+theme_bw()
ShellWeights<-mDW %>% group_by(Mussel) %>% filter(Weight==min(Weight))
names(ShellWeights)[4]="shellW"

OrigWeights<-mDW %>% group_by(Mussel) %>% filter(Date==min(Date))
names(OrigWeights)[4]="origW"
summary(OrigWeights)
#k=(1/time)[ln(massfinal/massinitial)]
library(lubridate)
WCdata <- mDW %>% group_by(Mussel) %>% full_join(select(OrigWeights, -Date)) %>% 
  full_join(select(ShellWeights, -Date)) %>% 
  filter(Date<ymd_hms("2018-07-20 00:00:00"), Date>ymd_hms("2018-07-02 19:00:00")) %>%
  mutate(Weightsh=Weight-shellW,
         lastW=lag(Weight,1),
         lastWsh=lag(Weightsh, 1),
         perOrig=round((Weightsh/(origW-shellW))*100,1),
         perLost=round(((lastW-Weightsh)/(origW-shellW))*100,1),
         daysSince=(Date-ymd_hms("2018-07-02 15:00:00")-5)/24,
         lostWeight=origW-Weight,
         dailLost=lastW-Weight,
         kstepW=(1/as.numeric(daysSince))*(log(Weight/origW)), #k of whole tissue
         kstepWsh=(1/as.numeric(daysSince))*(log(Weightsh/(origW-shellW))),#k of soft tissue
         Date.day=as.Date(Date))
decayrate<-WCdata %>% filter(Date<ymd_hms("2018-08-10 00:00:00") ,
                             Date>ymd_hms("2018-07-12 20:00:00")) %>%
  mutate(kdecayW=(1/as.numeric(daysSince))*(log(Weight/origW)), # k of whole tissue
         kdecayWS=(1/as.numeric(daysSince))*(log(Weightsh/(origW-shellW)))) %>% # soft tissue k
  select(-lastW)

decayrate %>% group_by(Tank, Date) %>% 
  summarize(meanKW=mean(kdecayW, na.rm=T),
            meanKWs=mean(kdecayWS, na.rm=T),
            meanperOrig=mean(perOrig, na.rm=T),
            meanperLost=100-meanperOrig)
mean(decayrate$kdecayW, na.rm=T) #decay of shell & soft tissue
mean(decayrate$kdecayWS, na.rm=T) #decay of soft tissue alone
x<-WCdata$daysSince
y1<-exp(1)^(-.01639*as.numeric(x))
y2<-exp(1)^(-.3358*as.numeric(x))
plot(x,y2, col="red", pch=16, ylab="% left", xlab="Days")#soft tissue alone 
points(x,y1, col="blue", pch=16) #shell&soft tissue, soft tissue is ~14%

ggplot(WCdata, aes(x=Date, y=perOrig))+
  ylab("percent original weight")+
  geom_line(aes(group=Mussel), alpha=.2)+
  geom_smooth(size=2)+theme_bw()+
  geom_vline(xintercept=ymd_hms("2018-07-11 00:00:00"))
