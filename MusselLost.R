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
            musselBMDen=round(sumBM/pi*(1.8288/2)^2,1))
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
summary(OrigWeights)
#k=(1/time)[ln(massfinal/massinitial)]
WCdata <- mDW %>% group_by(Mussel) %>% full_join(select(OrigWeights, -Date)) %>% 
  mutate(lastW=lag(Weight, 1),
         perOrig=round((Weight/origW)*100,1),
         perLost=round(((lastW-Weight)/origW)*100,1),
         daysSince=(Date-ymd_hms("2018-07-02 15:00:00")-5)/24,
         lostWeight=origW-Weight,
         dailLost=lastW-Weight,
         kstep=(1/as.numeric(daysSince))*(log(Weight/origW)),
         Date.day=as.Date(Date))
decayrate<-WCdata %>% filter(Date<ymd_hms("2018-08-10 00:00:00") &
                             Date>ymd_hms("2018-07-26 00:00:00")) %>%
  mutate(kdecay=(1/as.numeric(daysSince))*(log(Weight/origW))) %>%
  select(-lastW, perLost)
decayrate %>% group_by(Tank) %>% summarize(meanK=mean(kdecay, na.rm=T),
                                           meanperOrig=mean(perOrig, na.rm=T),
                                           meanperLost=100-meanperOrig)
regwc<-WCdata %>% group_by(Tank, Date.day) %>%
  summarize(meanKstep=mean(kstep, na.rm=T)) %>% 
  inner_join(physchem, by=c("Tank", "Date.day"="Date")) %>% ungroup() %>%
  mutate(TankF=factor(Tank))

decaymodel<-lme(meanKstep~Temp.C+Cond.uS+DO.mgL, random = ~1|Tank/Date.day, data=regwc)
decaynull<-lme(meanKstep~1, random=~1|Tank/Date.day, data=regwc)
anova(decaynull,decaymodel)
## NOT A BETTER MODEL THAN NULL

ggplot(WCdata, aes(x=Date, y=perOrig))+
  ylab("percent original weight")+
  geom_line(aes(group=Mussel), alpha=.2)+
  geom_smooth(size=2)+theme_bw()+
  geom_vline(xintercept=ymd_hms("2018-07-11 00:00:00"))
interval(ymd("2018-07-02"), ymd("2018-07-11")) %/% days()
  
ggplot(WCdata, aes(x=Date, y=perLost, group=Mussel))+
    ylab("percent lost daily")+geom_point()

#### Field Community Data ####
FieldM<-read_excel("CommunityDataDroughtExamples.xlsx") %>% select(-Before, -After) %>%
  gather(time, percentage, c("PerBefore", "PerAfter"))
tgplot<-ggplot(FieldM, aes(x=time, y=percentage, fill=ThemalGuild))+
  geom_bar(stat="identity")+facet_wrap(~River)+
  scale_fill_manual(name="Thermal Guild", label=c("Unknown","Sensitive","Tolerant"),
                    values=c("gold","#f8a800","#cf1c24"))+
  scale_y_continuous(label=c("0%","25%","50%","75%","100%"), name="Thermal Guild Abundance")+
  scale_x_discrete(label=c("Before","After"), name="")+fronteirstheme+
  theme(legend.position ="right",
        legend.text = element_text(size=rel(.45)),
        legend.title = element_text(size=rel(.5)),
        strip.text.x = element_text(size=rel(.6)),
        strip.background = element_rect(color="black", fill="white"))
complot<-ggplot(FieldM, aes(x=time, y=percentage, fill=Tribe))+
  geom_bar(stat="identity")+facet_wrap(~River) +
  scale_fill_manual(values=c("#005b68","#7dc3cb","#bbdddf",
                             "#d7b410","#8e8f27"))+
  scale_y_continuous(label=c("0%","25%","50%","75%","100%"), name="Community Composition")+
  scale_x_discrete(label=c("Before","After"), name="")+fronteirstheme+
  theme(legend.position ="right",
        legend.text = element_text(size=rel(.37)),
        legend.title = element_text(size=rel(.5)),
        strip.text.x = element_text(size=rel(.6)),
        strip.background = element_rect(color="black", fill="white"))
fieldplot<-plot_grid(tgplot, complot, ncol=1)
ggsave("Fig1.tiff", fieldplot, dpi=300, width=3.34, height=5)
