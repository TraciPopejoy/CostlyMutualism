library(tidyverse)

#### HOBO Temperatures #####
A3<-read.csv("./hobotemp/A3.csv", header = F)
C3<-read.csv("./hobotemp/C3.csv", header = F)
F4<-read.csv("./hobotemp/F4.csv", header = F)
E5<-read.csv("./hobotemp/E5.csv", header = F)
G1<-read.csv("./hobotemp/G1.csv", header = F)

colnames(A3)<-c("Date","Time","Temp.C","Intensity.Lux")
A3<-A3[-c(1:2),-c(5:9)]
A3$Tank<-"D"
colnames(C3)<-c("Date","Time","Temp.C","Intensity.Lux")
C3<-C3[-c(1:2),-c(5:9)]
C3$Tank<-"Q"
colnames(F4)<-c("Date","Time","Temp.C","Intensity.Lux")
F4<-F4[-c(1:2),-c(5:8)]
F4$Tank<-"L"
colnames(E5)<-c("Date","Time","Temp.C","Intensity.Lux")
E5<-E5[-c(1:2),-c(5:9)]
E5$Tank<-"G"
colnames(G1)<-c("Date","Time","Temp.C","Intensity.Lux")
G1<-G1[-c(1:2),-c(5:9)]
G1$Tank<-"W"

watertemp<-rbind(A3,C3, F4, E5, G1)
library(lubridate)
watertempGood<-watertemp %>% filter(Temp.C!="") %>% 
  mutate(GoodDate=mdy_hms(Time, tz="GMT"),
         Date=date(GoodDate),
         GoodTC=as.numeric(paste(Temp.C))) %>%
  filter(GoodDate < ymd_hms("2018-06-29 10:00:00") |
         GoodDate > ymd_hms("2018-06-30 09:00:00")) %>%
  filter(GoodDate < ymd_hms("2018-08-10 10:00:00"))

library(scales)
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}

ggplot(watertempGood, aes(x=GoodDate, y=GoodTC, color=Tank))+
    geom_line()+
    ylab("Water Temperature (deg C)")+xlab("Date")+
    scale_x_datetime(date_breaks="7 days", date_labels="%b %d")+
    geom_hline(yintercept = 35, 
               color="red",alpha=.3,cex=2)+
  xlim(ymd_hms("2018-06-18 00:00:00"), ymd_hms("2018-08-10 08:00:00"))+
    scale_y_continuous(labels=fmt_dcimals(2))+theme_bw()

library(darksky)
airTemp<-seq(ymd("2018-06-18"), ymd("2018-08-10"), "1 day") %>%
map(~get_forecast_for(33.9987, -96.7197, .x)) %>% 
  map_df("hourly")
airTemp$TempC<-(airTemp$temperature-32)*.5556
airTemp$TimeCST<-airTemp$time-18000


library(readxl)
physchem<-read_excel("./data/CostMutData.xlsx",sheet = "PhysioChem")
physchem[,"Date"]<-as.Date(physchem$Time)
head(physchem)
physchem %>% filter(Week==1 | Week==2) %>% group_by(Week) %>% summarize(mean(DO.mgL))
DOdiff<-physchem %>% mutate(DOdif=case_when(Week==1 ~ DO.mgL-7.42,
                                            Week==2 ~ DO.mgL-9.87)) %>%
  left_join(treat) %>%  select(Week, Tank, Time, DO.mgL, Treatment, DOdif) %>%
  filter(Week==1 | Week==2)
ggplot(DOdiff, aes(x=Treatment, y=DOdif, color=Treatment))+
  geom_point()+geom_hline(yintercept = 0)+facet_wrap(~Week)
summary(Anova(lm(DOdif~Treatment*Week, DOdiff)))

histTemp<-read_excel("./data/CostMutData.xlsx",sheet = "HistWeath") #accuweather
histTemp$TempC<-(histTemp$HistTemps-32)*.5556

ggplot()+
  geom_line(data=airTemp, aes(x=TimeCST, y=TempC), size=1) +
  geom_line(data=watertempGood[watertempGood$Tank=="L",],
            aes(x=GoodDate,y=GoodTC, group=Tank),color="blue",
            size=1.1)+
  geom_hline(yintercept = 35, color="red", size=1.1)+
  geom_point(data=physchem, aes(x=Time, y=Temp.C), color="blue")+
  #geom_point(aes(x=ymd_hm("2018-06-28 14:00"), y=37), color="black")+
  #geom_point(data=histTemp, aes(x=Date, y=TempC), color="red", size=3)+
  ylab(expression("Temperature" *~degree*C)) + xlab("Date")+
  ylim(22,37)+
  xlim(ymd_hms("2018-06-18 00:00:00"), ymd_hms("2018-06-29 08:00:00"))+theme_bw()
ggsave("Temperature.tiff",SiteDepth,width=7, height=4, dpi=300)

#### covariate table ####
head(physchem)


physchem$Week<-as.factor(physchem$Week)
pcgraph<-physchem %>% left_join(treat) %>% mutate(WaterV.cfs=WaterV.mLs*3.53147e-5)
ggplot(pcgraph, aes(x=Date, y=WaterV.cfs, color=Treatment))+
  geom_point(size=2, alpha=.6, position = position_jitter(width=1))+fungraph
physCOV<-physchem %>% left_join(treat) %>% 
  filter(Time < ymd("2018-06-30"), !is.na(Week)) %>% 
  group_by(Treatment, Week) %>% 
  summarize(meanDO=round(mean(DO.mgL, na.rm=T),2),
            meanCond=round(mean(Cond.uS, na.rm=T),0),
            meanTemp=round(mean(Temp.C),1),
            meanWV=round(mean(WaterV.mLs, na.rm=T),0))
tempDailyData<-watertempGood %>% group_by(Date) %>% 
  summarize(dailymax=max(GoodTC),
            dailymin=min(GoodTC),
            averagetemp=mean(GoodTC), 
            dailyrange=max(GoodTC)-min(GoodTC))
tempStats<-tempDailyData %>% filter(Date > ymd("2018-06-18") &
                                    Date < ymd("2018-06-29")) %>%
  summarise(avgMax=mean(dailymax),
            avgMin=mean(dailymin),
            avgTemp=mean(averagetemp),
            avgTempRang=mean(dailyrange),
            MIN=min(dailymin),
            MAX=max(dailymax),
            maxRANGE=max(dailyrange))
library(tidyverse)
library(dataRetrieval)
nwisQW <- readNWISqw(c("07338500","07337900", "07336200"),
                     "00010", 
                     startDate = "1960-01-01", 
                     endDate = "2015-12-30") %>%
  select(site_no, sample_dt, result_va)
library(lubridate)
View(nwisQW %>% mutate(Month=months(sample_dt),
                  Year=year(sample_dt)) %>% 
  filter(Month=="June" | Month=="July" | Month=="August") %>%
  group_by(site_no, Year) %>%
  summarize(count=n(),
            avgT=mean(result_va, na.omit=T),
            maxT=max(result_va, na.omit=T)))

                    

### week date table
unique(physchem$Week)
WeekDate<-physchem %>% select(Week, Time) %>% group_by(Week) %>%
            slice(1) %>% mutate(SamplingDate=date(Time))

ggplot()+
  geom_line(data=watertempGood[watertempGood$Tank=="L",], 
            aes(x=GoodDate,y=GoodTC),color="blue") +
  geom_vline(xintercept=ymd_hms("2018-06-12 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-06-17 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-06-29 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-07-06 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-07-13 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-07-27 12:00:00"), size=3, alpha=.4)+
  geom_vline(xintercept=ymd_hms("2018-08-10 12:00:00"), size=3, alpha=.4)+
  geom_text()
  ylab("Temperature degCelcius") + xlab("Date")+
  xlim(ymd_hms("2018-06-12 00:00:00"), ymd_hms("2018-08-10 08:00:00"))+
  theme_bw()
