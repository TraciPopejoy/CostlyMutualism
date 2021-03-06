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
  purrr::map(~get_forecast_for(33.9987, -96.7197, .x)) %>% 
  purrr::map_df("hourly")
airTemp$TempC<-(airTemp$temperature-32)*.5556
airTemp$TimeCST<-airTemp$time-18000

ggplot(airTemp, aes(x=TimeCST, y=temperature))+
  geom_line()+
  geom_vline(xintercept = ymd_hms("2018-07-02 12:00:00"), 
             linetype=2)+
  geom_vline(xintercept = ymd_hms("2018-07-06 12:00:00"), 
             color="red", alpha=.3, size=2)+
  geom_vline(xintercept = ymd_hms("2018-07-13 12:00:00"), 
             color="red", alpha=.3, size=2)+
  geom_vline(xintercept = ymd_hms("2018-07-27 12:00:00"), 
             color="red", alpha=.3, size=2)
  

library(readxl)
physchem<-read_excel("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx",sheet = "PhysioChem")
physchem[,"Date"]<-as.Date(physchem$Time)
head(physchem)
physchem %>% filter(Week!=5) %>% summarize(mean(hms(Time)))
physchem %>% filter(Week==1 | Week==2) %>% group_by(Week) %>% summarize(mean(DO.mgL))
DOdiff<-physchem %>% mutate(DOdif=case_when(Week==1 ~ DO.mgL-7.42,
                                            Week==2 ~ DO.mgL-9.87)) %>%
  left_join(treat) %>%  select(Week, Tank, Time, DO.mgL, Treatment, DOdif) %>%
  filter(Week==1 | Week==2)
ggplot(DOdiff, aes(x=Treatment, y=DOdif, color=Treatment))+
  geom_point()+geom_hline(yintercept = 0)+facet_wrap(~Week)
summary(Anova(lm(DOdif~Treatment*Week, DOdiff)))

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
  ylab("Temperature degCelcius") + xlab("Date")+
  xlim(ymd_hms("2018-06-12 00:00:00"), ymd_hms("2018-08-10 08:00:00"))+
  theme_bw()
  
oxymodel<-physchem %>% left_join(treat) %>%
  select(Week, Tank, Temp.C, Cond.uS, DO.mgL, WaterV.mLs, Date, NewTreat)
  
summodel<-oxymodel %>% mutate(timeseg=case_when(Week==0 | Week==1 | Week==2 ~ "before",
                           Week==3 | Week==4 | Week==5 | Week==6 | Week==7 | Week==8  ~"after"))%>%
  group_by(Date) %>%
  summarize(Temp.C.mean=round(mean(Temp.C, na.rm=T),1),
            Temp.C.sd=round(sd(Temp.C, na.rm=T),1),
            Cond.uS.mean=round(mean(Cond.uS, na.rm=T),0),
            Cond.uS.sd=round(sd(Cond.uS, na.rm=T),0),
            DO.mgL.mean=round(mean(DO.mgL, na.rm=T),2),
            DO.mgL.sd=round(sd(DO.mgL, na.rm=T),2), 
            wvel.mLs.mean=round(mean(WaterV.mLs,na.rm=T),0),
            wvel.mLs.sd=round(sd(WaterV.mLs,na.rm=T),0))
summary(lm(DO.mgL~Temp.C, data=oxymodel[oxymodel$Date==ymd("2018-06-17") | oxymodel$Date==ymd("2018-06-29"),]))
ttt<-watertempGood %>% mutate(modelDO=0.6129*GoodTC-9.155,
                         timecor=substr(Time,10,20))%>%
  group_by(Date) %>% 
  summarize(dailymax=max(modelDO),
            dailymin=min(modelDO),
            averagetemp=mean(modelDO), 
            dailyrange=max(modelDO)-min(modelDO)) %>% left_join(summodel) %>% group_by(DO.mgL.mean) %>%
  summarise(Date=mean(Date),
            avgMax=mean(dailymax),
            avgMin=mean(dailymin),
            avgDO=mean(averagetemp),
            avgDORang=mean(dailyrange),
            MIN=min(dailymin),
            MAX=max(dailymax),
            maxRANGE=max(dailyrange)) %>% select(Date, DO.mgL.mean, avgMin, avgDO, MIN, MAX)
watertempGood %>% mutate(modelDO=0.6129*GoodTC-9.155,
                         timecor=substr(Time,10,20)) %>% filter(timecor=="12:00:00 PM") %>%
  filter(Date==ymd("2018-06-17") | Date==ymd("2018-06-29")) %>% 
  select(Date, timecor, modelDO, GoodTC)
physchem %>% filter(Tank=="L" | Tank=="D" | Tank=="Q") %>% filter(Week==1)%>% select(Week, DO.mgL, Temp.C)
tempcor<-physchem %>% mutate(hour=substr(Time,1,13)) %>%
  select(Week, Tank, Time, Temp.C, DO.mgL, Date, hour)
tempcor2<-watertempGood %>% mutate(hour=substr(paste(GoodDate),1,13))%>%
  left_join(tempcor, by=c("Tank","hour")) %>% filter(!is.na(Temp.C.y))
summary(lm(GoodTC~Temp.C.y, data=tempcor2))
TempCor<-watertempGood %>% mutate(tempMOD=GoodTC*0.7171+9.3428,
                                  modelDO=0.6129*tempMOD-9.155,
                                  timecor=substr(Time,10,20)) %>% filter(timecor=="12:00:00 PM") %>%
  filter(Date==ymd("2018-06-17") | Date==ymd("2018-06-29")) %>% 
  select(Date, timecor, modelDO, GoodTC, tempMOD)

mesotable<-physchem %>% left_join(treat) %>% group_by(Date, NewTreat) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02"))) %>% 
  summarize(Day=Day[1],
            Temp.C.m=mean(Temp.C, na.rm=T), Temp.C.sd=sd(Temp.C, na.rm=T),
            Cond.uS.m=mean(Cond.uS, na.rm=T),Cond.uS.sd=sd(Cond.uS, na.rm=T),
            DO.mgL.m=mean(DO.mgL, na.rm=T),DO.mgL.sd=sd(DO.mgL, na.rm=T),
            WV.mLs.m=mean(WaterV.mLs, na.rm=T),WV.mLs.sd=sd(WaterV.mLs, na.rm=T)) %>%
  arrange(NewTreat)
write.csv(mesotable, "physicalchem.csv")

mesotable<-physchem %>% left_join(treat) %>% group_by(Date, NewTreat) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02"))) %>% 
  summarize(Day=Day[1], mT=mean_se(Temp.C)[[1]],Tse=(mean_se(Temp.C)[[3]]-mean_se(Temp.C)[[2]])/2, 
            mC=mean_se(Cond.uS)[[1]],Cse=(mean_se(Cond.uS)[[3]]-mean_se(Cond.uS)[[2]])/2, 
            mDO=mean_se(DO.mgL)[[1]],DOse=(mean_se(DO.mgL)[[3]]-mean_se(DO.mgL)[[2]])/2, 
            mWV=mean_se(WaterV.mLs)[[1]],WVse=(mean_se(WaterV.mLs)[[3]]-mean_se(WaterV.mLs)[[2]])/2) %>%
  arrange(NewTreat)
write.csv(mesotable, "physicalchem.csv")
