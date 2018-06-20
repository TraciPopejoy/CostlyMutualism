library(tidyverse)

#### HOBO Temperatures #####
hobofiles<-list.files(pattern = "*.csv", path="./hobotem/")
A3<-read.csv("./hobotem/A3_0.csv", header = F)
C3<-read.csv("./hobotem/C3_0.csv", header = F)
F4<-read.csv("./hobotem/F4_0.csv", header = F)
E5<-read.csv("./hobotem/E5_0.csv", header = F)
G1<-read.csv("./hobotem/G1_0.csv", header = F)

colnames(A3)<-c("Date","Time","Temp.C","Intensity.Lux")
A3<-A3[-c(1:2),-c(5:8)]
A3$Tank<-"D"
colnames(C3)<-c("Date","Time","Temp.C","Intensity.Lux")
C3<-C3[-c(1:2),-c(5:8)]
C3$Tank<-"Q"
colnames(F4)<-c("Date","Time","Temp.C","Intensity.Lux")
F4<-F4[-c(1:2),-c(5:7)]
F4$Tank<-"L"
colnames(E5)<-c("Date","Time","Temp.C","Intensity.Lux")
E5<-E5[-c(1:2),-c(5:8)]
E5$Tank<-"G"
colnames(G1)<-c("Date","Time","Temp.C","Intensity.Lux")
G1<-G1[-c(1:2),-c(5:8)]
G1$Tank<-"W"

watertemp1<-rbind(A3,C3)
watertemp2<-rbind(watertemp1,F4)
watertemp3<-rbind(watertemp2,E5)
watertemp<-rbind(watertemp3,G1)
watertemp<-watertemp[watertemp$Temp.C!="",]

library(lubridate)
watertemp$GoodDate<-mdy_hms(watertemp[,2], tz="GMT")
watertemp$Date<-date(watertemp$GoodDate)
watertemp$GoodTC<-as.numeric(paste(watertemp$Temp.C))
head(watertemp)
ggplot(watertemp, aes(x=GoodDate, y=GoodTC))+geom_point()+facet_grid(~Tank)


library(scales)
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}

ggplot(watertemp, aes(x=GoodDate, y=GoodTC, color=Tank))+
    geom_line()+
    ylab("Water Temperature (deg C)")+xlab("Date")+
    scale_x_datetime(date_breaks="7 days", date_labels="%b %d")+
    geom_hline(yintercept = 31, 
               color="red",alpha=.3,cex=2)+
    scale_y_continuous(labels=fmt_dcimals(2))+theme_bw()

library(darksky)
airTemp<-seq(Sys.Date()-13, Sys.Date(), "1 day") %>%
map(~get_forecast_for(33.9987, -96.7197, .x)) %>% 
  map_df("hourly")
airTemp$TempC<-(airTemp$temperature-32)*.5556


library(readxl)
physchem<-read_excel("./data/CostMutData.xlsx",sheet = "PhysioChem")
head(physchem)

ggplot()+geom_line(data=watertemp, aes(x=GoodDate,y=GoodTC, color=Tank))+
  geom_line(data=airTemp, aes(x=time, y=TempC)) +
  geom_point(data=physchem, aes(x=Time, y=Temp.C))+
  xlim(ymd_hms("2018-06-16 00:00:00"), ymd_hms("2018-06-20 00:00:00"))


ggsave("SiteDepth.tiff",SiteDepth,width=7, height=4, dpi=300)
