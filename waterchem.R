library(readxl); library(tidyverse); library(ggplot2)
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData") #treatment data
nitrogen<-read_xlsx("./data/WaterNutrients.xlsx", sheet="Ammonia") #ammonia absorbance
#isolate nitrogen standards
nitstand<-nitrogen %>% filter(Type=="standard")

ggplot(nitstand, aes(x=ActualNH3.ug.L, y=x640nm))+
  geom_point()+geom_smooth(method="lm", formula=y~x)
#build a linear regression we use to predict NH4 concentration based on absorbance
niteq<-summary(lm(x640nm~ActualNH3.ug.L,data=nitstand))
niteq
Na=niteq$coefficients[2,1]
Nb=niteq$coefficients[1,1]

#predict ammonium from absorbance readings and delimit water type (filter/unfiltered)
NH3raw<-nitrogen %>% mutate(PredNH3=((x640nm-Nb)/Na)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))
#clean up and separate ammonium concentrations by filtered/unfiltered
#now have dataframe with tank, week, filtered NH4 and unfiltered NH4
NH3WaterNuts<-NH3raw %>% 
  select(-Type, -Location, -entry.position, -ActualNH3.ug.L) %>%
  mutate(Week=case_when(WaterType=="WCNF"~substring(Water.Sample, 10),
                        WaterType=="WCN "~substring(Water.Sample, 9,9))) %>%
  spread(WaterType, PredNH3) %>% select(-x640nm) %>% group_by(Tank, Week) %>%
  summarize(FilterdNH3ugL=mean(`WCNF`, na.rm=T), 
            UnFiltNH3ugL=mean(`WCN `, na.rm=T))

#### phosphorus readings ####
phosphorus<-read_xlsx("./data/WaterNutrients.xlsx", sheet="SRP")
#isolate the phosphorus standards
phosstand<-phosphorus %>% filter(Type=="standard")%>% filter(ActualSRP.ug.L!=640)

ggplot(phosstand, aes(x=ActualSRP.ug.L, y=x885nm))+
  geom_point()+geom_smooth(method="lm")
#build a linear model to predict SRP concentrations from absorbance
phoeq<-summary(lm(x885nm~ActualSRP.ug.L,data=phosstand))
phoeq
Pa=phoeq$coefficients[2,1]
Pb=phoeq$coefficients[1,1]
#predict SRP from absorbance readings and delimit water type (filter/unfiltered)
SRPraw<-phosphorus %>% mutate(PredSRPugL=((x885nm-Pb)/Pa)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))

#clean up and separate SRP concentrations by filtered/unfiltered
#now have dataframe with tank, week, filtered SRP and unfiltered SRP
SRPWaterNuts<-SRPraw %>%
  select(-Type, -Location, -entry.position, -ActualSRP.ug.L) %>%
  mutate(Week=case_when(WaterType=="WCNF"~substring(Water.Sample, 10),
                        WaterType=="WCN "~substring(Water.Sample, 9,9))) %>%
  spread(WaterType, PredSRPugL) %>% select(-x885nm) %>% group_by(Tank, Week) %>%
  summarize(FilterdSRPugL=mean(`WCNF`, na.rm=T), 
            UnFiltSRPugL=mean(`WCN `, na.rm=T))
#using physchem to add dates to the table
head(physchem) #in physiochem.R script
#physchem<-read_excel("./data/CostMutData.xlsx",sheet = "PhysioChem")
#physchem[,"Date"]<-as.Date(physchem$Time)
library(lubridate)
# joining with physchem and calculating molar ratio between the samples
# this is the data frame used for graphs and models
WaterNutrients<-left_join(NH3WaterNuts,SRPWaterNuts) %>% 
  mutate(Filt.element.ratio=FilterdNH3ugL/FilterdSRPugL*(30.97/14.01),
         Unfilt.element.ratio=UnFiltNH3ugL/UnFiltSRPugL*(30.97/14.01),
         Week.c=as.numeric(paste(Week))) %>%
  left_join(physchem, by=c("Week.c"="Week", "Tank"))%>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02"))) %>%
  left_join(treat[,c(1,7)]) %>% filter(Week!=4) %>%
  select(-Temp.C, -Cond.uS, -DO.mgL, -WaterV.mLs, -Time, -Chl1, -Chl2, -WCFilter1,
         -FilterVolume1, -WCFilter2, -FilterVolume2) %>% #removing irrelevant columns
  ungroup()
#### graphing water column nutrients ####
fronteirstheme<-theme(axis.title.y=element_text(size=rel(.6)),
                      axis.title.x=element_text(size=rel(.6)),
                      axis.text.y=element_text(size=rel(.6)),
                      axis.text.x=element_text(size=rel(.5)),
                      legend.direction = "vertical", legend.position =c(0,.74),
                      legend.text = element_text(size=rel(.5)),
                      legend.title= element_text(size=rel(.5)))
#col6<-c("darkgrey","#5389a6","forestgreen") #col6 in Producers.R
col6<-c("black","blue3","yellow3") #col6 in Producers.R
nh3filt<-ggplot(WaterNutrients,
                aes(x=Day, y=FilterdNH3ugL, color=NewTreat)) + 
  stat_summary(fun.y = mean, geom = "line")+
  scale_x_continuous(breaks = sort(unique(WaterNutrients$Day)), name="")+
  scale_color_manual(values=col6, name="Treatment")+
  scale_fill_manual(values=col6,
                  guide=guide_legend(override.aes=list(shape=c(23,22,21))), name="Treatment")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab(expression("NH"[3]*"-N "*mu*"g "%*%L^-1))+
  fronteirstheme+theme(legend.background = element_rect(fill=NA),
                       legend.position = c(0, .746),
                       axis.title.x=element_text(size=rel(0)))
srpfil<-ggplot(WaterNutrients,
               aes(x=Day, y=FilterdSRPugL, color=NewTreat)) + 
  stat_summary(fun.y = mean, geom = "line")+
  scale_x_continuous(breaks = sort(unique(WaterNutrients$Day)), name="Sampling Day")+
  scale_y_continuous(breaks=c(50,100,200,300,400,500,575))+
  scale_color_manual(values=col6, guide=F)+
  scale_fill_manual(values=col6,guide=F)+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab(expression("SRP "*mu*"g "%*%L^-1))+fronteirstheme
library(cowplot)
wnplot<-plot_grid(nh3filt, srpfil, ncol=1,labels = NULL)
ggsave("DeathFigures/Fig2.tiff", wnplot, width=3.34, height=5, dpi=300)

##### LINEAR MODELS #####
## models
library(lme4); library(emmeans); library(lmerTest) 
## WaterNutrients table has a row for each tank * week and 
## columns of nutrient concentrations (filtered & unfiltered)
#nitrogen
WNmodFnh3<-lmer(log10(FilterdNH3ugL)~NewTreat * Day + (1|Tank), data=WaterNutrients)
anova(WNmodFnh3)
summary(WNmodFnh3)

WNfNh<-emmeans(WNmodFnh3, pairwise~NewTreat, adjust="tukey")
CLD(WNfNh, alpha=.05, Letters=letters, adjust="tukey")

#phosphorus
WNmodFsrp<-lmer(FilterdSRPugL~NewTreat * Day + (1|Tank), data=WaterNutrients)
anova(WNmodFsrp)
summary(WNmodFsrp)

##### assumptions
#nitrogen
hist(residuals(WNmodFnh3),col="darkgrey") #normal distribution?
plot(fitted(WNmodFnh3), residuals(WNmodFnh3)) #heteroscadastic
qqnorm(resid(WNmodFnh3)); qqline(resid(WNmodFnh3)) #large values poorly predicted
#phosphorus
hist(residuals(WNmodFnh3),col="darkgrey") #normal distribution?
plot(fitted(WNmodFsrp), residuals(WNmodFsrp)) #heteroscadastic
qqnorm(resid(WNmodFsrp)); qqline(resid(WNmodFsrp))