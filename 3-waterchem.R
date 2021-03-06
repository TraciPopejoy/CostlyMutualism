library(Hmisc);library(readxl); library(ggplot2); library(tidyverse)
treat<-read_excel("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx",sheet="TankData") #treatment data
nitrogen<-read_xlsx("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx", sheet="Ammonia") #ammonia absorbance
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
  dplyr::summarize(FilterdNH3ugL=mean(`WCNF`, na.rm=T), 
            UnFiltNH3ugL=mean(`WCN `, na.rm=T))

#### phosphorus readings ####
phosphorus<-read_xlsx("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx", sheet="SRP")
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
  dplyr::summarize(FilterdSRPugL=mean(`WCNF`, na.rm=T), 
            UnFiltSRPugL=mean(`WCN `, na.rm=T))
#using physchem to add dates to the table
#head(physchem) #in physiochem.R script
physchem<-read_excel("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx",sheet = "PhysioChem")
physchem[,"Date"]<-as.Date(physchem$Time)
library(lubridate)
# joining with physchem and calculating molar ratio between the samples
# this is the data frame used for graphs and models
WaterNutrients<-left_join(NH3WaterNuts,SRPWaterNuts, by=c("Tank","Week")) %>% 
  mutate(Filt.element.ratio=FilterdNH3ugL/FilterdSRPugL*(30.97/14.01),
         Unfilt.element.ratio=UnFiltNH3ugL/UnFiltSRPugL*(30.97/14.01),
         Week.c=as.numeric(paste(Week))) %>%
  left_join(physchem, by=c("Week.c"="Week", "Tank"))%>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")),
         logNH3=log10(FilterdNH3ugL)) %>%
  left_join(treat[,c(1,7)], by="Tank") %>% filter(Week!=4) %>%
  select(-Temp.C, -Cond.uS, -DO.mgL, -WaterV.mLs, -Time, -Chl1, -Chl2, -WCFilter1,
         -FilterVolume1, -WCFilter2, -FilterVolume2) %>% #removing irrelevant columns
  ungroup() %>% mutate(TankF=as.factor(Tank))
#### graphing water column nutrients & chlorophyll ####
library(cowplot)
#theme for powerpoint
ppt<-theme(axis.title.y=element_text(size=rel(1.5)),
                      axis.title.x=element_text(size=rel(1.5)),
                      axis.text.y=element_text(size=rel(1.5)),
                      axis.text.x=element_text(size=rel(1.3)),
                      legend.text = element_text(size=rel(1.1)),
                      legend.title= element_text(size=rel(1.1)))
#theme for journal
fronteirstheme<-theme(axis.title.y=element_text(size=rel(.6)),
                      axis.title.x=element_text(size=rel(.6)),
                      axis.text.y=element_text(size=rel(.6)),
                      axis.text.x=element_text(size=rel(.5)),
                      legend.direction = "vertical", legend.position =c(0,.74),
                      legend.text = element_text(size=rel(.5)),
                      legend.title= element_text(size=rel(.5)))
col6<-c("black","blue3","yellow3") #col6 in Producers.R
nh3filt<-ggplot(WaterNutrients,
                aes(x=Day, y=FilterdNH3ugL, 
                    color=NewTreat,fill=NewTreat, shape=NewTreat)) + 
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1.75),fun.args=list(mult=1))+
  scale_x_continuous(breaks = sort(unique(WaterNutrients$Day)), name="")+
  scale_color_manual(values=col6, name="Treatment", labels=c("Control","Dead Mussels","Live Mussels"))+
  scale_fill_manual(values=col6,
                  guide=guide_legend(override.aes=list(shape=c(23,22,21))), 
                  name="Treatment", labels=c("Control","Dead Mussels","Live Mussels"))+
  stat_summary(fun.y=mean,geom="point", position=position_dodge(width=1.75),
               size=2)+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab(expression("NH"[4]*"-N "*mu*"g "%*%L^-1))+
  fronteirstheme+
  theme(legend.position = c(0,.65),
        axis.title.x = element_text(size=rel(0)))
srpfil<-ggplot(WaterNutrients,
               aes(x=Day, y=FilterdSRPugL, color=NewTreat,
                   fill=NewTreat, shape=NewTreat)) + 
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1.75), fun.args=list(mult=1))+
  scale_x_continuous(breaks = sort(unique(WaterNutrients$Day)), name="")+
  scale_y_continuous(breaks=c(25,100,200,300,400,500,600))+
  scale_color_manual(values=col6, guide=F)+
  scale_fill_manual(values=col6,guide=F)+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=1.75),
               size=2)+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept = 0, linetype="dashed") +
  ylab(expression("SRP "*mu*"g "%*%L^-1))+fronteirstheme+
  theme(axis.title.x = element_text(size=rel(0)))
#wnplot<-plot_grid(nh3filt, srpfil, ncol=1,labels = NULL)

#Chlorophyll
head(ChlSummary) #found in Producers.R
wcchl<-ggplot(ChlSummary[ChlSummary$WaterColChlA.ug.L>0,], 
              aes(x=Day, y=WaterColChlA.ug.L, 
                  shape=NewTreat,fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1.75), fun.args=list(mult=1))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=1.75),
               size=2)+
  scale_fill_manual(guide=F,values=col6,
                     
                    name="Treatment", labels=c("Control","Dead Mussels","Live Mussels"))+
  scale_color_manual(values=col6, guide=F)+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_y_continuous(trans="log10", breaks=c(1,3,10,30,100,300))+
  scale_x_continuous(breaks = unique(ChlSummary$Day), name="")+
  ylab(expression("WC Chl.a "*mu*"g "%*%L^-1))+
  fronteirstheme+
  theme(axis.title.x = element_text(size=rel(0)))
benchl<-ggplot(ChlSummary[ChlSummary$BenthicChlA.ug.cm>0,], 
               aes(x=Day, y=BenthicChlA.ug.cm, 
                   shape=NewTreat, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1.75), fun.args=list(mult=1))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=1.75),
               size=2)+
  scale_fill_manual(values=col6, guide=F)+ 
  scale_color_manual(values=col6, guide=F)+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(breaks = unique(ChlSummary$Day), name="Sampling Day")+
  scale_y_continuous(breaks=c(0,2.5,5,7.5,10,12.5))+
  ylab(expression("B Chl.a "*mu*"g "%*%cm^-2))+
  fronteirstheme
#chlplot<-plot_grid(wcchl,benchl, ncol=1, labels="")
plot_grid(nh3filt, srpfil,wcchl,benchl, ncol=1)

ggsave("DeathFigures/Fig3.tiff", width=6, height=7, dpi=300)

##### Water Nutrients Statistics #####
library(lme4); library(emmeans); library(lmerTest) 
## WaterNutrients table has a row for each tank * week and 
## columns of nutrient concentrations (filtered & unfiltered)
nrow(WaterNutrients)
18*7 #one sample per tank per sampling day
#nitrogen
WNmodFnh3R<-lmer(logNH3 ~ NewTreat * Day +(1|Tank), data=WaterNutrients)
#devWNH3Tanks<-allFit(WNmodFnh3R)
# over fitting the model but what you going to do with reviewers
ranova(WNmodFnh3R) #Tank not significant
anova(WNmodFnh3R)
summary(WNmodFnh3R)

#tukey's post hoc
WNfNh<-emmeans(WNmodFnh3R, pairwise~NewTreat*Day, adjust="tukey")
CLD(WNfNh, alpha=.05, Letters=letters, adjust="tukey")
#nitrogen
hist(residuals(WNmodFnh3R),col="darkgrey") #normal distribution?
plot(fitted(WNmodFnh3R), residuals(WNmodFnh3R)) #heteroscadastic
qqnorm(resid(WNmodFnh3R)); qqline(resid(WNmodFnh3R)) #large values poorly predicted

#phosphorus
WNmodFsrp<-lmer(FilterdSRPugL ~ NewTreat * Day + (1|TankF), data=WaterNutrients)
anova(WNmodFsrp)
summary(WNmodFsrp)
#tukey's post hoc
WNsrpT<-emmeans(WNmodFsrp, pairwise~NewTreat*Day, adjust="tukey")
CLD(WNsrpT, alpha=.05, Letters=letters, adjust="tukey")

# assumptions:phosphorus
hist(residuals(WNmodFnh3),col="darkgrey") #normal distribution?
plot(fitted(WNmodFsrp), residuals(WNmodFsrp)) #heteroscadastic
qqnorm(resid(WNmodFsrp)); qqline(resid(WNmodFsrp))


WaterNutrients %>% group_by(NewTreat, Tank,Day) %>%
  filter(Day==-3|Day==39) %>%
  select(Tank, FilterdSRPugL) %>%
  spread(Day, FilterdSRPugL) %>%
  summarise_if(is.numeric, mean, is.na=F) %>%
  mutate(PerDif=((`39`-`-3`)/`-3`)*100,
         rawDif=`39`-`-3`) %>%
  group_by(NewTreat)%>% select(-Tank)%>% 
  summarise_all(mean)

##### importance of nutrient release
## This Paper
# NH3 ug m-2 - multiplied by L in tank
mean(c(100, 135, 120))*635 #75141.7
# SRP ug m-2 - multiplied by L in tank
mean(c(125,160,150))*635 #92075
#approximate volume of water in mesocosms
.35* (.76^2)*pi *1000 #in L (946 max)

## Atkinson et al. 2018 
# Nitrogen mg m-2 h-1
mean(c(1.3, 6, 5.8, 1.6, 2, 3.9, 5, 3.6))

3.65*1000 #um m-2 h-1
75141/3650
# Phosphorus mg m-2 h-1
mean(c(.3, .6, .7, .2, .3,.6, .7, .2, .6))
.47*1000 #um m-2 h-1
92075/470
## Wegner et al 2018
# P release from mussel shells to be 0.05 g m−2 y−1.
0.05*1e6/365 #ug P m-2 day-1
92075/137
