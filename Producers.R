library(readxl); library(lubridate); library(tidyverse)
##### Metabolism Analysis #####
Met<-read_excel("./data/CostMutMetabolism.xlsx",sheet = 1)
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")

discA=(pi*(27.5/2)^2)/100 #area of discs in cm2

MetC<-Met %>% mutate(NEPtime=Stime-Ftime,
                     NEPhr=(SDO-FDO)/as.numeric(NEPtime)/discA,
                     ERtime=FoTime-Ttime,
                     ERhr=(FoDO-TDO)/as.numeric(NEPtime)/discA,
                     GPPhr=NEPhr+abs(ERhr),
                     Date=date(Stime))
library(tidyverse)
MetDO<- Met %>% left_join(treat, by="Tank") %>%
  mutate(Date=as.Date(Ftime))

Metstats<-MetC %>% group_by(Tank, Date) %>% summarize(meanNEP=mean(NEPhr),
                                                meanER=mean(ERhr),
                                                meanGPP=mean(GPPhr)) %>%
  full_join(treat, by="Tank") %>%
  mutate(DatF=format(Date, format="%b-%d"),
         Day=as.numeric(Date-ymd("2018-07-02"))) %>% 
  select(-Treatment, -nLiveMussels,-Notes, -InfectionRound, -Excretion)

Metgraph<-Metstats %>% select(Tank,Date,meanNEP,meanER,meanGPP) %>% 
  gather(variable, value,-Tank,-Date) %>% full_join(treat) %>% 
  select(NewTreat, Tank, Date, variable, value)

### "per unit surface area per unit time"
### GPP is mg/L per surface area per hour

ggplot(Metgraph, aes(x=variable,y=value, fill=NewTreat))+
  geom_boxplot()+
  ylab("Dissolved Oxygen mg/L per cm2 per hour")+
  xlab(NULL)+
  scale_x_discrete(labels=c("ER","GPP","NEP"))+
  facet_wrap(~Date)+theme_bw()
### metabolism statistics
library(car); library(lme4); library(lmerTest)
hist(Metstats$meanGPP,col="darkgrey") #skewed

Met0<-lmer(meanGPP~Day + (1|Tank), data=Metstats, REML=F)
Met1<-lmer(meanGPP~NewTreat + Day + (1|Tank), data=Metstats, REML=F)
anova(Met0,Met1) #treatment improves the fit, model is better than random
anova(Met1)
summary(Met1)
ranova(Met1)

##### assumptions
hist(residuals(Met1),col="darkgrey") #approximates normal
plot(fitted(Met1), residuals(Met1))  #approximates heteroskodastity
qqnorm(resid(Met1))

library(emmeans)
Metlsd1<-emmeans(Met1, pairwise~NewTreat, adjust="tukey")
CLD(Metlsd1, alpha=.05, Letters=letters, adjust="tukey")
ggplot(Metstats, aes(x=NewTreat, y=meanGPP))+geom_boxplot()

library(scales)
NEP<-ggplot(Metgraph[Metgraph$variable=="meanNEP",], 
            aes(x=Date,y=value, fill=NewTreat, group=interaction(Date,NewTreat)))+
  geom_boxplot()+
  ylab("Net Ecosystem Production")+
  scale_fill_aaas(name="Treatment")+
  scale_x_date(breaks = unique(Metgraph$Date), labels = date_format("%b-%d"))+
  theme(legend.position="top")
ER<-ggplot(Metgraph[Metgraph$variable=="meanER",], 
           aes(x=Date,y=value, fill=NewTreat, group=interaction(Date,NewTreat)))+
  geom_boxplot()+
  ylab("Respiration")+
  scale_fill_aaas(guide=F)+
  scale_x_date(breaks = unique(Metgraph$Date), labels = date_format("%b-%d"))+
  theme(axis.text.x=element_text(angle = 30, hjust=.7))
GPP<-ggplot(Metstats, 
            aes(x=Date,y=meanGPP, fill=NewTreat, group=interaction(Date,NewTreat)))+
  geom_boxplot()+
  ylab(expression("Gross DO Production ug "%*%L^-1)) +
  scale_fill_grey(start=0.3, end=.9, name="Treatment")+ 
  scale_x_date(breaks = unique(Metstats$Date), 
               labels = c('day 4\nJuly 06', 'day 11\nJuly 13', 
                          'day 25\nJuly 27', 'day 39\nAug 10'))+
  fronteirstheme+
  theme(axis.title.y=element_text(size=rel(.7)),
        axis.title.x=element_text(size=rel(.7)),
        axis.text.y=element_text(size=rel(.7)),
        axis.text.x=element_text(size=rel(.7)),
        legend.direction = "vertical", legend.position =c(0,.8),
        legend.text = element_text(size=rel(.65)),
        legend.title= element_text(size=rel(.7)))
ggsave("Fig2.tiff", GPP, width=3.34, height= 3.34, dpi=300)
library(cowplot) #has a default theme that I am using apparently
bottom_row <- plot_grid(ER, GPP, labels = c('B', 'C'), align = 'h')
plot_grid(NEP, bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1.2,1))

##### Chlorophyll Analysis #####
Chl<-read_excel("./data/CostMutData.xlsx",sheet = "CHL")
physchem<-read_excel("./data/CostMutData.xlsx",sheet = "PhysioChem") %>% 
  mutate(Date=date(Time))

KeyTile<-physchem %>% select(Date, Tank, Chl1, Chl2) %>% 
  gather(col, ChlSample, -Date, -Tank) %>% select(-col)

ChlTile<-Chl %>% inner_join(KeyTile) %>% 
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  select(-Notes) %>% group_by(Tank, Date) %>% 
  summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Benthic") %>%
  left_join(treat)

KeyFila<-physchem %>% select(Date, Tank, WCFilter1, FilterVolume1)
KeyFilb<-physchem %>% select(Date, Tank, WCFilter2, FilterVolume2)
names(KeyFilb)[c(3,4)]<-c("WCFilter1","FilterVolume1")
KeyFil<-rbind(KeyFila, KeyFilb)
  
ChlFilter<-Chl %>% inner_join(KeyFil, by=c("ChlSample"="WCFilter1")) %>%
  mutate(ChlAdensity.ug=case_when(
    Notes == "half dilution" ~ 26.7*((fir664-fir750)-(sec665-sec750))*(20/FilterVolume1)*1,
    Notes == "third dilution" ~ 3*26.7*((fir664-fir750)-(sec665-sec750))*(30/FilterVolume1)*1,
    TRUE ~ 26.7*((fir664-fir750)-(sec665-sec750))*(10/FilterVolume1)*1)) %>%
  select(-Notes)%>%
  group_by(Tank, Date) %>% 
  summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Water Column")

ChlSummary<- full_join(ChlFilter, ChlTile, by=c("Tank","Date")) %>% 
  select(-Excretion, -InfectionRound, -Notes, -nLiveMussels) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")))
names(ChlSummary)[c(3,5)]<-c("WaterColChlA.ug.mL","BenthicChlA.ug.cm")
ChlSummary<-ChlSummary %>% select(-Compartment.x, -Compartment.y) %>%
  mutate(DayF=as.factor(Day),
         WaterColChlA.ug.L=WaterColChlA.ug.mL*1000)
#completely missing chorophyl samples for Q benthos
ChlSummary[ChlSummary$Tank=="Q" & ChlSummary$Date==ymd("2018-06-17"),5:6]<-"Control"

wcchl<-ggplot(ChlSummary[ChlSummary$WaterColChlA.ug.L>0,], 
              aes(x=Date, y=WaterColChlA.ug.L, fill=NewTreat))+
  geom_boxplot(aes(group=interaction(Date,NewTreat))) +
  geom_vline(xintercept=ymd("2018-07-02"))+
  scale_fill_grey(start=0.3, end=.9, name="Treatment")+
  scale_y_continuous(trans="log10", breaks=c(1,3,10,30,100,1000))+
  scale_x_date(breaks = unique(ChlSummary$Date), 
               labels = c('day -20\nJune 12', 'day -15\nJune 17',
                          'day -3\nJune 29',
                          'day 4\nJuly 06', 'day 11\nJuly 13', 
                          'day 25\nJuly 27', 'day 39\nAug 10'),
               name="")+
  ylab(expression(atop("Water Column Chl. a", paste("ug "%*%L^-1))))+
  fronteirstheme+
  theme(legend.direction="horizontal",legend.position = c(0,1))
benchl<-ggplot(ChlSummary[ChlSummary$BenthicChlA.ug.cm>0,], 
               aes(x=Date, y=BenthicChlA.ug.cm, fill=NewTreat))+
  geom_boxplot(aes(group=interaction(Date,NewTreat))) +
  geom_vline(xintercept=ymd("2018-07-02")) +
  scale_fill_grey(start=0.3, end=.9,guide=F)+
  scale_x_date(breaks = unique(ChlSummary$Date), 
               labels = c('day -20\nJune 12', 'day -15\nJune 17',
                          'day -3\nJune 29',
                          'day 4\nJuly 06', 'day 11\nJuly 13', 
                          'day 25\nJuly 27', 'day 39\nAug 10'))+
  ylab(expression(atop("Benthic Chlorophyll a", paste("ug "%*%cm^-2))))+
  fronteirstheme
chlplot<-plot_grid(wcchl,benchl, ncol=1, labels="", rel_heights = c(1.1,1))
ggsave("Fig3.tiff",chlplot, width=6, height=3.34, dpi=300)

library(car);library(lme4);library(lmerTest)
#watercolumn
Wchl0<-lmer(WaterColChlA.ug.mL~DayF + (1|Tank), data=ChlSummary, REML=F)
Wchl1<-lmer(WaterColChlA.ug.mL~NewTreat+ DayF + (1|Tank), data=ChlSummary, REML=F)
anova(Wchl1)
summary(Wchl1)
ranova(Wchl1)
anova(Wchl0,Wchl1)

#benthic
Bchl0<-lmer(BenthicChlA.ug.cm~DayF + (1|Tank), data=ChlSummary, REML=F)
Bchl1<-lmer(BenthicChlA.ug.cm~NewTreat + DayF + (1|Tank), data=ChlSummary, REML=F)
anova(Bchl1)
summary(Bchl1)
ranova(Bchl1)
anova(Bchl0,Bchl1)

##### assumptions
hist(residuals(Bchl1),col="darkgrey") #skewed
plot(fitted(Bchl1), residuals(Bchl1)) 

library(emmeans)
BCHLlsd1<-emmeans(Bchl1, pairwise~NewTreat, adjust="tukey")
CLD(BCHLlsd1, alpha=.05, Letters=letters, adjust="tukey")

#### correlation between metabolism and chlorophyll ####
discCor<-inner_join(Chl,MetC, by=c("ChlSample"="ChlName")) %>% 
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  left_join(treat) %>%
  select(Tank, ChlSample, ChlAdensity.ug, GPPhr, NewTreat,Date)
dcMOD<-lm(ChlAdensity.ug~GPPhr, data=discCor)
summary(dcMOD)
discCor$resid<-residuals(dcMOD)

ggplot(discCor, aes(x=ChlAdensity.ug, y=GPPhr))+
  geom_smooth(method="lm", alpha=.5)+
  geom_point(aes(color=NewTreat), size=1.5)+
  scale_color_aaas()+theme_bw()
ggplot(discCor, aes(x=NewTreat, y=resid))+geom_boxplot()

##### Ash-Free Dry Mass Analysis #####
AFDMraw<-read_excel("./data/CostMutData.xlsx",sheet = "AFDM")
filtweight<-read_excel("./data/CostMutData.xlsx",sheet = "PreFilter")

head(KeyFil)

AFDMfil<-AFDMraw %>% inner_join(KeyFil, by=c("Filter"="WCFilter1")) %>%
  left_join(filtweight)%>% filter(State!="BAD"|is.na(State) | Date!=ymd("2018-07-20")) %>%
  mutate(Matter=(DryWeight-TinWeight-Pre.Weight)/FilterVolume1,
         OrgM.gml=(DryWeight-AshedWeight)/FilterVolume1,
         InOrgM.gml=(Matter-(AshedWeight-DryWeight))/FilterVolume1) %>%
  select(-Reason, -Tin)%>%
  group_by(Tank, Date) %>% 
  summarize(meanMatter.gl=mean(Matter, na.rm=T)*1000,
            meanOrgM.gl=mean(OrgM.gml, na.rm=T)*1000,
            meanInOrgM.gl=mean(InOrgM.gml, na.rm=T)*1000,
            Compartment="Water Column") %>%
  left_join(treat) %>% 
  select(-nLiveMussels, -Excretion, -InfectionRound, -Notes) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")))

afdmplot<-ggplot(AFDMfil, aes(x=Date, y=meanOrgM.gl, fill=NewTreat))+
  geom_boxplot(aes(group=interaction(Date, NewTreat)))+
  scale_fill_grey(start=.3, end=0.9, name="Treatment")+
  ylab(expression("Organic Matter g "%*%L^-1))+
  geom_vline(xintercept=ymd("2018-07-02"))+
  scale_y_log10(breaks=c(0,.1,.25,.5,1))+
  scale_x_date(breaks=sort(unique(AFDMfil$Date)),
               labels = c('day -20\nJune 12', 'day -15\nJune 17', 'day -3\nJune 29', 
                          'day 4\nJuly 06', 'day 11\nJuly 13', 'day 25\nJuly 27', 'day 39\nAug 10'))+
  fronteirstheme+
  theme(axis.title.y=element_text(size=rel(.5)),axis.text.y=element_text(size=rel(.5)),
        axis.text.x=element_text(hjust=.75, size=rel(.45)))
afdmplot
ggsave("Fig5.tiff", afdmplot, width = 3.34, height=2.25)

hist(AFDMfil$meanOrgM.gl)
hist(log10(AFDMfil$meanOrgM.gl))
AFDMfil <- AFDMfil %>% mutate(log10Org=log10(meanOrgM.gl))

afdm0<-lmer(log10Org~Day + (1|Tank), data=AFDMfil, REML=F)
afdm1<-lmer(log10Org~NewTreat + Day + (1|Tank), data=AFDMfil, REML=F)
anova(afdm0,afdm1) #treatment DOES NOT improves the fit, model is not better than random
anova(afdm1)
summary(afdm1)
ranova(afdm1)

##### assumptions
hist(residuals(afdm1),col="darkgrey") #very not normal
plot(fitted(afdm1), residuals(afdm1))  #heteroskodastic
qqnorm(resid(afdm1))