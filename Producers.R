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
library(car); library(lme4); library(lmerTest);library(emmeans)
Met1<-lmer(meanGPP~NewTreat + Day + (1|Tank), data=Metstats, REML=F)
anova(Met1)
summary(Met1)

Metlsd1<-emmeans(Met1, pairwise~NewTreat, adjust="tukey")
CLD(Metlsd1, alpha=.05, Letters=letters, adjust="tukey")

##### assumptions
hist(residuals(Met1),col="darkgrey") #approximates normal
plot(fitted(Met1), residuals(Met1))  #approximates heteroskodastity
qqnorm(resid(Met1));qqline(resid(Met1))

library(scales); library(ggsci)
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

col6<-c("darkgrey","#5389a6","forestgreen")
### fronteirstheme found in waterchem.R
GPP<-ggplot(Metstats, 
            aes(x=Day,y=meanGPP, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+
  ylab(expression("Gross DO Production mg "%*%cm^-2*" hr"^-1)) +
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  scale_x_continuous(breaks = unique(Metstats$Day), name="Sampling Day")+
  fronteirstheme+
  theme(axis.title.y=element_text(size=rel(.7)),
        axis.title.x=element_text(size=rel(.7)),
        axis.text.y=element_text(size=rel(.7)),
        axis.text.x=element_text(size=rel(.7)),
        legend.direction = "vertical", legend.position =c(0,.8),
        legend.text = element_text(size=rel(.65)),
        legend.title= element_text(size=rel(.7)))
ggsave("DeathFigures/Fig3.tiff", GPP, width=3.34, height= 3.34, dpi=300)
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
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/FilterVolume1)*1) %>%
  select(-Notes)%>%
  group_by(Tank, Date) %>% 
  summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Water Column")

ChlSummary<- inner_join(ChlFilter, ChlTile, by=c("Tank","Date"))%>% 
  select(-Excretion, -InfectionRound, -Notes, -nLiveMussels) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")))
names(ChlSummary)[c(3,5)]<-c("WaterColChlA.ug.mL","BenthicChlA.ug.cm")
ChlSummary<-ChlSummary %>% select(-Compartment.x, -Compartment.y) %>%
  mutate(DayF=as.factor(Day),
         WaterColChlA.ug.L=WaterColChlA.ug.mL*1000)
#completely missing chorophyl samples for Q benthos
ChlSummary[ChlSummary$Tank=="Q" & ChlSummary$Date==ymd("2018-06-17"),5:6]<-"Control"

wcchl<-ggplot(ChlSummary[ChlSummary$WaterColChlA.ug.L>0,], 
              aes(x=Day, y=WaterColChlA.ug.L, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_y_continuous(trans="log10", breaks=c(1,3,10,30,100,1000))+
  scale_x_continuous(breaks = unique(ChlSummary$Day), name="")+
  ylab(expression(atop("Water Column Chl. a", paste("ug "%*%L^-1))))+
  fronteirstheme+
  theme(legend.direction="horizontal",legend.position = c(0,1),
        axis.title.x = element_text(size=rel(0)))
benchl<-ggplot(ChlSummary[ChlSummary$BenthicChlA.ug.cm>0,], 
               aes(x=Day, y=BenthicChlA.ug.cm, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+
  ylab(expression("Gross DO Production ug "%*%L^-1)) +
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(breaks = unique(ChlSummary$Day), name="Sampling Days")+
  ylab(expression(atop("Benthic Chlorophyll a", paste("ug "%*%cm^-2))))+
  fronteirstheme
chlplot<-plot_grid(wcchl,benchl, ncol=1, labels="")
ggsave("DeathFigures/Fig4.tiff",chlplot, width=9, height=4, dpi=300)

#### Chlorophy statistics
library(car);library(lme4);library(lmerTest);library(emmeans)
#watercolumn
Wchl1<-lmer(log10(WaterColChlA.ug.L)~ NewTreat + Day + (1|Tank), data=ChlSummary, REML=F)
anova(Wchl1)
summary(Wchl1)
ranova(Wchl1)
#assumptions
hist(residuals(Wchl1),col="darkgrey") #normally distibuted?
plot(fitted(Wchl1), residuals(Wchl1)) #heteroscadastic
qqnorm(resid(Wchl1)); qqline(resid(Wchl1))

#benthic
Bchl1<-lmer(log10(BenthicChlA.ug.cm)~NewTreat + Day + (1|Tank), data=ChlSummary, REML=F)
anova(Bchl1)
summary(Bchl1)

BCHLlsd1<-emmeans(Bchl1, pairwise~NewTreat, adjust="tukey")
CLD(BCHLlsd1, alpha=.05, Letters=letters, adjust="tukey")
#assumptions
hist(residuals(Bchl1),col="darkgrey") #normally distributed?
plot(fitted(Bchl1), residuals(Bchl1)) #heteroscadastic
qqnorm(resid(Bchl1)); qqline(resid(Bchl1))

#### correlation between metabolism and chlorophyll ####
discCor<-inner_join(Chl,MetC, by=c("ChlSample"="ChlName")) %>% 
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  left_join(treat) %>%
  select(Tank, ChlSample, ChlAdensity.ug, GPPhr, NewTreat,Date)
dcMOD<-lm(ChlAdensity.ug~GPPhr, data=discCor)
summary(dcMOD)
discCor$resid<-residuals(dcMOD)

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
afdmSub<-AFDMfil %>% filter(Date==ymd("2018-07-06") | Date==ymd("2018-06-29"))

afdmplot<-ggplot(afdmSub, aes(x=Day, y=meanOrgM.gl, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  ylab(expression("Organic Matter g "%*%L^-1))+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_y_log10(breaks=c(0,.1,.25,.5,1))+
  scale_x_continuous(breaks=sort(unique(AFDMfil$Day)), name="Sampling Day")+
  fronteirstheme+
  theme(axis.title.y=element_text(size=rel(.5)),axis.text.y=element_text(size=rel(.5)),
        axis.text.x=element_text(hjust=.75, size=rel(.45)))
afdmplot
ggsave("Fig5.tiff", afdmplot, width = 5, height=3.34)


hist(AFDMfil$meanOrgM.gl)
hist(log10(AFDMfil$meanOrgM.gl))
AFDMfil <- AFDMfil %>% mutate(log10Org=log10(meanOrgM.gl))

afdm0<-glmer(meanOrgM.gl~as.factor(Day) + (1|Tank), data=AFDMfil, family=poisson)
afdm1<-glmer(meanOrgM.gl~NewTreat + Day + (1|Tank), data=AFDMfil, family=Gamma)
anova(afdm0,afdm1) #treatment DOES NOT improves the fit, model is not better than random
anova(afdm1)
summary(afdm1)
ranova(afdm1)

##### assumptions
hist(residuals(afdm1),col="darkgrey") #very not normal
plot(fitted(afdm1), residuals(afdm1))  #heteroskodastic
qqnorm(resid(afdm1))

#### correlation with chlorophyll ####
filcor<-left_join(ChlFilter, AFDMfil)
summary(lm(meanOrgM.gl~mChlA.ug.cm, data=filcor))

autoin<-Chl %>% inner_join(KeyFil, by=c("ChlSample"="WCFilter1")) %>%
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/FilterVolume1)*1) %>%
  select(-Notes) %>% left_join(AFDMraw, by=c("ChlSample"="Filter"))%>% 
  left_join(filtweight, by=c("ChlSample"="Filter"))%>%
  filter(State!="BAD"|is.na(State) | Date!=ymd("2018-07-20")) %>%
  mutate(Matter=(DryWeight-TinWeight-Pre.Weight)/FilterVolume1,
         OrgM.gml=(DryWeight-AshedWeight)/FilterVolume1,
         InOrgM.gml=(Matter-(AshedWeight-DryWeight))/FilterVolume1) %>%
  select(-Reason, -Tin) %>% group_by(Tank, Date) %>%
  summarize(AutoI=mean(ChlAdensity.ug/OrgM.gml, na.rm=T)) %>% left_join(treat) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")))
ggplot(autoin, aes(x=Date, y=AutoI, fill=NewTreat,color=NewTreat))+
  stat_summary(aes(shape=NewTreat), position=position_dodge(width=1))+
  stat_summary(geom="line") +
  geom_vline(xintercept = ymd("2018-07-02"), linetype="dashed")+
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)

auto1<-lmer(AutoI~NewTreat + Day + (1|Tank), data=autoin)
anova(auto1)

ggplot(AFDMfil, aes(x=Day, y=meanMatter.gl, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1.75))+
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  ylab(expression("Organic Matter g "%*%L^-1))+
  geom_vline(xintercept=0, linetype="dashed")

##### assumptions
hist(residuals(afdm1),col="darkgrey") #very not normal
plot(fitted(afdm1), residuals(afdm1))  #heteroskodastic
qqnorm(resid(afdm1))

##### Decomposers - cotton strips #####
cotton<-read_excel("./data/Traci_Popejoy_tensile_data_2018.xlsx")

decomp<-cotton %>% left_join(treat) %>% 
  select(Tank, NewTreat, Day.decomp, Day.postMM, Tensile.lbs) %>%
  mutate(og.lost=mean(as.matrix(cotton[cotton$Tank=="CTRL",7])) - Tensile.lbs,
         per.ratio=Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])),
         lost.str=1-Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])),
         lost.str.p=1-Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7]))*100,
         per.lost=(1-Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7]))*100)/Day.decomp)

cottonplt<-ggplot(decomp[!is.na(decomp$NewTreat),], 
                  aes(x=Day.decomp, y=Tensile.lbs, fill=NewTreat, color=NewTreat)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), position=position_dodge(width=1))+
  ylab(expression("Gross DO Production ug "%*%L^-1)) +
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))),
                    labels=c("Control", "Dead","Live"))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  scale_x_continuous(breaks=c(0,11,21,32), labels=c("initial",18,28,38),name="Sampling Day")+
  scale_y_continuous(limits=c(0,70), name="Tensile Strength (lbs)")+
  fronteirstheme+
  theme(legend.position = c(.25,.98), legend.direction = "horizontal")
ggsave("DeathFigures/Fig5.tiff", cottonplt,width=3.34, height=3, dpi=300)
mean(as.matrix(decomp[is.na(decomp$NewTreat),5]), na.rm=T)
mean_se(as.matrix(decomp[is.na(decomp$NewTreat),5]))

### Decomposition Statistics
dc1<-lmer(og.lost~NewTreat + Day.decomp + (1|Tank), 
          data=decomp[decomp$Tank!="CTRL",], REML=F)
anova(dc1)
summary(dc1)

dc1st<-emmeans(dc1, pairwise~NewTreat, adjust="tukey")
CLD(dc1st, alpha=.05, Letters=letters, adjust="tukey")

ggplot(decomp, aes(x=Day.decomp, y=og.lost, fill=NewTreat)) +
  geom_boxplot(aes(group=interaction(Day.decomp, NewTreat)))

##### assumptions
hist(residuals(dc1),col="darkgrey") #normality of residuals
plot(fitted(dc1), residuals(dc1)) #heteroscadastic
qqnorm(resid(dc1)); qqline(resid(dc1))
