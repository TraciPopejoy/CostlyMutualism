library(readxl)
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
nitrogen<-read_xlsx("./data/WaterNutrients.xlsx", sheet="Ammonia")

library(tidyverse)
nitstand<-nitrogen %>% filter(Type=="standard")

library(ggplot2)
ggplot(nitstand, aes(x=ActualNH3.ug.L, y=x640nm))+
  geom_point()+geom_smooth(method="lm", formula=y~x)

niteq<-summary(lm(x640nm~ActualNH3.ug.L,data=nitstand))
niteq
Na=niteq$coefficients[2,1]
Nb=niteq$coefficients[1,1]

NH3raw<-nitrogen %>% mutate(PredNH3=((x640nm-Nb)/Na)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))

NH3Filtered<-NH3raw %>% filter(WaterType=="WCNF") %>% 
  select(-Type, -Location, -entry.position, -ActualNH3.ug.L) %>%
  mutate(Week=substring(Water.Sample, 10))
NH3UnFilt<-NH3raw %>% filter(WaterType=="WCN ") %>% 
  select(-Type, -Location, -entry.position, -ActualNH3.ug.L) %>%
  mutate(Week=substring(Water.Sample, 9,9))
NH3WaterNuts<- left_join(NH3Filtered, NH3UnFilt, by=c("Tank","Week")) %>%
  select(Tank,Week, PredNH3.x, PredNH3.y)
names(NH3WaterNuts)<-c("Tank","Week","FilterdNH3ugL","UnFiltNH3ugL")
  
phosphorus<-read_xlsx("./data/WaterNutrients.xlsx", sheet="SRP")
library(tidyverse)
phosstand<-phosphorus %>% filter(Type=="standard")%>% filter(ActualSRP.ug.L!=640)

library(ggplot2)
ggplot(phosstand, aes(x=ActualSRP.ug.L, y=x885nm))+
  geom_point()+geom_smooth(method="lm")

phoeq<-summary(lm(x885nm~ActualSRP.ug.L,data=phosstand))
phoeq
Pa=phoeq$coefficients[2,1]
Pb=phoeq$coefficients[1,1]

SRPraw<-phosphorus %>% mutate(PredSRPugL=((x885nm-Pb)/Pa)) %>% filter(Type=="sample") %>%
  mutate(Tank=substring(Water.Sample,1,1),
         WaterType=substring(Water.Sample,3,6))

SRPFiltered<-SRPraw %>% filter(WaterType=="WCNF") %>% 
  select(-Type, -Location, -entry.position, -ActualSRP.ug.L) %>%
  mutate(Week=substring(Water.Sample, 10))
SRPUnFilt<-SRPraw %>% filter(WaterType=="WCN ") %>% 
  select(-Type, -Location, -entry.position, -ActualSRP.ug.L) %>%
  mutate(Week=substring(Water.Sample, 9,9))
SRPWaterNuts<- left_join(SRPFiltered, SRPUnFilt, by=c("Tank","Week")) %>%
  select(Tank,Week, PredSRPugL.x, PredSRPugL.y)
names(SRPWaterNuts)<-c("Tank","Week","FilterdSRPugL","UnFiltSRPugL")

WaterNutrients<-left_join(NH3WaterNuts,SRPWaterNuts) %>% 
  mutate(Filt.element.ratio=FilterdNH3ugL/FilterdSRPugL*(30.97/14.01),
         Unfilt.element.ratio=UnFiltNH3ugL/UnFiltSRPugL*(30.97/14.01),
         Week.c=as.numeric(paste(Week))) %>%
  left_join(physchem, by=c("Week.c"="Week", "Tank"))%>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02"))) %>%
  left_join(treat[,c(1,7)])

WNmelted<-WaterNutrients %>% gather(variable, value, -c(Tank,Week, NewTreat,Date)) %>%
  mutate(Week.c=as.numeric(paste(Week)))%>%
  select(Tank, Week, NewTreat, variable, value, Date)
WNmelted[WNmelted$variable=="Filt.element.ratio", "VarType"]="ratio"
WNmelted[WNmelted$variable=="Unfilt.element.ratio", "VarType"]="ratio"

ggplot(na.omit(WNmelted[WNmelted$VarType=="ratio" &
                        WNmelted$Tank!="T",]), 
       aes(x=Week, y=value, color=NewTreat)) + 
  geom_boxplot() + facet_grid(~variable) +theme_bw()

ggplot(WNmelted[is.na(WNmelted$VarType),], 
       aes(x=Week, y=value, color=NewTreat)) + 
  geom_boxplot() + facet_grid(~variable)+theme_bw()

fronteirstheme<-theme(axis.title.y=element_text(size=rel(.6)),
                      axis.title.x=element_text(size=rel(.6)),
                      axis.text.y=element_text(size=rel(.6)),
                      axis.text.x=element_text(size=rel(.5)),
                      legend.direction = "vertical", legend.position =c(0,.74),
                      legend.text = element_text(size=rel(.5)),
                      legend.title= element_text(size=rel(.5)))

nh3filt<-ggplot(WaterNutrients, 
       aes(x=Date, y=FilterdNH3ugL, color=NewTreat)) + 
  #geom_line(aes(group=NewTreat), alpha=.5) + 
  geom_smooth(aes(group=NewTreat, fill=NewTreat), span=.6, level=.95)+
  scale_x_date(breaks = sort(unique(WaterNutrients$Date)), 
               labels = c('day -20', 'day -15', 'day -3','day 4', 
                          'day 11', 'day 17', 'day 25', 'day 39'),
               name="")+
  scale_color_grey(start=0, end=.9,name="Treatment")+
  scale_fill_grey(start=0, end=0.9, guide=F)+
  stat_summary(aes(fill=NewTreat))+
  geom_vline(xintercept = ymd("2018-07-02"), size=2)+
  ylab(expression("Filtered NH "[3]*" ug "%*%L^-1))+xlab("")+
  fronteirstheme+theme(legend.background = element_rect(fill=NA),
                       legend.position = c(0, .746))
nh3un<-ggplot(WaterNutrients, 
                aes(x=Date, y=UnFiltNH3ugL, color=NewTreat)) + 
  geom_smooth(aes(group=NewTreat, fill=NewTreat), span=.6, level=.95)+
  scale_x_date(breaks = sort(unique(WaterNutrients$Date)), 
               labels = c('day -20', 'day -15', 'day -3','day 4', 
                          'day 11', 'day 17', 'day 25', 'day 39'),
               name="")+
  scale_color_grey(start=0, end=.9, guide=F)+
  scale_fill_grey(start=0, end=0.9, guide=F)+
  stat_summary(aes(color=NewTreat))+
  geom_vline(xintercept = ymd("2018-07-02"), size=2)+
  ylab(expression("Total NH "[3]*" ug "%*%L^-1))+xlab("")+fronteirstheme
srpfil<-ggplot(WaterNutrients,
                aes(x=Date, y=FilterdSRPugL, color=NewTreat)) + 
  geom_smooth(aes(group=NewTreat, fill=NewTreat), span=.6, level=.95)+
  scale_x_date(breaks = sort(unique(WaterNutrients$Date)), 
               labels = c('day -20', 'day -15', 'day -3','day 4', 
                          'day 11', 'day 17', 'day 25', 'day 39'),
               name="")+
  scale_color_grey(start=0, end=.9, guide=F)+
  scale_fill_grey(start=0, end=0.9,guide=F)+
  stat_summary(aes(color=NewTreat))+
  geom_vline(xintercept = ymd("2018-07-02"), size=2)+
  ylab(expression("Filtered SRP ug "%*%L^-1))+xlab("")+fronteirstheme
srpun<-ggplot(WaterNutrients,
                aes(x=Date, y=UnFiltSRPugL, color=NewTreat)) + 
  geom_smooth(aes(group=NewTreat, fill=NewTreat), span=.6, level=.95)+
  scale_x_date(breaks = sort(unique(WaterNutrients$Date)), 
               labels = c('day -20\nJune 12', 'day -15\nJune 17',
                          'day -3\nJune 29',
                          'day 4\nJuly 06', 'day 11\nJuly 13', 
                          'day 17\nJuly20',
                          'day 25\nJuly 27', 'day 39\nAug 10'))+
  scale_color_grey(start=0, end=.9, guide=F)+
  scale_fill_grey(start=0, end=0.9, guide=F)+
  stat_summary(aes(color=NewTreat))+
  geom_vline(xintercept = ymd("2018-07-02"), size=2)+
  ylab(expression("Total SRP ug "%*%L^-1))+fronteirstheme
wnplot<-plot_grid(nh3filt, nh3un, srpfil,srpun, ncol=1,labels = NULL)
ggsave("Fig1.tiff", wnplot, width=3.34, height=7, dpi=300)

#### LINEAR MODELS ####
## data assumptions

## models
library(lme4); library(emmeans); library(lmerTest) 
WNmodFsrp<-lmer(FilterdSRPugL~NewTreat + Day+ (1|Tank), data=WaterNutrients)
anova(WNmodFsrp)

WNmodUsrp<-lmer(UnFiltSRPugL~NewTreat + Day +(1|Tank), data=WaterNutrients)
anova(WNmodUsrp)

WNmodFnh3<-lmer(FilterdNH3ugL~NewTreat + Day + (1|Tank), data=WaterNutrients)
anova(WNmodFnh3)
summary(WNmodFnh3)
ranova(WNmodFnh3)

WNfNh<-emmeans(WNmodFnh3, pairwise~NewTreat, adjust="tukey")
CLD(WNfNh, alpha=.05, Letters=letters, adjust="tukey")

WNmodUnh3<-lmer(UnFiltNH3ugL~NewTreat+ Day+ (1|Tank), data=WaterNutrients)
anova(WNmodUnh3)

##### assumptions
hist(residuals(MetMod),col="darkgrey") #skewed
plot(fitted(MetMod), residuals(MetMod)) 


WaterNutrients %>% filter(Week==3 | Week==2) %>% group_by(Tank,Week) %>%
  select(Tank,Week,FilterdNH3ugL) %>% mutate(id=1:n()) %>%
  spread(Week, FilterdNH3ugL) %>% 
  mutate(change=(`3`-`2`)/`3`*100) %>% left_join(treat, by="Tank") %>%
  select(Tank, NewTreat, `3`,`2`,change) %>% 
  filter(!is.na(change), Tank!="W") %>% group_by(NewTreat) %>%
  summarize(meanChange=mean(change),
            max=max(change),
            min=min(change))

library(vegan)
WCcom<-WNmelted %>% mutate(T.W=paste(Tank, Week, sep="."),
                    nob=seq(1:nrow(.)))%>% 
  group_by(T.W, NewTreat) %>%
  filter(is.na(VarType))%>%
  select(-Date, -VarType) %>% spread(variable, value)%>%
  summarize_if(is.numeric, mean, na.rm=T) %>% select(-nob) %>%
  mutate(Tank=substring(T.W, 1,1), Week=substring(T.W, 3,3))
WCcom1<-WCcom[!is.na(WCcom$UnFiltSRPugL),]
WCrda<-rda(WCcom1[,3:6]~WCcom1$NewTreat, data=WCcom1)

WCnmds<-metaMDS(!is.na(WCcom[,3:6]),k=2)
ordiplot(WCnmds,type="n")
orditorp(WCnmds,display="species",col="red",air=0.01)
orditorp(WCnmds,display="sites",cex=1.25,air=0.01)  
ordihull(WCnmds,groups=WCcom$NewTreat,draw="polygon",col="grey90",label=T)
orditorp(WCnmds,display="sites",col=c(rep(rainbow(8), 18)),air=0.01,cex=1.25)

         