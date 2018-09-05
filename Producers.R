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

Metstats<-MetC %>% group_by(Tank, Date) %>% summarize(meanNEP=mean(NEPhr),
                                                meanER=mean(ERhr),
                                                meanGPP=mean(GPPhr)) %>%
  full_join(treat)

Metgraph<-Metstats %>% select(Tank,Date,meanNEP,meanER,meanGPP) %>% 
  gather(variable, value,-Tank,-Date) %>% full_join(treat) %>% 
  select(NewTreat, Tank, Date, variable, value)

### "per unit surface area per unit time"
### GPP is mg/L per surface area per hour

library(colorspace)
CP<-sequential_hcl(3, h=c(180,70), c = 100, l = c(50, 90), power = 1)
CP<-c("black",CP[1],CP[2])

ggplot(Metgraph, aes(x=variable,y=value, color=NewTreat))+
  geom_boxplot()+
  ylab("Dissolved Oxygen mg/L per cm2 per hour")+
  xlab(NULL)+scale_color_manual(values=CP, name="Treatment")+
  scale_x_discrete(labels=c("ER","GPP","NEP"))+
  facet_wrap(~Date)+theme_bw()

MetMod<-lm(meanGPP~NewTreat, data=Metstats)
summary(MetMod)
anova(MetMod)

library(lsmeans)
leastm<-lsmeans(MetMod, "NewTreat",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

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
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  select(-Notes) %>%
  group_by(Tank, Date) %>% 
  summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Water Column")%>%
  left_join(treat) 

ChlSummary<- full_join(ChlFilter, ChlTile) %>% 
  select(-Excretion, -InfectionRound, -Notes, -nLiveMussels)

ggplot(ChlSummary, aes(x=Date, y=mChlA.ug.cm, color=NewTreat))+
  geom_boxplot() + 
  facet_grid(~Compartment+Date, space="free",scales="free")+
  theme(axis.text.x=element_text(angle = 30))

##### Ash-Free Dry Mass Analysis #####
AFDM$DryMass<-AFDM$Tin.filter.dry-AFDM$Tare
AFDM$OrganicM<-AFDM$Tin.filter.dry-AFDM$Tin.filter.combusted
AFDM$InOrganicM<-AFDM$DryMass-AFDM$OrganicM
AFDM$AFDMdensity<-AFDM$OrganicM/discA

