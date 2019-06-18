library(readxl); library(lubridate); library(tidyverse)
##### Metabolism Analysis #####
Met<-read_excel("./data/CostMutMetabolism.xlsx",sheet = 1) #metabolism data
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData") #treatment data

discA=(pi*(27.5/2)^2)/100 #area of discs in cm2

MetC<-Met %>% mutate(NEPtime=Stime-Ftime,
                     NEPhr=(SDO-FDO)/as.numeric(NEPtime)/discA,
                     ERtime=FoTime-Ttime,
                     ERhr=(FoDO-TDO)/as.numeric(NEPtime)/discA,
                     GPPhr=NEPhr+abs(ERhr),
                     Date=date(Stime))
MetDO<- Met %>% left_join(treat, by="Tank") %>%
  mutate(Date=as.Date(Ftime))

Metstats<-MetC %>% group_by(Tank, Date) %>% dplyr::summarize(meanNEP=mean(NEPhr),
                                                meanER=mean(ERhr),
                                                meanGPP=mean(GPPhr)) %>%
  full_join(treat, by="Tank") %>%
  mutate(DatF=format(Date, format="%b-%d"),
         Day=as.numeric(Date-ymd("2018-07-02"))) %>% 
  dplyr::select(-Treatment, -nLiveMussels,-Notes, -InfectionRound, -Excretion) %>%
  ungroup()

Metgraph<-Metstats %>% dplyr::select(Tank,Date,meanNEP,meanER,meanGPP) %>% 
  gather(variable, value,-Tank,-Date) %>% full_join(treat, by="Tank") %>% 
  dplyr::select(NewTreat, Tank, Date, variable, value)

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
# Metstats dataframe has tank, date, treatment, and metabolism metrics
# 3 treatments, 4 sampling days over 18 tanks
nrow(Metstats)
Met0<-lmer(meanGPP~(1|Tank), data=Metstats, REML=F)
Met1<-lmer(meanGPP~NewTreat * Day + (1|Tank), data=Metstats, REML=F)
anova(Met1)
summary(Met1)

#tukeys post hoc test
Metlsd1<-emmeans(Met1, pairwise~NewTreat) #object contains contrasts & sig
CLD(Metlsd1, alpha=.05, Letters=letters, adjust="tukey") #letters on dif groups

##### assumptions
hist(residuals(Met1),col="darkgrey") #approximates normal
plot(fitted(Met1), residuals(Met1))  #approximates heteroskodastity
qqnorm(resid(Met1));qqline(resid(Met1))

##### graphing metabolism
library(scales); library(ggsci); library(cowplot)
col6<-c("black","blue3","yellow3")
### fronteirstheme found in waterchem.R
GPP<-ggplot(Metstats, 
            aes(x=Day,y=meanGPP, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1.75), fun.args=list(mult=1))+
  stat_summary(fun.y=mean, geom="point",
               aes( shape=NewTreat), position=position_dodge(width=1.75),
               size=4)+
  ylab(expression("Gross DO Production mg "%*%cm^-2*" hr"^-1)) +
  scale_fill_manual(values=col6, name="Treatment", 
                    guide=guide_legend(override.aes=list(shape=c(23,22,21))))+ 
  scale_color_manual(values=col6, name="Treatment")+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  scale_x_continuous(breaks = unique(Metstats$Day), name="Sampling Day")+
  theme(#axis.title.y=element_text(size=rel(.7)),
        #axis.title.x=element_text(size=rel(.7)),
        #axis.text.y=element_text(size=rel(.7)),
        #axis.text.x=element_text(size=rel(.7)),
        legend.direction = "vertical", legend.position =c(0,.8))
        #legend.text = element_text(size=rel(.65)),
        #legend.title= element_text(size=rel(.7)))

## Figure not used enough in manuscript; not including per reviewer request
#ggsave("DeathFigures/Fig3.tiff", GPP, width=3.34, height= 3.34, dpi=300)
#bottom_row <- plot_grid(ER, GPP, labels = c('B', 'C'), align = 'h')
#plot_grid(NEP, bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1.2,1))


##### Chlorophyll Analysis #####
Chl<-read_excel("./data/CostMutData.xlsx",sheet = "CHL") #chlorophyll data
physchem<-read_excel("./data/CostMutData.xlsx",sheet = "PhysioChem") %>% 
  mutate(Date=date(Time))

#building a key to link tiles to the tanks they were taken from
KeyTile<-physchem %>% dplyr::select(Date, Tank, Chl1, Chl2) %>% 
  gather(col, ChlSample, -Date, -Tank) %>% dplyr::select(-col)

# prediction chlorophyl abundance from 665 and 664 readings(to account for peophyton)
# chlorophyl abundance divided by samping area (glass fritted disc size)
ChlTile<-Chl %>% inner_join(KeyTile, by="ChlSample") %>% 
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  dplyr::select(-Notes) %>% group_by(Tank, Date) %>% 
  dplyr::summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Benthic") %>%
  left_join(treat, by="Tank")

#building a key to link filter numbers to the tanks they were taken from
KeyFila<-physchem %>% dplyr::select(Date, Tank, WCFilter1, FilterVolume1)
KeyFilb<-physchem %>% dplyr::select(Date, Tank, WCFilter2, FilterVolume2)
names(KeyFilb)[c(3,4)]<-c("WCFilter1","FilterVolume1")
KeyFil<-rbind(KeyFila, KeyFilb)

# prediction chlorophyl abundance from 665 and 664 readings(to account for peophyton)
# chlorophyl abundance divided by the volume of water filtered
ChlFilter<-Chl %>% inner_join(KeyFil, by=c("ChlSample"="WCFilter1")) %>%
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/FilterVolume1)*1) %>%
  dplyr::select(-Notes)%>%
  group_by(Tank, Date) %>% 
  dplyr::summarize(mChlA.ug.cm=mean(ChlAdensity.ug),
            Compartment="Water Column")

# joining benthic and water column samples into one table
ChlSummary<- inner_join(ChlFilter, ChlTile, by=c("Tank","Date"))%>% 
  dplyr::select(-Excretion, -InfectionRound, -Notes, -nLiveMussels) %>%
  mutate(Day=as.numeric(Date-ymd("2018-07-02")))
names(ChlSummary)[c(3,5)]<-c("WaterColChlA.ug.mL","BenthicChlA.ug.cm")
ChlSummary<-ChlSummary %>% dplyr::select(-Compartment.x, -Compartment.y) %>%
  mutate(DayF=as.factor(Day),
         WaterColChlA.ug.L=WaterColChlA.ug.mL*1000) %>%
  ungroup()
#completely missing chorophyl samples for Q benthos
ChlSummary[ChlSummary$Tank=="Q" & ChlSummary$Date==ymd("2018-06-17"),5:6]<-"Control"
# ChlSummary is used for graphing and data analysis

##### graphing chlorophyll concentration found in waterchem.R

#### Chlorophyll statistics
library(car);library(lme4);library(lmerTest);library(emmeans)
# ChlSummary has tank, date, and water column and benthic chlorophyl concentrations
nrow(ChlSummary)
18*7 #one sample per tank per date
#watercolumn
Wchl1<-lmer(log10(WaterColChlA.ug.L)~ NewTreat * Day + (1|Tank), data=ChlSummary, REML=F)
#devWChl<-allFit(Wchl1)
anova(Wchl1)
summary(Wchl1)
ranova(Wchl1)

#assumptions
hist(residuals(Wchl1),col="darkgrey") #normally distibuted?
plot(fitted(Wchl1), residuals(Wchl1)) #heteroscadastic?
qqnorm(resid(Wchl1)); qqline(resid(Wchl1))

#benthic
Bchl1<-lmer(log10(BenthicChlA.ug.cm)~NewTreat * Day + (1|Tank), data=ChlSummary, REML=F)
anova(Bchl1)
summary(Bchl1)

BCHLlsd1<-emmeans(Bchl1, pairwise~NewTreat, adjust="tukey")
CLD(BCHLlsd1, alpha=.05, Letters=letters, adjust="tukey")
#assumptions
hist(residuals(Bchl1),col="darkgrey") #normally distributed?
plot(fitted(Bchl1), residuals(Bchl1)) #heteroscadastic
qqnorm(resid(Bchl1)); qqline(resid(Bchl1))

#### correlation between metabolism and chlorophyll ####
# benthic chlorophyll concentrations and metabolism were run on 
# the same glass fritted discs. Metabolism accounts for both heterotrophic and
# photosynthetic microbes while chlorophyll only accounts for photosynthetic.
discCor<-inner_join(Chl,MetC, by=c("ChlSample"="ChlName")) %>% 
  mutate(ChlAdensity.ug=26.7*((fir664-fir750)-(sec665-sec750))*(10/discA)*1) %>%
  left_join(treat) %>%
  select(Tank, ChlSample, ChlAdensity.ug, GPPhr, NewTreat,Date)
dcMOD<-lm(ChlAdensity.ug~GPPhr, data=discCor)
summary(dcMOD)
discCor$resid<-residuals(dcMOD)

##### Ash-Free Dry Mass Analysis #####
# did not include this because correlated with other measures 
# also silicone on the tiles kept igniting and spreading ash throughout
# the muffle furnace making some of this data dubious


##### Decomposers - cotton strips #####
cotton<-read_excel("./data/Traci_Popejoy_tensile_data_2018.xlsx")

# using cotton strips to determine organic matter decomposition in the tanks
# quantification methods follow Tiegs et al 2013
decomp<-cotton %>% left_join(treat, by="Tank") %>% 
 dplyr::select(Tank, NewTreat, Day.decomp, Day.postMM, Tensile.lbs) %>%
  mutate(og.lost=mean(as.matrix(cotton[cotton$Tank=="CTRL",7])) - Tensile.lbs,
         per.ratio=Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])),
         lost.str=1-(Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7]))),
         lost.str.otr=(Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])))*100,
         lost.str.p=1-(Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])))*100,
         per.lost=abs(1-(((Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])))*100)/Day.decomp)),
         kstrip=(1/Day.decomp)*(log(Tensile.lbs/mean(as.matrix(cotton[cotton$Tank=="CTRL",7])))))
#k=(1/time)[ln(massfinal/massinitial)]
decompG<-decomp %>% filter(!is.na(decomp$NewTreat)) %>%
  mutate(NTF=factor(NewTreat))
ggplot(decompG, aes(x=Day.decomp, y=per.lost, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1), fun.args=list(mult=1))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=1), size=3)+
  scale_fill_manual(values=col6, name="Treatment", 
                    labels=c("Control", "Dead Mussels","Live Mussels"))+ 
  scale_color_manual(values=col6, name="Treatment",
                     guide=guide_legend(override.aes=list(shape=c(23,22,21))),
                     labels=c("Control", "Dead Mussels","Live Mussels"))+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  scale_x_continuous(breaks=c(11,21,32), labels=c(18,28,38),name="Sampling Day")+
  ylab("Percent Tensile Loss / Day")

cottonplt<-ggplot(decompG, 
                  aes(x=Day.decomp, y=Tensile.lbs, 
                      shape=NewTreat, fill=NewTreat, color=NewTreat)) +
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=1), fun.args=list(mult=1))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=1), size=4)+
  scale_fill_manual(values=col6, name="Treatment", 
                    labels=c("Control", "Dead Mussels","Live Mussels"))+ 
  scale_color_manual(values=col6, name="Treatment",
                     guide=guide_legend(override.aes=list(shape=c(23,22,21))),
                     labels=c("Control", "Dead Mussels","Live Mussels"))+ 
  scale_shape_manual(name = "Group", values = c(23, 22, 21), guide=F)+
  scale_x_continuous(breaks=c(11,21,32), labels=c(18,28,38),name="Sampling Day")+
  scale_y_continuous(name="Tensile Strength (lbs)")+
  geom_hline(yintercept=65.6) +fronteirstheme+
  #theme(legend.position = c(0.05,.99), legend.direction = "horizontal",
  theme(legend.position = c(0.07,.2), legend.direction = "vertical")
        axis.text.x=element_text(size=rel(.6))
ggsave("DeathFigures/Fig4.tiff", cottonplt,width=3.34, height=3, dpi=300)
mean_sdl(as.matrix(decomp[is.na(decomp$NewTreat),5]), mult=1)

### Decomposition Statistics
# decomp has tank, date of experiment, day of decomp and tensile data
# control strips were treated like other strips to account for natural variation 
# in strip tensile strength, but were not placed in any streams
# thus we excluded them from the statistics
dc1<-lmer(Tensile.lbs~ NTF + Day.decomp+(1|Tank), 
          data=decompG, REML=F)
anova(dc1)
summary(dc1)

#tukey's post hoc
dc1st<-emmeans(dc1, pairwise~NTF, adjust="tukey")
CLD(dc1st, alpha=.05, Letters=letters, adjust="tukey")
decompG %>% dplyr::group_by(NTF) %>% 
  filter(kstrip!=-Inf) %>% #tensile strength=0, can't take ln()
  dplyr::summarize(meanKstrip=mean(kstrip, na.rm=T))

##### assumptions
hist(residuals(dc1),col="darkgrey") #normality of residuals
plot(fitted(dc1), residuals(dc1)) #heteroscadastic
qqnorm(resid(dc1)); qqline(resid(dc1))
