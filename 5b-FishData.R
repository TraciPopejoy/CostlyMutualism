library(readxl)
fish<-read_excel("./data/Fish.xlsx", sheet="FishData")
treat<-read_excel("./data/Mesocosm_WaterQuality_MicrobialAct.xlsx",sheet="TankData")
library(tidyverse)
library(lubridate)
FishData <-left_join(fish,treat, by="Tank") %>% 
  mutate(LengthChange=DeadLength.mm - StandLength.mm,
         WeightChange=DeadWeight.g - Weight.g,
         InfectionDensALT=InfectionDensity.gloch,
         DaysSurvived=interval(ymd("2018-06-18"),ymd(DiedALT)) %/% days())
FishData$InfectionRound<-as.factor(FishData$InfectionRound)
#FishData is individual fish in rows with subsequent variables
FishData %>% group_by(alive)%>%tally()

FishData %>% 
  summarise(meanW=mean(Weight.g),sdW=sd(Weight.g),
            meanL=mean(StandLength.mm), sdL=sd(StandLength.mm))

ggplot(FishData[is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()
ggplot(FishData[!is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()

#graph of infection densitites based on infection round
ggplot(FishData[FishData$Infected=="Y",], 
       aes(x=InfectionRound, y=InfectionDensity.gloch))+geom_boxplot()+
  ylab("Glochidia per fish")+xlab("Infection round")

#Infectiondata provides summary statistics of infection density
Infectiondata<-FishData %>% filter(Infected=="Y") %>% group_by(InfectionRound) %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T),
            meanSurvival=mean(DaysSurvived, na.rm = T),
            sdSurvival=sd(DaysSurvived, na.rm=T))

FishData %>% filter(Infected=="Y") %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T))

#this adds the mean density to fish that weren't dissected for their gills
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==1,22]<-10
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==2,22]<-3
FishData[is.na(FishData$InfectionDensity.gloch) & FishData$Infected=="Y" & 
           FishData$InfectionRound==3,22]<-8
#plots survival against infection density
ggplot(FishData, aes(x=DaysSurvived, y=InfectionDensALT, color=Infected))+
  geom_point()+geom_smooth(method="lm")+ylab("Infection Density (gloch per fish)")+
  xlab("Days Survived")+scale_color_manual(values=c("Y"="red","N"="black"))
#plots histogram (?) that shows when fish died and their infection density
ggplot(FishData, aes(x=Died,y=Infected, fill=InfectionDensALT))+
  geom_bar(stat="identity")+facet_grid(~Treatment+Infected)+
  scale_fill_continuous("Infection Density",high="red")+
  theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
        title=element_text(size=11))
#plots infection density with post mortem weight gain
ggplot(FishData[FishData$Infected!="N",], aes(x=InfectionDensALT, y=WeightChange))+
  geom_point()+geom_smooth(method="lm")+facet_grid(~Treatment+Infected, space="free")+
  ylab("Weight Change (g)")+xlab("Infection Density (gloch per fish)")

#DeadPerDayTreat shows the number of fish dead in each treatment/tank per day
DeadPerDayTreat<-FishData %>% group_by(Treatment, Infected, Died) %>% 
  summarize(n=n()) %>% filter(Died >"2018-06-18")

#point graph with infection density and treatment type with weight change
ggplot(FishData, aes(x=Infected, y=WeightChange, color=InfectionDensALT))+
  geom_jitter(size=3, width=.2)+ stat_summary(color="black")+
  theme_bw()+ylab("Post-mortem weight change (grams)")+facet_wrap(~Treatment)
#graph of how many of each type of fish died per day (histogram)
ggplot(DeadPerDayTreat, aes(x=Died, y=n, fill=Infected))+
  geom_bar(stat="identity", position="dodge")+
  ylab("Number of Dead Fish per day")+
  xlab("Day Fish Died, placed Jun 18")+
  facet_wrap(~Treatment)+theme_bw()

FishData %>% filter(Died > ymd("2018-06-19 UTC")) %>%
  group_by(Treatment) %>% 
  dplyr::summarize(meanDW=mean(DeadWeight.g, na.rm=T),
                   sdDW=sd(DeadWeight.g, na.rm=T),
                   meanWC=mean(WeightChange, na.rm=T),
                   sdWC=sd(WeightChange, na.rm=T),
                   meanDL=mean(DeadLength.mm, na.rm=T),
                   sdDL=sd(DeadLength.mm, na.rm=T))

library(lmerTest)
FDW.mod<-lmer(DeadWeight.g~Treatment+(1|Tank), data=FishData)
summary(FDW.mod)
anova(FDW.mod)
#assumptions
hist(residuals(FDW.mod),col="darkgrey") #normal distribution
qqnorm(resid(FDW.mod)); qqline(resid(FDW.mod)) 

FDL.mod<-lmer(DeadLength.mm~Treatment+(1|Tank), data=FishData)
summary(FDL.mod)
anova(FDL.mod)
#assumptions
hist(residuals(FDL.mod),col="darkgrey") #normal distribution
qqnorm(resid(FDL.mod)); qqline(resid(FDL.mod)) 

FWC.mod<-lmer(WeightChange~Treatment+(1|Tank), data=FishData)
summary(FWC.mod)
anova(FWC.mod)
#assumptions
hist(residuals(FWC.mod),col="darkgrey") #normal distribution
qqnorm(resid(FWC.mod)); qqline(resid(FWC.mod)) 

### CumDeath is the tank, type of fish, and treatment type and 
### the number of cumualtive dead fish per day
CumDeath<-FishData  %>% 
  group_by(Treatment,Tank,Infected,DiedALT) %>% 
  filter(DiedALT > "2018-06-19") %>% summarise(nDied=n()) %>% 
  mutate(cumDeath=cumsum(nDied),
         Alive=5-cumDeath)
#line graph of number of dead fish in a infection type
ggplot(CumDeath, aes(x=DiedALT, y=cumDeath))+
  geom_smooth(aes(color=Infected), method="lm")+geom_point(position="jitter")+
  ylab("Number of Dead Fish")+
  xlab("Day Fish Died, placed Jun 18")+theme_bw()

#when the fish died based on TANK treatment type
ggplot(CumDeath, aes(x=DiedALT, y=Alive))+
  geom_smooth(aes(color=Treatment), method="lm")+geom_point(position="jitter")+
  ylab("Number of Alive Fish")

FishData$TreatmentType<-paste(FishData$Treatment, FishData$Infected)
ggplot(FishData, aes(x=Died, y=StandLength.mm, color=TreatmentType))+
  geom_point()+
  geom_smooth(method="lm", level=.2)

#### Invertebrate Analysis ####
library(ggsci)
#source('4-invertebrates.R')
#InvBioMass found in 4-invertebrates.R
#calculates biomass of invertebrates from entire experiment
InvBioMass %>% filter(Week==1 | Week==2) %>% group_by(Taxa) %>% 
  summarise(sumC=sum(Count),sumBM=sum(TotalBM)) %>% arrange(desc(sumC))
InvFishBMSum<-InvBioMass %>% filter(Week==1 | Week==2) %>%
  group_by(TxW, Tank, Week, Treatment, Treatment) %>% 
  summarize(TotalBiomass=sum(TotalBM, na.rm=T)/1000,
            Richness=length(unique(Taxa)),
            Count=sum(Count, na.rm=T),
            BMDensity=TotalBiomass/(3*.0103),
            CountDensity=Count/(3*.0103),
            AvgLength=mean(Length.mm, na.rm=T)) %>%
  mutate(Avg.BM=TotalBiomass/AvgLength,
         Week.n=as.numeric(Week),
         Day=case_when(Week=="1"~-1,
                       Week=="2"~11),
         Treat.F=factor(Treatment, levels=c("Control","Mussel")))
InvFishBMSum %>% group_by(Week, Treatment) %>% tally() #check we have what we need

TinvBMplot<-ggplot(InvFishBMSum, 
                   aes(x=Week,y=BMDensity,shape=Treatment,fill=Treatment,
                       group=Treatment))+
  stat_summary(fun.y = mean, geom = "line",
               position=position_dodge(width=.2),
               size=1.2,aes(color=Treatment))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=.2), fun.args=list(mult=1),
               size=1.2,aes(color=Treatment))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=.2),
               size=3)+
  scale_color_manual(values=c("black","grey"))+
  scale_fill_manual(values=c("black","grey"))+
  scale_shape_manual(values=c(21,24))+
  ylab(expression("Invert. biomass g"%*%m^-2))+
  xlab("Day")+
  scale_x_discrete(labels=c("-1","11"))+
  theme(legend.justification=c(1,1),legend.position ="none")
TinvCount<-ggplot(InvFishBMSum, 
                  aes(x=Week,y=CountDensity,fill=Treatment,shape=Treatment,
                      group=Treatment))+
  stat_summary(fun.y = mean, geom = "line",
               position=position_dodge(width=.2),
               size=1.2,aes(color=Treatment))+
  stat_summary(fun.data = mean_sdl, geom="linerange", 
               position=position_dodge(width=.2), fun.args=list(mult=1),
               size=1.2,aes(color=Treatment))+
  stat_summary(fun.y=mean, geom="point", position=position_dodge(width=.2),
               size=3)+
  scale_color_manual(values=c("black","grey"), guide=F)+
  scale_fill_manual(values=c("black","grey"), guide=F)+
  scale_shape_manual(values=c(21,24), guide=F)+
  ylab(expression("Invertebrate # m"^-2))+xlab("Day")+
  scale_x_discrete(labels=c("-1","11"))

InvBMslope<-InvFishBMSum %>% ungroup() %>%
  dplyr::select(Tank,Week, TotalBiomass) %>% spread(Week, TotalBiomass) %>%
  mutate(InvBMrate=(`2`-`1`)/12)
slopessss<-InvBMslope %>% left_join(treat)
InvSlopplot<-ggplot(slopessss, 
                    aes(x=Treatment, y=InvBMrate, color=Treatment))+
  geom_point(size=3, position=position_jitter(.25))+scale_color_jco(guide=F)+
  labs(y=expression(paste(Delta, "Biomass g")%*%day^ -1))+
  geom_hline(yintercept=0, lty=2) 
library(cowplot)
legend1<-get_legend(TinvBMplot+
                      theme(legend.position = c(.75,1.3),
                            legend.direction = "horizontal"))
fishInvPlot<-plot_grid( TinvCount,TinvBMplot, ncol=2,
                       labels = c("(a)","(b)"))
plot_grid(fishInvPlot,legend1, nrow=2, rel_heights = c(1,.1))
ggsave("FishFig1.tiff", dpi=600, width = 6, height = 3.4, units="in")


InvFishBMSum %>% 
  group_by(Treatment, Week) %>% 
  summarize(meanBM=mean(BMDensity),meanDen=mean(CountDensity))

library(lmerTest)
Bcount<-lmer(log10(Count)~Treat.F+Day+(1|Tank), data=InvFishBMSum)
summary(Bcount)
anova(Bcount)
#assumptions
hist(residuals(Bcount),col="darkgrey") #normal distribution
qqnorm(resid(Bcount)); qqline(resid(Bcount)) 

Bbio<-lmer(log10(BMDensity)~Treatment+Day+(1|Tank), data=InvFishBMSum)
summary(Bbio)
anova(Bbio)
#assumptions
hist(residuals(Bbio),col="darkgrey") #normal distribution
qqnorm(resid(Bbio)); qqline(resid(Bbio)) 

InvFishBMSum %>% group_by(Treatment, Week) %>% summarize(median(Avg.BM),
                                                         mean(Avg.BM),
                                                         median(AvgLength),
                                                         mean(CountDensity))
sizeInv<-InvBioMass %>% filter(Week==1 | Week==2, Length.mm>=1) 
sizeInv.exp<-sizeInv[rep(row.names(sizeInv), sizeInv$Count), c(2:6,14)]
bmdist<-ggplot(sizeInv.exp, aes(x=Week, y=BM))+
  #geom_violin(fill="black")+
  geom_boxplot()+
  #scale_y_continuous(trans="log")+
  facet_grid(~Treatment)
sizedist<-ggplot(sizeInv.exp, aes(x=Week, y=Length.mm))+
   # geom_violin(fill="black")+
  geom_boxplot(alpha=.5)+
  #scale_y_continuous(trans="log")+
  facet_grid(~Treatment)
library(cowplot)
plot_grid(bmdist, sizedist)
##### COX Regression #####
#need to have a list with dataframes of time, status, and variables to consider
FishSurv1<-FishData %>% filter(Died >= ymd("2018-06-20")) %>% 
  dplyr::select(Died, ID,Tank,Infected, Treatment, alive) %>%
  gather(variable, value, c(-Died, -ID,-Tank, -Treatment, -alive)) %>% 
  group_by(Died, ID) %>%
  mutate(status=n(), 
         status.char="dead",
         timeNumA=interval(ymd("2018-06-18"),ymd(Died)) %/% days()) %>% 
  dplyr::select(-variable)
FishSurv1[FishSurv1$alive=="ALIVE" &
          !is.na(FishSurv1$alive), 7]<-0

#InvBMslope found above
head(InvBMslope)
InvBModel<-tibble(Date=rep(unique(FishSurv1$Died), 18)) %>%  arrange(Date) %>% 
  add_column(Tank=rep(sort(unique(FishSurv1$Tank)),10)) %>% arrange(Tank) %>%
  mutate(timeNumA=interval(ymd("2018-06-18"),ymd(Date)) %/% days()) %>% 
  left_join(InvBMslope) %>%
  mutate(Est.I.bm=InvBMrate*timeNumA + `1`) #G is pretty bad

#bad idea to try and look at community change from 1 to 2         

FST<-left_join(FishSurv1, InvBModel) %>%  
  dplyr::select(-status.char, -alive)

library(survival)
Surv(FST$timeNumA, FST$status)

giantmodel<-coxph(Surv(timeNumA, status)~Treatment+value+`1`+Est.I.bm, data=FST)

#watertempGood found in physiochem.R
head(watertempGood)
TempModel<- watertempGood %>% filter(Date <= ymd("2018-06-30") &
                                       Date >= ymd("2018-06-18")) %>%
  select(Date, GoodTC, Tank) %>% group_by(Tank,Date) %>% filter(GoodTC >=30) %>%
  summarize(Count.30over=n()) %>% mutate(Date.x=as.POSIXct.Date(Date)+24*60*60,
                                         Cum.30over=cumsum(Count.30over)) 

FSTreduced<- FST %>% filter(Tank=="D"|Tank=="G"|Tank=="L"|Tank=="W"|Tank=="Q") %>%
  left_join(TempModel, by=c("Tank","Date"="Date.x")) %>%
  select(-Date, -Date.y) %>% filter (!is.na(status))

library(survival)
M1<-coxph(Surv(timeNumA, status)~Treatment, data=FST)
summary(M1)
M2<-coxph(Surv(timeNumA, status)~value, data=FST)
summary(M2)
M3<-coxph(Surv(timeNumA, status)~Est.I.bm, data=FST)
summary(M3)
M4<-coxph(Surv(timeNumA, status)~Treatment+value, data=FST)
summary(M4)
M5<-coxph(Surv(timeNumA, status)~Treatment+Est.I.bm, data=FST)
summary(M5)
M6<-coxph(Surv(timeNumA, status)~Est.I.bm+value, data=FST)
summary(M6)
M7<-coxph(Surv(timeNumA, status)~Treatment*value, data=FST)
summary(M7)
M8<-coxph(Surv(timeNumA, status)~Treatment*Est.I.bm, data=FST)
summary(M8)
M9<-coxph(Surv(timeNumA, status)~Est.I.bm*value, data=FST)
summary(M9)
M10<-coxph(Surv(timeNumA, status)~Treatment*Est.I.bm*value, data=FST)
summary(M10)
M0<-coxph(Surv(timeNumA, status)~1, data=FST)
summary(M0)
MX<-coxph(Surv(timeNumA, status)~Tank, data=FST)
summary(MX)

FullModels<-data.frame(BIC(MX,M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10),
                 variables=c("Tank","Null",
                             "Treatment",
                             "Infection",
                             "Inv BM",
                             "Treatment + Infection",
                             "Treatment + Inv BM",
                             "Inv BM + Infection",
                             "Treatment * Infection",
                             "Treatment * Inv BM",
                             "Inv BM * Infection",
                             "Treatment * Inv BM * Infection"))
FullModels <- FullModels %>%
  mutate(model=c("MX","M0","M1","M2","M3","M4","M5","M6","M7","M8","M9","M10"),
         deltaBIC=BIC - as.numeric(min(FullModels[FullModels$model!="M0",2])),
         rankBIC=rank(BIC)) #have to run this twice to get it to work?
View(FullModels)

loglikely<-tibble(model=c("MX","M1","M2","M3","M4","M5","M6","M7","M8","M9","M10")) %>%
  group_by(model)%>%
  mutate(loglik=round(as.numeric(get(model)$score),2),
         rsq=round(summary(get(model))$rsq[1],2),
         conc=round(summary(get(model))$concordance[1],2),
         logp.val=round(summary(get(model))$logtest[3],4))

modelTableM<-full_join(FullModels, loglikely, by="model") %>% arrange(deltaBIC) %>%
  select(model, variables, df, BIC, deltaBIC, loglik, rsq, logp.val)
library(survminer)
TsurvP<-ggsurvplot(survfit(Surv(timeNumA, status)~Treatment, data=FST),
           pval = F, conf.int = F,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1, # Change line type by groups
           palette=c('black','grey'),
           #conf.int.alpha=0.4,
           #conf.int.style="step",
           surv.median.line = "hv",  # Specify median survival
           legend.title="Treatment",
           legend.labs=c("Control","Mussel"),
           xlab="Time (days)")
IsurvP<-ggsurvplot(survfit(Surv(timeNumA, status)~value, data=FST),
           pval = F, conf.int = F,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = c(2,3), # Change line type by groups
           palette=c("darkgrey","lightgrey"),
           surv.median.line = "hv",
           legend.title="Infect. Status",
           legend.labs=c("Not Infected","Infected"),
           xlab="Time (days)",ylab="")
#library(cowplot)
survplots<-plot_grid(TsurvP$plot, IsurvP$plot, labels = c("(a)","(b)"))
ggsave("FishFigures/FishFig3.tiff",survplots, width=7, height=3.5, dpi=300)

R1<-coxph(Surv(timeNumA, status)~Treatment+value+Cum.30over, data=FSTreduced)
summary(R1)
R2<-coxph(Surv(timeNumA, status)~Cum.30over, data=FSTreduced)
summary(R2)
R3<-coxph(Surv(timeNumA, status)~Treatment+Cum.30over+Est.I.bm, data=FSTreduced)
summary(R3)
R4<-coxph(Surv(timeNumA, status)~Cum.30over+Est.I.bm, data=FSTreduced)
summary(R4)
R5<-coxph(Surv(timeNumA, status)~Count.30over, data=FSTreduced)
summary(R5)
R0<-coxph(Surv(timeNumA, status)~1, data=FSTreduced)
summary(R0)

TempR<-data.frame(BIC(R0,R1,R2,R3, R4),
                 variables=c("Null",
                             "Treatment + Infection + Temperature",
                             "Temperature",
                             "Treatment + Temperature + Inv BM",
                             "Temperature + Inv BM"))
TempR<-TempR %>%
  mutate(model=c("R0","R1","R2","R3", "R4"), 
         deltaBIC=BIC-min(TempR[TempR$model!="R0",2]))

TempLL<-tibble(model=c("R1","R2","R3","R4")) %>%  group_by(model)%>%
  mutate(loglik=round(as.numeric(get(model)$score),2),
         rsq=round(summary(get(model))$rsq[1],2),
         conc=round(summary(get(model))$concordance[1],2),
         logp.val=round(summary(get(model))$logtest[3],6))

TempModelTable<-full_join(TempR, TempLL, by="model") %>% arrange(deltaBIC) %>%
  select(model, variables, df, BIC, deltaBIC, loglik, rsq, logp.val)
mt<-rbind(modelTableM,TempModelTable)
write.csv(mt, "modeltable1.csv")

survdiff(Surv(timeNumA, status)~Treatment, data=FST)

summary(M1)
summary(M3)
summary(M5)
summary(R3)

ggsurvplot(survfit(Surv(timeNumA, status)~InvBMrate, data=FST),
           pval = TRUE, conf.int = F,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1, # Change line type by groups
           surv.median.line = "hv")

##### covariate table #####
head(physchem)
head(bigComdata)
head(ChlSummary)
head(InvFishBMSum)

covtab<-physchem %>% filter(Week==1 | Week==2) %>%
  mutate(wkc=as.character(Week)) %>% left_join(ChlSummary) %>%
  left_join(InvFishBMSum, by=c("wkc"="Week","Tank","Treatment"))%>% left_join(bigComdata) %>%
  select(Treatment, Week, Tank,Temp.C, Cond.uS, DO.mgL, WaterV.mLs, WaterColChlA.ug.mL, 
         BenthicChlA.ug.cm,BMDensity,CountDensity,Richness,AvgLength.mm,
         musselBMDen,FAvgLength.mm, FishBMDen) %>% group_by(Treatment, Week)%>%
  mutate(WaterColChlA.ug.L = WaterColChlA.ug.mL*1000) %>% select(-WaterColChlA.ug.mL)
###q week 1 doesn't have treatment
covtab[covtab$Week=="1" & covtab$Tank=="Q","Treatment"]<-"Control"
covtab1<-covtab %>% group_by(Treatment, Week) %>%
  summarize_if(is.numeric, funs(mean, sd), na.rm=T) 
write.csv(covtab1, "fishcovariate.csv")

covtabFull<-physchem %>% filter(Week==1 | Week==2) %>% left_join(ChlSummary) %>% 
  left_join(bigComdata) %>%
  select(Treatment, Week, Tank,Temp.C, Cond.uS, DO.mgL, WaterV.mLs,
         WaterColChlA.ug.L, BenthicChlA.ug.cm, AvgLength.mm,
         musselBMDen,FAvgLength.mm, FishBMDen)

summary(lmer(Temp.C~Treatment*Week+(1|Tank), data=covtabFull))#week
summary(lmer(Cond.uS~Treatment*Week+(1|Tank), data=covtabFull))#week marginal
summary(lmer(DO.mgL~Treatment*Week+(1|Tank), data=covtabFull))#week
summary(lmer(WaterColChlA.ug.L~Treatment+Week+(1|Tank), data=covtabFull))
anova(lmer(WaterColChlA.ug.L~Treatment+Week+(1|Tank), data=covtabFull))
summary(lmer(BenthicChlA.ug.cm~Treatment+Week+(1|Tank), data=covtabFull))
anova(lmer(BenthicChlA.ug.cm~Treatment+Week+(1|Tank), data=covtabFull))
#treatment, week and interaction

library(ggsci);library(cowplot)
benthic<-ggplot(covtabFull, aes(x=as.character(Week), y=BenthicChlA.ug.cm, fill=Treatment))+
  geom_boxplot()+
  scale_fill_jco(name="Tank Treatment")+
  ylab(expression("Benthic Chlorophyl a ug "%*%cm^-2))+
  xlab("Time")+
  scale_x_discrete(labels=c("start","finish"))+
  theme(legend.justification=c(1,1),legend.position = c(.45,1))
wcchl<-ggplot(covtabFull, aes(x=as.character(Week), y=WaterColChlA.ug.L, fill=Treatment))+
  geom_boxplot()+
  scale_fill_jco(guid=F)+
  ylab(expression("Water Column Chlorophyl a ug "%*%L^-1))+
  xlab("Time")+
  scale_x_discrete(labels=c("start","finish"))
plot_grid(benthic, wcchl, labels = "AUTO")



###combining invert data
fishINV<-WCinv %>% filter(Week==1 | Week==2) %>%
  mutate(TxW=paste(Tank,Week, sep=".")) %>% 
  left_join(InvFishBMSum, by=c("Treatment","TxW"))%>% 
  select(TxW, AvgDensityL.n.perL,CountDensity, Week.x, Treatment) %>% 
  gather(inv, count, -TxW,-Tank.x,-Week.x,-Treatment)
ggplot(fishINV, aes(x=Week.x, y=count, color=Treatment, group=Tank.x)) +
  geom_point() +  geom_line() + facet_wrap(~inv, scales="free")
