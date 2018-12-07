library(readxl)
fish<-read_excel("./data/CostMutData.xlsx", sheet="FishData")

treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
library(tidyverse)
library(lubridate)
FishData <-left_join(fish,treat, by="Tank") %>% 
  mutate(LengthChange=DeadLength.mm - StandLength.mm,
         WeightChange=DeadWeight.g - Weight.g,
         InfectionDensALT=InfectionDensity.gloch,
         DaysSurvived=interval(ymd("2018-06-18"),ymd(FishData$DiedALT)) %/% days())
FishData$InfectionRound<-as.factor(FishData$InfectionRound)
#FishData is individual fish in rows with subsequent variables

FishData %>% group_by(Treatment) %>% 
  summarise(meanW=mean(Weight.g),sdW=sd(Weight.g),
            meanL=mean(StandLength.mm), sdL=sd(StandLength.mm))


ggplot(FishData[is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()
ggplot(FishData[!is.na(FishData$InfectionDensity.gloch),], aes(x=Died))+geom_histogram()
ggplot(FishData, aes(x=DaysSurvived, y=InfectionDensALT, color=Infected))+
  geom_point()+geom_smooth(method="lm")
ggplot(FishData, aes(x=Died))+
  geom_histogram()+facet_grid(~Treatment+Infected)+
  fungraph+theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1))
ggplot(FishData[FishData$Infected=="Y",], 
       aes(x=InfectionRound, y=InfectionDensity.gloch))+geom_boxplot()+
  ylab("Glochidia per fish")+xlab("Infection round")+fungraph

#Infectiondata provides summary statistics of infection density
Infectiondata<-FishData %>% filter(Infected=="Y") %>% group_by(InfectionRound) %>% 
  summarize(meanInfection=mean(InfectionDensity.gloch, na.rm=T),
            sdInfection=sd(InfectionDensity.gloch, na.rm=T),
            meanSurvival=mean(DaysSurvived, na.rm = T),
            sdSurvival=sd(DaysSurvived, na.rm=T))

#Infectiondata provides summary statistics of infection density
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
  xlab("Days Survived")+scale_color_manual(values=c("Y"="red","N"="black"))+fungraph
#plots histogram (?) that shows when fish died and their infection density
ggplot(FishData, aes(x=Died,y=Infected, fill=InfectionDensALT))+
  geom_bar(stat="identity")+facet_grid(~Treatment+Infected)+
  scale_fill_continuous("Infection Density",high="red")+fungraph+
  theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
        title=element_text(size=11))
#plots infection density with post mortem weight gain
ggplot(FishData[FishData$Infected!="N",], aes(x=InfectionDensALT, y=WeightChange))+
  geom_point()+geom_smooth(method="lm")+facet_grid(~Treatment+Infected, space="free")+fungraph+
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

### CumDeath is the tank, type of fish, and treatment type and 
### the number of cumualtive dead fish per day
CumDeath<-FishData  %>% 
  group_by(Treatment,Tank, DiedALT) %>% 
  filter(DiedALT > "2018-06-19") %>% summarise(nDied=n()) %>% 
  mutate(cumDeath=cumsum(nDied),
         Alive=10-cumDeath)
#line graph of number of dead fish in a treatment and infection type
ggplot(CumDeathTEST, aes(x=DiedALT, y=cumDeath))+
  geom_smooth(aes(color=Infected), method="lm")+geom_point(position="jitter")+
  ylab("Number of Dead Fish")+
  xlab("Day Fish Died, placed Jun 18")+theme_bw()

#when the fish died based on BOTH treatment type
ggplot(CumDeathTEST, aes(x=DiedALT, y=Alive, color=Treatment))+
  #geom_jitter(alpha=.4, size=3, width=.2)+
  geom_smooth(method="lm", level=.5)+
  ylab("Number of Alive Fish")+fungraph

#when the fish died based on TANK treatment type
ggplot(CumDeathTEST, aes(x=DiedALT, y=Alive, color=Treatment))+
  #geom_jitter(alpha=.4, size=3, width=.2)+
  geom_smooth(method="lm", level=.5)+
  ylab("Number of Alive Fish")+fungraph

FishData$TreatmentType<-paste(FishData$Treatment, FishData$Infected)
ggplot(FishData, aes(x=Died, y=Weight.g, color=TreatmentType))+
  geom_point()+geom_smooth(method="lm", level=.2)
ggplot(FishData, aes(x=Died, y=StandLength.mm, color=TreatmentType))+
  geom_point()+
  geom_smooth(method="lm", level=.2)

graphing<-TankMean %>% select(-meanInfect, -meanWeightChange) %>% 
  spread(Infected, meanDaysSurv)%>% gather(variable, value, -c(Tank, Treatment))
ggplot(graphing, 
       aes(x=variable, y=value, color=Treatment))+
  geom_point(aes(color=Treatment), size=2,alpha=.3, position=position_dodge(.5))+
  stat_summary(aes(color=Treatment),size=2, alpha=1, position=position_dodge(.5))+
  theme_bw()+
  scale_x_discrete("",labels=c("N"="Not Infected","Y"="Infected"))+
  scale_color_manual(values=c("Mussel"="grey","Control"="black"))+
  ylab("Average Fish Survival (days)")+fishTheme
  #geom_line(aes(group=Tank), alpha=.2)

##### COX Regression #####
#need to have a list with dataframes of time, status, and variables to consider
FishSurv1<-FishData %>% filter(Died >= ymd("2018-06-20")) %>% 
  select(Died, ID,Tank,Infected, Treatment, alive) %>%
  gather(variable, value, c(-Died, -ID,-Tank, -Treatment, -alive)) %>% 
  group_by(Died, ID) %>%
  mutate(status=n(), 
         status.char="dead",
         timeNumA=interval(ymd("2018-06-18"),ymd(Died)) %/% days()) %>% select(-variable)
FishSurv1[FishSurv1$alive=="ALIVE" &
          !is.na(FishSurv1$alive), 7]<-0

InvBMslope<-InvBMSum %>% ungroup() %>%
  filter(Week==1 | Week==2) %>%
  select(Tank,Week, TotalBiomass) %>% spread(Week, TotalBiomass) %>%
  mutate(InvBMrate=(`2`-`1`)/12)
InvBModel<-tibble(Date=rep(unique(FishSurv1$Died), 18)) %>%  arrange(Date) %>% 
  add_column(Tank=rep(sort(unique(FishSurv1$Tank)),10)) %>% arrange(Tank) %>%
  mutate(timeNumA=interval(ymd("2018-06-18"),ymd(Date)) %/% days()) %>% 
  left_join(InvBMslope) %>% 
  mutate(Est.I.bm=InvBMrate*timeNumA + `1`)
         
head(comdistance)
comModel<-tibble(Date=rep(unique(FishSurv1$Died), 18)) %>%  arrange(Date) %>% 
  add_column(Tank=rep(sort(unique(FishSurv1$Tank)),10)) %>% arrange(Tank) %>%
  mutate(timeNumA=interval(ymd("2018-06-18"),ymd(Date)) %/% days()) %>% 
  left_join(comdistance) %>% 
  mutate(Est.I.COM=InvComSlope*timeNumA) %>% select(-`1`,-`2`,-`1 n`, -`2 n`)

FST<-left_join(FishSurv1, InvBModel) %>%  
  left_join(comModel) %>%
  select(-status.char, -alive)

library(survival)
giantmodel<-coxph(Surv(timeNumA, status)~Treatment+value+`1.x`+Est.I.bm, data=FST)

TempModel<- watertempGood %>% filter(Date <= ymd("2018-06-30") &
                                       Date >= ymd("2018-06-18")) %>%
  select(Date, GoodTC, Tank) %>% group_by(Tank,Date) %>% filter(GoodTC >=30) %>%
  summarize(Count.30over=n()) %>% mutate(Date.x=as.POSIXct.Date(Date)+24*60*60,
                                         Cum.30over=cumsum(Count.30over)) 

FSTreduced<- FST %>% filter(Tank=="D" | Tank=="G" | Tank=="L" | Tank=="W" | Tank=="Q") %>%
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

FullModels<-data.frame(BIC(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10),
                 variables=c("Null",
                             "Treatment",
                             "Infection",
                             "Inv BM",
                             "Treatment + Infection",
                             "Treatment + Inv BM",
                             "Inv BM + Infection",
                             "Treatment * Infection",
                             "Treatment * Inv BM",
                             "Inv BM * Infection",
                             "Treatment * Inv BM * Infection")) %>%
  mutate(model=c("M0","M1","M2","M3","M4","M5","M6","M7","M8","M9","M10"),
         deltaBIC=BIC - min(FullModels[FullModels$model!="M0",2]))

loglikely<-tibble(model=c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10")) %>%
  group_by(model)%>%
  mutate(loglik=round(as.numeric(get(model)$score),2),
         rsq=round(summary(get(model))$rsq[1],2),
         conc=round(summary(get(model))$concordance[1],2),
         logp.val=round(summary(get(model))$logtest[3],4))

modelTableM<-full_join(FullModels, loglikely, by="model") %>% arrange(deltaBIC) %>%
  select(model, variables, df, BIC, deltaBIC, loglik, rsq, logp.val)
library(survminer)
TsurvP<-ggsurvplot(survfit(Surv(timeNumA, status)~Treatment, data=FST),
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1, # Change line type by groups
           palette="jco",
           surv.median.line = "hv",  # Specify median survival
           legend.title="Tank Treatment",
           legend.labs=c("Control","Mussel"),
           xlab="Time (days)")
IsurvP<-ggsurvplot(survfit(Surv(timeNumA, status)~value, data=FST),
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = 1, # Change line type by groups
           palette=c("black","red"),
           surv.median.line = "hv",
           legend.title="Infection Status",
           legend.labs=c("Not Infected","Infected"),
           xlab="Time (days)",
           ylab="")
BsurvP<-ggsurvplot(survfit(Surv(timeNumA, status)~Treatment+value, data=FST),
                   pval = TRUE, conf.int = TRUE,
                   risk.table = F, # Add risk table
                   risk.table.col = "strata", # Change risk table color by groups
                   linetype = 1, # Change line type by groups
                   palette=c("darkgrey","red"),
                   surv.median.line = "hv",
                   legend.title="Fish Status",
                   legend.labs=c("Control, Not Infected","Control,Infected",
                                 "Mussel, Not Infected","Mussel, Infected"),
                   xlab="",
                   legend=c(2,3))
library(cowplot)
plot_grid(TsurvP$plot, IsurvP$plot, labels = "AUTO")

modelplotLIST<-list(NHyp=survfit(M0, data=FST),
                    M8fit=survfit(M8, data=FST),
                    M5fit=survfit(M5, data=FST))
ggsurvplot_combine(modelplotLIST, FST)

g2 <- ggplotGrob(TsurvP)
g3 <- ggplotGrob(IsurvP)
min_ncol <- min(ncol(g2), ncol(g3))
g <- gridExtra::rbind.gtable(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
g$widths <- grid::unit.pmax(g2$widths, g3$widths)
grid::grid.newpage()
grid::grid.draw(g)



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
                             "Temperature + Inv BM")) %>%
  mutate(model=c("R0","R1","R2","R3", "R4"), deltaBIC=BIC-min(TempR[TempR$model!="R0",2]))

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

##### covariate table #####
head(physchem)
head(bigComdata)
head(ChlSummary)
head(InvBMSum)

covtab<-physchem %>% filter(Week==1 | Week==2) %>%
  mutate(wkc=as.character(Week)) %>% left_join(ChlSummary) %>%
  left_join(InvBMSum, by=c("wkc"="Week","Tank","Treatment","NewTreat"))%>% left_join(bigComdata) %>%
  select(Treatment, Week, Tank,Temp.C, Cond.uS, DO.mgL, WaterV.mLs, WaterColChlA.ug.mL, 
         BenthicChlA.ug.cm,BMDensity,CountDensity,Richness,AvgLength.mm,
         musselBMDen,FAvgLength.mm, FishBMDen) %>% group_by(Treatment, Week)%>%
  mutate(WaterColChlA.ug.L = WaterColChlA.ug.mL*1000) %>% select(-WaterColChlA.ug.mL) %>%
  summarize_if(is.numeric, funs(median, IQR), na.rm=T) 
###q week 1 doesn't have treatment
write.csv(covtab, "fishcovariate.csv")

covtabFull<-physchem %>% filter(Week==1 | Week==2) %>% left_join(ChlSummary) %>% left_join(bigComdata) %>%
  select(Treatment, Week, Tank,Temp.C, Cond.uS, DO.mgL, WaterV.mLs, WaterColChlA.ug.cm, BenthicChlA.ug.cm, AvgLength.mm,
         musselBMDen,FAvgLength.mm, FishBMDen)
friedman.test(y~A|B)
summary(aov(lm(Temp.C~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(Cond.uS~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(DO.mgL~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(WaterV.mLs~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(WaterColChlA.ug.L~Treatment+(1|Week), data=covtabFull))) #significant
summary(aov(lm(BenthicChlA.ug.cm~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(FAvgLength.mm~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(FishBMDen~Treatment+(1|Week), data=covtabFull)))
summary(aov(lm(BMDensity~Treatment+(1|Week), data=covtabFull)))

library(ggsci)
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
