library(reshape2); library(ggplot2); library(tidyverse); library(vegan)
library(readxl)

#############     Invert Data in Benthic Samples    ##############
#bring in enclosure and treatment data, trait/taxa list, and length weight regressions
head(treat)
BiomassReg<-read_excel("./data/Macroinv Power Law Coeffs TBP.xlsx")
#THIS IS THE INSECT DATA
Inv<-read_excel("./data/CostMutData.xlsx", sheet="BenthicMacroInvert",
                col_types=c("guess","text","guess","numeric","numeric",
                            "guess","guess")) %>%
  select(-Notes) %>%
  left_join(treat) %>% select(-nLiveMussels, -Notes, -InfectionRound, -Excretion) %>% 
  mutate(TxW=as.factor(paste(Tank,Week, sep=".")))
#clean the data frame
sort(unique(Inv$Taxa)) #check to make sure no misspellings

invkey<-read_excel("./data/Macroinv Power Law Coeffs TBP.xlsx", sheet="Key")
InvClean<-Inv %>% right_join(invkey) %>% filter(Macro.Meio=="macro")
#how many individuals of each taxa in each enclosure/time; long format
RawCounts<-InvClean %>% group_by(TxW, Treatment, NewTreat, Tank, Week) %>% 
  summarise(TotalInv=sum(Count))

ggplot(RawCounts[RawCounts$Week==1 | 
                 RawCounts$Week==2,], 
       aes(x=Week, y=TotalInv, fill=Treatment))+
  geom_boxplot()+facet_wrap(~Week, scale="free")+theme_bw()
RawCounts %>% group_by(Week) %>% summarize(n())

ggplot(RawCounts, aes(x=Week, y=TotalInv, fill=NewTreat))+
  geom_boxplot()+facet_wrap(~Week, scale="free")+theme_bw()

#############    Invert Biomass Calculation     #############
#apply appropriate biomass regressions to each length

biomass<-function(Fam, Ord, Length) {
  tmass<-NULL
  #if not ID'd to family, use order level regressions
  if(is.na(match(Fam, BiomassReg$Family))){ 
    plcoeffs<-BiomassReg[BiomassReg$Order == Ord &
                           !is.na(BiomassReg$Order == Ord),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==Fam&
                           BiomassReg$Order == Ord&
                           !is.na(BiomassReg$Family==Fam), ] }
  #which regressions were actually built for insects this size
  #idx2<- c(Length>=plcoeffs[,19] & Length<=plcoeffs[,20]) 
  #idx2[is.na(idx2)]<-T #if no size range listed, use the regression anyways
  d1<-plcoeffs[]
  mass<-d1$a*(Length^d1$b) #power law to determine biomass
  tmass<-cbind(tmass,mass) #get vector of mass values from different eq.
}
InvBioMass<-InvClean %>% select(-Notes,-ID.person) %>%
  group_by(TxW, Taxa, Family, Length.mm) %>% 
  mutate(BM=mean(biomass(Family,Order, Length.mm)),
         TotalBM=Count*BM)

InvBMSum<-InvBioMass %>% group_by(TxW, Tank, Week, NewTreat, Treatment) %>% 
    summarize(TotalBiomass=sum(TotalBM, na.omit=T)/1000,
              Richness=length(unique(Taxa)),
              Count=sum(Count, na.omit=T),
              BMDensity=TotalBiomass/(3*.0103),
              CountDensity=Count/(3*.0103),
              AvgLength=mean(Length.mm, na.rm=T)) %>%
  mutate(Avg.BM=TotalBiomass/AvgLength,
         Week.n=as.numeric(Week))

ggplot(InvBMSum, aes(x=Week, y=TotalBiomass, fill=NewTreat))+
  geom_boxplot()+facet_wrap(~Week, scale="free")+theme_bw()

invgraph<-InvBMSum %>% left_join(physchem, by=c("Week.n"="Week","Tank")) %>%
  filter(Week.n>.5) %>%
  select(-Temp.C, -Cond.uS, -DO.mgL, -WaterV.mLs, -Chl1, -Chl2)

library(lubridate)
Invplot<-ggplot(invgraph, aes(x=Date,y=BMDensity, 
                     fill=NewTreat, group=interaction(Date,NewTreat)))+
  geom_boxplot()+geom_vline(xintercept=ymd("2018-07-02"))+
  ylab(expression("Invertebrate Biomass mg "%*%m^-2)) +
  scale_fill_grey(start=0.3, end=0.9, name="Treatment")+fronteirstheme
ggsave("Fig4.tiff",Invplot, width=3.34, height=3.34, dpi=300)


###FISH
TinvBMplot<-ggplot(InvBMSum[InvBMSum$Week==1 |
                InvBMSum$Week==2, ], 
       aes(x=Week, y=BMDensity, fill=Treatment))+
  geom_boxplot()+
  scale_fill_jco(name="Tank Treatment")+
  ylab(expression("Invertebrate Biomass g"%*%m^-2))+
  xlab("Time")+
  scale_x_discrete(labels=c("start","finish"))+
  theme(legend.justification=c(1,1),legend.position = c(.75,1))
TinvCount<-ggplot(InvBMSum[InvBMSum$Week==1 |
                              InvBMSum$Week==2,], 
                   aes(x=Week, y=CountDensity, fill=Treatment))+
  geom_boxplot()+scale_fill_jco(guide=F)+
  ylab(expression("Invertebrates m"^-2))+xlab("Time")+
  scale_x_discrete(labels=c("start","finish"))
slopessss<-InvBMslope %>% left_join(treat)
InvSlopplot<-ggplot(slopessss, 
                  aes(x=Treatment, y=InvBMrate, color=Treatment))+
  geom_point(size=3, position="jitter")+scale_color_jco(guide=F)+
  labs(y="delta ( Invertebrate Biomass )")+
  geom_hline(yintercept=0, lty=2)  
plot_grid(TinvBMplot, TinvCount, InvSlopplot, ncol=3,labels = "AUTO")

fishinv<-lm(Count~Treatment+Week, data=InvBMSum[InvBMSum$Week==1 |
                                           InvBMSum$Week==2, ])
summary(fishinv)
anova(fishinv)

library(lsmeans)
leastm<-lsmeans(fishinv, "Week", adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

BiomGraph<-InvBMSum %>% 
  filter(Week==1 | Week==2) %>% 
  gather(variable, value, -c(TxW,Tank,Week,Treatment,NewTreat)) %>% 
  filter(variable!="Richness" & variable!="BMDensity")

ggplot(BiomGraph[BiomGraph$variable="TotalBiomass",], 
       aes(x=Week, y=value, fill=Treatment))+
  geom_boxplot()+facet_wrap(~variable, scale="free")+fishTheme

ggplot(BiomGraph, aes(x=Week, y=value, color=NewTreat))+
  geom_point(position=position_dodge(width=.4))+
  facet_wrap(~variable, scale="free")+theme_bw()
  fishTheme

fishinv<-lm(TotalBiomass~Treatment, data=InvBMSum)
summary(fishinv)
anova(fishinv)

library(lsmeans)
leastm<-lsmeans(fishinv, "Treatment",cadjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

###

Surv.Inv<-BiomGraph %>% filter(variable=="TotalBiomass" & Week==1) %>% 
  spread(variable, value) %>% left_join(TankMean) %>% 
  select(-meanWeightChange, -meanInfect)
summary(aov(meanDaysSurv~Treatment*TotalBiomass, data=Surv.Inv))

model <- lm(meanDaysSurv ~ Treatment * TotalBiomass, data=Surv.Inv)
drop1(model, test="F")


###### Invertebrates in Water Column Samples ######
#area sampled in cm3
library(tidyverse)
library(readxl)
WCinvertRaw<-read_excel("./data/CostMutData.xlsx",sheet = "WaterColumnInv")

WCinv<-WCinvertRaw %>% mutate(Sample=paste(Tank, Week, CountN, sep="."),
                              VolSampled=(14*19)*Depth.mm, #volume of petri dish = volume sampled nsamples * area of square * depth
                              VolTotal=(8.3/2)^2*pi*Depth.mm,#volume of petri dish = area of dish * depth
                              VolumePull=540*(4*4*pi), #volume of the tank sampled
                              DensityNL=(Count/VolSampled * VolTotal/VolumePull)/1e-6) %>% 
  group_by(Sample, Tank, Week) %>%
  summarise(InvDen=sum(DensityNL)) %>% group_by(Tank, Week) %>% 
  summarise(AvgDensityL=mean(InvDen))

WCcom<-WCinvertRaw %>% mutate(Sample=paste(Tank, Week, CountN, sep=".")) %>% 
  select(-Depth.mm, -Size.mm, -CountN,-Week,-Tank) %>% 
  spread(Taxa, Count) %>% replace_na(list(ChironomidaeL=0,
                                          Copepoda=0,
                                          Daphnia=0,
                                          DipteraAdult=0,
                                          Dytiscidae=0,
                                          Oligocheata=0))
##### citations
View(semi_join(BiomassReg,InvBioMass[,8:10], by="Family"))
semi_join(BiomassReg,InvBioMass[,8:10], by="Family")[,c("Primary Source","Secondary Source")] %>%
  filter(is.na(`Secondary Source`)) %>% group_by(`Primary Source`) %>% tally()
#plus Benke et al 1999