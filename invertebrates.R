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
library(lubridate)
tempinv<-InvBMSum %>% select(-Treatment, -Avg.BM,-Week.n) %>%
  mutate(Date=case_when(Week==0~ymd("2018-06-12"),
                        Week==1~ymd("2018-06-18"),
                        Week==2~ymd("2018-06-29"),
                        Week==3~ymd("2018-07-06"),
                        Week==4~ymd("2018-07-13")))
#write_excel_csv(tempinv, "DeathInvertebrate1025.xls")

ggplot(InvBMSum, aes(x=Week, y=TotalBiomass, fill=NewTreat))+
  geom_boxplot()+facet_wrap(~Week, scale="free")+theme_bw()

invgraph<-InvBMSum %>% left_join(physchem, by=c("Week.n"="Week","Tank")) %>%
  filter(Week.n>.5) %>%
  select(-Temp.C, -Cond.uS, -DO.mgL, -WaterV.mLs, -Chl1, -Chl2)

Invplot<-ggplot(invgraph, aes(x=Date,y=BMDensity, 
                     fill=NewTreat, group=interaction(Date,NewTreat)))+
  geom_boxplot()+geom_vline(xintercept=ymd("2018-07-02"))+
  ylab(expression("Invertebrate Biomass mg "%*%m^-2)) +
  scale_fill_grey(start=0.3, end=0.9, name="Treatment")+fronteirstheme
ggsave("Fig4.tiff",Invplot, width=3.34, height=3.34, dpi=300)

#### NMDS ####
InvComMat<-InvBioMass %>% 
  group_by(TxW, Taxa) %>% summarize(sumC=sum(Count)) %>%
  spread(Taxa, sumC) %>% left_join(InvBMSum[,1:5]) %>%
  replace(., is.na(.),0)

library(vegan)
nmds<-metaMDS(InvComMat[,-c(1,23:26)])
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- InvComMat$TxW  # create a column of site names, from the rownames of data.scores
data.scores$grp <- InvComMat$NewTreat  #  add the grp variable created earlier
data.scores$week<- InvComMat$Week
data.scores$time<- case_when(data.scores$week==0 ~ "before",
                             data.scores$week==1 ~ "before",
                             data.scores$week==2 ~ "before",
                             data.scores$week==3 ~ "after",
                             data.scores$week==4 ~ "after")
data.scores$grptime<-paste(data.scores$grp, data.scores$time, sep=".")
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

grp.a <- data.scores[data.scores$grp == "Control", ][chull(data.scores[data.scores$grp == 
                                                                         "Control", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "Dead", ][chull(data.scores[data.scores$grp == 
                                                                      "Dead", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$grp == "Live", ][chull(data.scores[data.scores$grp == 
                                                                      "Live", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.before <- data.scores[data.scores$time=="before", ][chull(data.scores[data.scores$time=="before",
                                                                          c("NMDS1", "NMDS2")]), ]
grp.after <- data.scores[data.scores$time=="after", ][chull(data.scores[data.scores$time=="after",
                                                                        c("NMDS1", "NMDS2")]), ]


grp.cb <- data.scores[data.scores$grptime=="Control.before", ][chull(data.scores[data.scores$grptime=="Control.before",
                                                                                 c("NMDS1", "NMDS2")]), ]
grp.ca <- data.scores[data.scores$grptime=="Control.after", ][chull(data.scores[data.scores$grptime=="Control.after",
                                                                                c("NMDS1", "NMDS2")]), ]
grp.db <- data.scores[data.scores$grptime=="Dead.before", ][chull(data.scores[data.scores$grptime=="Dead.before",
                                                                              c("NMDS1", "NMDS2")]), ]
grp.da <- data.scores[data.scores$grptime=="Dead.after", ][chull(data.scores[data.scores$grptime=="Dead.after",
                                                                             c("NMDS1", "NMDS2")]), ]
grp.lb <- data.scores[data.scores$grptime=="Live.before", ][chull(data.scores[data.scores$grptime=="Live.before",
                                                                              c("NMDS1", "NMDS2")]), ]
grp.la <- data.scores[data.scores$grptime=="Live.after", ][chull(data.scores[data.scores$grptime=="Live.after",
                                                                             c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c)  #combine grp.a and grp.b
hull.data.time <- rbind(grp.before, grp.after)  #combine grp.a and grp.b
hull.data.both <- rbind(grp.cb, grp.ca, grp.da, grp.db, grp.lb, grp.la)  #combine grp.a and grp.b
nmdstheme<-theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=rel(.7)), # remove x-axis labels
        axis.title.y = element_text(size=rel(.7)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

ctrlnmds<-ggplot() + 
  geom_polygon(data=hull.data.both[hull.data.both$grp=="Control",],
               aes(x=NMDS1,y=NMDS2,fill=grptime,group=grptime),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,
            aes(x=NMDS1,y=NMDS2,label=species),alpha=0, size=3) +  # add the species labels
  geom_text(data=data.scores[data.scores$grp=="Control",],aes(x=NMDS1,y=NMDS2, label=site),size=4) + # add the point markers
  scale_fill_manual(values=pal_futurama()(6)[1:2], guide=F)+
  coord_equal() + nmdstheme
deadnmds<-ggplot() + 
  geom_polygon(data=hull.data.both[hull.data.both$grp=="Dead",],
               aes(x=NMDS1,y=NMDS2,fill=grptime,group=grptime),alpha=0.20) +
  geom_text(data=species.scores,
            aes(x=NMDS1,y=NMDS2,label=species),alpha=0, size=3) +  # add the species labels
  geom_text(data=data.scores[data.scores$grp=="Dead",],aes(x=NMDS1,y=NMDS2,label=site),size=4) + # add the point markers
  scale_fill_manual(values=pal_futurama()(6)[3:4], guide=F)+
  coord_equal() + nmdstheme
livenmds<-ggplot() + 
  geom_polygon(data=hull.data.both[hull.data.both$grp=="Live",],
               aes(x=NMDS1,y=NMDS2,fill=grptime,group=grptime),alpha=0.20) +
  geom_text(data=species.scores,
            aes(x=NMDS1,y=NMDS2,label=species),alpha=0, size=3) +  # add the species labels
  geom_text(data=data.scores[data.scores$grp=="Live",],aes(x=NMDS1,y=NMDS2,label=site),size=4) + # add the point markers
  scale_fill_manual(values=pal_futurama()(6)[5:6], guide=F)+
  coord_equal() + nmdstheme
spnmds<-ggplot() + 
  geom_polygon(data=hull.data.both, aes(x=NMDS1,y=NMDS2,fill=grptime,group=grptime),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,
            aes(x=NMDS1,y=NMDS2,label=species), size=3) +
  scale_fill_manual(values=pal_futurama()(6), name="Treat.Time")+nmdstheme+
  theme(legend.position = c(.88,.7), legend.text = element_text(size=rel(.55)),
        legend.title= element_text(size=rel(.7)),
        legend.background = element_rect(fill=alpha(0.1)))
plot_grid(ctrlnmds,livenmds,deadnmds,spnmds, labels=c("Control","Live","Dead", ""))
ggsave("./Figures/NMDSwhole.png")


#### FISH MANUSCRIPT ####
library(ggsci)
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
InvBMslope<-InvBMSum %>% ungroup() %>%
  filter(Week==1 | Week==2) %>%
  select(Tank,Week, TotalBiomass) %>% spread(Week, TotalBiomass) %>%
  mutate(InvBMrate=(`2`-`1`)/12)
slopessss<-InvBMslope %>% left_join(treat)
InvSlopplot<-ggplot(slopessss, 
                  aes(x=Treatment, y=InvBMrate, color=Treatment))+
  geom_point(size=3, position="jitter")+scale_color_jco(guide=F)+
  labs(y=expression("delta ( Invertebrate Biomass )  "%*%day^ -1))+
  geom_hline(yintercept=0, lty=2)  
library(cowplot)
fishInvPlot<-plot_grid(TinvBMplot, TinvCount, InvSlopplot, ncol=3,labels = "AUTO")
ggsave("Fig1.tiff", fishInvPlot, dpi=300)

library(lmerTest)
fishinv<-lmer(Count~Treatment+(1|Week)+(1|Tank), data=InvBMSum[InvBMSum$Week==1 |
                                           InvBMSum$Week==2, ])
summary(fishinv)
anova(fishinv)

fishbio<-lmer(BMDensity~Treatment+(1|Week)+(1|Tank), data=InvBMSum[InvBMSum$Week==1 |
                                                                 InvBMSum$Week==2, ])
summary(fishbio)
anova(fishbio)

#### NMDS ####
InvComMatfish<-InvBioMass %>% filter(Week=="1" | Week=="2") %>%
  group_by(TxW, Taxa) %>% summarize(sumC=sum(Count)) %>%
  spread(Taxa, sumC) %>% left_join(InvBMSum[,1:5]) %>%
  replace(., is.na(.),0)

library(vegan)
nmds<-metaMDS(InvComMatfish[,-c(1,21:24)])
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- InvComMatfish$TxW  # create a column of site names, from the rownames of data.scores
data.scores$grp <- InvComMatfish$Treatment  #  add the grp variable created earlier
data.scores$week<- InvComMatfish$Week
data.scores$time<- case_when(data.scores$week==1 ~ "before",
                             data.scores$week==2 ~ "after")
data.scores$grptime<-paste(data.scores$grp, data.scores$time, sep=".")
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

grp.a <- data.scores[data.scores$grp == "Control", ][chull(data.scores[data.scores$grp == 
                                                                         "Control", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "Mussel", ][chull(data.scores[data.scores$grp == 
                                                                      "Mussel", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
nmdsFISH<-ggplot() + 
  geom_polygon(data=hull.data, aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,
            aes(x=NMDS1,y=NMDS2,label=species),alpha=0, size=3) +  # add the species labels
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2, label=site),size=4) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2, label=species),
            size=3, alpha=.7, color="darkgray") + # add the point markers
  scale_fill_manual(values=pal_jco()(2))+
  coord_equal() + nmdstheme
nmdsFISH

comdistance<-data.scores %>% left_join(InvComMatfish, by=c("site"="TxW"))%>% 
  select(NMDS1,NMDS2, week, Tank) %>% mutate(week2=paste(week,"n")) %>%
  spread(week,NMDS1) %>% spread(week2, NMDS2) %>% group_by(Tank) %>%
  summarise_all(funs(.[which(!is.na(.))])) %>% left_join(treat) %>%
  select(Tank, Treatment, `1`,`2`,`1 n`,`2 n`) %>% 
  mutate(InvComDistance=sqrt(abs((`2`-`1`)^2+(`2 n`-`1 n`)^2)),
         InvComSlope=InvComDistance/12)
ggplot(comdistance, aes(x=Treatment, y=InvComDistance)) + geom_boxplot()
summary(lm(InvComDistance~Treatment, data=comdistance))

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