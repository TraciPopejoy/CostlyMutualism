library(reshape2); library(ggplot2); library(tidyverse); 
library(vegan); library(readxl)

#############     Invert Data in Benthic Samples    ##############
#bring in enclosure and treatment data, trait/taxa list, and length weight regressions
treat<-read_excel("./data/CostMutData.xlsx",sheet="TankData")
head(treat)
BiomassReg<-read_excel("./data/Macroinv Power Law Coeffs TBP.xlsx") #regressions excel sheet
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

InvClean<-Inv %>% left_join(invkey) %>% filter(Macro.Meio=="macro")
InvClean %>% left_join(invkey, by="Taxa") %>% filter(is.na(Family.x))

#############    Invert Biomass Calculation     #############
#apply appropriate biomass regressions to each length
biomass<-function(Fam, Ord, Length) {
  #tmass<-NULL
  #if not ID'd to family, use order level regressions
  if(is.na(match(Fam, BiomassReg$Family))){ 
    plcoeffs<-BiomassReg[BiomassReg$Order == Ord &
                           !is.na(BiomassReg$Order == Ord),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==Fam&
                           BiomassReg$Order == Ord&
                           !is.na(BiomassReg$Family==Fam), ] }
  #which regressions were actually built for insects this size
  idx2<- c(Length>=plcoeffs[,19] & Length<=plcoeffs[,20]) 
  idx2[is.na(idx2)]<-T #if no size range listed, use the regression anyways
  d1<-plcoeffs[idx2,] #row by row????
  mass<-d1$a*(Length^d1$b) #power law to determine biomass
  #tmass<-cbind(tmass,mass)#get vector of mass values from different eq.
  mean(mass, na.rm=T)
  #print(paste(Fam,Ord)) #physidae order is messed up
  #add error message and a nA option
} #try #trycatch #stackoverflow
#for loop with a print

InvBioMass<-InvClean %>% select(-Notes,-ID.person) %>%
  rowwise() %>%
  mutate(BM=biomass(Family,Order, Length.mm), 
         TotalBM=Count*BM)
library(lubridate)
InvBMSum<-InvBioMass %>% 
  group_by(TxW, Tank, Week, NewTreat, Treatment) %>% 
    summarize(TotalBiomass=sum(TotalBM, na.rm=T)/1000,
              Richness=length(unique(Taxa)),
              Count=sum(Count, na.omit=T),
              BMDensity=TotalBiomass/(3*.0103),
              CountDensity=Count/(3*.0103),
              AvgLength=mean(Length.mm, na.rm=T)) %>%
  mutate(Avg.BM=TotalBiomass/AvgLength,
         Date=case_when(Week==0~ymd("2018-06-12"),
                        Week==1~ymd("2018-06-18"),
                        Week==2~ymd("2018-06-29"),
                        Week==3~ymd("2018-07-06"),
                        Week==4~ymd("2018-07-13")))
InvBMSum %>% group_by(Week, NewTreat) %>% tally() %>% spread(Week,n)
#write_excel_csv(InvBMSum, "DeathInvertebrate1025.xls")

ggplot(InvBMSum[InvBMSum$Week!=0,], 
       aes(x=Date, y=BMDensity, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1.5))+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), 
               position=position_dodge(width=1.5))+
  ylab(expression("Invertebrate Biomass mg "%*%m^-2))+
  geom_vline(xintercept=ymd("2018-07-02"))+theme_bw()
ggplot(InvBMSum[InvBMSum$Week!=0,], 
       aes(x=Date, y=CountDensity, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1.5))+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), 
               position=position_dodge(width=1.5))+
  ylab(expression("Invertebrates "%*%m^-2))+
  geom_vline(xintercept=ymd("2018-07-02"))+theme_bw()
ggplot(InvBMSum[InvBMSum$Week!=0,], 
       aes(x=Date, y=Richness, fill=NewTreat, color=NewTreat))+
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1.5))+
  stat_summary(aes(fill=NewTreat, shape=NewTreat), 
               position=position_dodge(width=1.5))+
  ylab("Taxonomic Richness")+
  geom_vline(xintercept=ymd("2018-07-02"))+theme_bw()

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
  scale_fill_manual(values=pal_futurama()(6), name="Treat.Time", guide=F)+nmdstheme+
  theme(legend.position = c(.88,.7), legend.text = element_text(size=rel(.55)),
        legend.title= element_text(size=rel(.7)),
        legend.background = element_rect(fill=alpha(0.1)))
library(cowplot);library(ggsci)
plot_grid(ctrlnmds,livenmds,deadnmds,spnmds, labels=c("Control","Live","Dead", ""))
ggsave("./Figures/NMDSwhole.png")

#exploring differences in life stages
IBMlifestage<-InvBioMass %>% filter(Macro.Meio=="macro" & !is.na(LifeStage)) %>%
  group_by(Tank, Week, TxW, NewTreat, LifeStage) %>%
  summarize(LSTotal=sum(Count))
ggplot(IBMlifestage,aes(x=Tank,y=LSTotal, fill=LifeStage))+
  geom_bar(stat="identity")+facet_wrap(~Week+NewTreat, scales="free_x")+
  theme_bw()
ggplot(IBMlifestage[IBMlifestage$LifeStage=="pupa",], 
       aes(x=Week, y=LSTotal, color=NewTreat))+
  geom_boxplot(aes(group=interaction(Week, NewTreat)))+
  geom_point(position=position_dodge(width=.75))+
  #stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=1.5))+
  #stat_summary(aes(fill=NewTreat, shape=NewTreat), 
  #             position=position_dodge(width=1.5))+
  ylab("Invertebrate Count")+
  geom_vline(xintercept=ymd("2018-07-02"))+theme_bw()

#exploring Order composition
IBMorder<-InvBioMass %>% filter(Macro.Meio=="macro" & Week !=0) %>%
  group_by(Tank, Week, TxW, NewTreat, Order) %>%
  summarize(OTotal=sum(Count))
ggplot(IBMorder[IBMorder$NewTreat=="Dead",],aes(x=Tank,y=OTotal, fill=Order))+
  geom_bar(stat="identity")+facet_wrap(~Week)+
  theme_bw()
ggplot(IBMorder[IBMorder$NewTreat=="Live",],aes(x=Tank,y=OTotal, fill=Order))+
  geom_bar(stat="identity")+facet_wrap(~Week)+
  theme_bw()
ggplot(IBMorder[IBMorder$NewTreat=="Control",],aes(x=Tank,y=OTotal, fill=Order))+
  geom_bar(stat="identity")+facet_wrap(~Week)+
  theme_bw()

###### Invertebrates in Water Column Samples ######
#area sampled in cm3
library(tidyverse);library(readxl)
WCinvertRaw<-read_excel("./data/CostMutData.xlsx",sheet = "WaterColumnInv",
                        col_types = c("text","numeric","text","text","numeric","numeric","numeric"))
WCinvInd<-WCinvertRaw %>% 
  mutate(Sample=paste(Tank, Week, CountN, sep="."),
        VolSampled=(14*19)*Depth.mm, #volume of petri dish = volume sampled nsamples * area of square * depth
        VolTotal=(8.3/2)^2*pi*Depth.mm,#volume of petri dish = area of dish * depth
        VolumePull=540*(4*4*pi), #volume of the tank sampled
        DensityNL=(Count/VolSampled * VolTotal/VolumePull)/1e-6) %>% 
  left_join(invkey)%>% filter(!is.na(Family) & Order!="misc" & 
                                Order!="Copepoda" & Family!="Hydracarina") %>%
  rowwise()%>%
  mutate(BMest.mg=mean(biomass(Family, Order, Size.mm))*Count,
         BMestDens=(BMest.mg/VolSampled * VolTotal/VolumePull)/1e-6) %>% 
  filter(!is.na(Macro.Meio))
WCinv<- WCinvInd %>% group_by(Tank, Week, Macro.Meio, Sample) %>%
  summarise(InvDen.n.perL=sum(DensityNL),
                              sumBMtrial=sum(BMestDens)) %>% 
  summarise(AvgDensityL.n.perL=mean(InvDen.n.perL, na.rm=T),
            meanBM.mg.L=mean(sumBMtrial, na.rm=T)) %>%
  left_join(treat) %>% 
  select(Tank, Week, AvgDensityL.n.perL, meanBM.mg.L, Treatment, Macro.Meio) %>%
  gather(measurement, value, -Tank, -Week, -Treatment, -Macro.Meio) 

WCinv2<-WCinv %>%
  mutate(part.treat=paste(Treatment,Macro.Meio, sep="."),
         treat.week=paste(Treatment,Week, sep="."),
         fme=factor(measurement, 
                    levels=c("AvgDensityL.n.perL","meanBM.mg.L"),
                    labels=c('"Invertebrates L"^-1', 
                             '"Invertebrate Biomass mg L"^-1'))) %>% 
  filter(Week==1|Week==2)
                         
ggplot(WCinv2,
       aes(x=part.treat, y=value, fill=Treatment))+
  geom_boxplot()+facet_wrap(~Week+measurement, scales="free_y")
mypal<-c("#0073C2","#42B4FF","#EFC000","#FFE16B")
show_col(mypal)
wcdots<-ggplot(WCinv2, aes(x=Week,y=value, color=part.treat, group=part.treat))+
  stat_summary(fun.y = mean, geom = "line", position=position_dodge(width=.2), size=1.1)+
  stat_summary(position=position_dodge(width=.2), size=1.2)+
  scale_y_continuous(name="")+
  scale_x_continuous(name="Sampling Day",breaks=c(1,2),labels=c("Day -1", "Day 11"))+
  scale_color_manual(values=mypal, name="Treatment-\nInvertebrate Size", 
                    labels=c("Control-macrofauna   ","Control-meiofauna   ",
                             "Mussel-macrofauna","Mussel-macrofauna"))+
  facet_wrap(~fme, scales="free_y", strip.position="left", 
             labeller = label_parsed)+
  theme(strip.background =element_rect(fill=NA),
        strip.placement="outside")
ggdraw(wcdots) + 
  draw_label("A", x=.1, y=.95) + 
  draw_label("B", x=.42, y=.95)
ggsave("FishWCinvDots.tiff",dpi=300, width=7, height=3)

wcnplot<-ggplot()+
  geom_boxplot(data=WCinv2[WCinv2$Macro.Meio=="meio" & WCinv2$measurement=="AvgDensityL.n.perL",],
               aes(x=Week,y=value, fill=part.treat, group=interaction(Week,part.treat)),
               position=position_dodge(width=.5))+
  geom_boxplot(data=WCinv2[WCinv2$Macro.Meio=="macro"& WCinv2$measurement=="AvgDensityL.n.perL",],
               aes(x=Week,y=value, fill=part.treat, group=interaction(Week,part.treat)),
               position=position_dodge(width=.5))+
  scale_fill_manual(values=mypal, guide=F)+
  scale_y_continuous(name=expression("Invertebrates L"^-1))+
  scale_x_continuous(name="",breaks=c(1,2),labels=c("Day -1", "Day 11"))+
  facet_wrap(~Treatment, strip.position = "bottom")+
  theme(strip.background =element_rect(fill=NA),
        strip.placement = "outside")
wcbmplot<-ggplot()+
  geom_boxplot(data=WCinv2[WCinv2$Macro.Meio=="meio" & WCinv2$measurement=="meanBM.mg.L",],
               aes(x=Week,y=value, fill=part.treat, group=interaction(Week,part.treat)),
               position=position_dodge(width=.5))+
  geom_boxplot(data=WCinv2[WCinv2$Macro.Meio=="macro"& WCinv2$measurement=="meanBM.mg.L",],
               aes(x=Week,y=value, fill=part.treat, group=interaction(Week,part.treat)),
               position=position_dodge(width=.5))+
  scale_fill_manual(values=mypal, name="Treatment-Invertebrate Size", 
                    labels=c("Control-macrofauna   ","Control-meiofauna   ",
                             "Mussel-macrofauna","Mussel-macrofauna"))+
  scale_y_continuous(name=expression("Invertebrate Biomass mg L"^-1))+
  scale_x_continuous(name="",breaks=c(1,2),labels=c("Day -1", "Day 11"))+
  facet_wrap(~Treatment, strip.position = "bottom")+
  guides(fill=guide_legend(ncol=2))+
  theme(strip.background =element_rect(fill=NA),
        strip.placement = "outside",
        legend.position = "none")
library(cowplot)
legend<-get_legend(wcbmplot+theme(legend.position =c(.25,.8)))
prow<-plot_grid(wcnplot,wcbmplot, labels="AUTO")
plot_grid(prow,legend, nrow=2, rel_heights = c(1,.12))
ggsave("FishWCinv.tiff", dpi=300, width=6, height=3.8)

WCinv2
WCmod1<-aov(lm(value~Treatment+Week, data=WCinv2[WCinv2$measurement=="AvgDensityL.n.perL",]))
summary(WCmod1)
WCmod2<-aov(lm(value~Treatment+Week, data=WCinv2[WCinv2$measurement=="meanBM.mg.L",]))
summary(WCmod2)


WCinvMData<- WCinvInd %>% group_by(Tank, Week, Sample) %>%
    summarise(InvDen.n.perL=sum(DensityNL),
              sumBMtrial=sum(BMest.mg)) %>% 
    summarise(AvgDensityL.n.perL=mean(InvDen.n.perL, na.rm=T),
              meanBM.mg=mean(sumBMtrial, na.rm=T)) %>%
    left_join(treat) %>% 
    select(Tank, Week, AvgDensityL.n.perL, meanBM.mg, Treatment)

library(lmerTest)
WCinvMacroM<-lmer(meanBM.mg~Treatment+Week+(1|Tank), data=WCinvMData)
summary(WCinvMacroM)
anova(WCinvMacroM)
    
WCinvMData %>% group_by(Week, Treatment) %>% tally() %>% spread(Treatment, n)

##### citations
View(semi_join(BiomassReg,InvBioMass[,8:10], by="Family"))
semi_join(BiomassReg,InvBioMass[,8:10], by="Family")[,c("Primary Source","Secondary Source")] %>%
  filter(is.na(`Secondary Source`)) %>% group_by(`Primary Source`) %>% tally()
#plus Benke et al 1999