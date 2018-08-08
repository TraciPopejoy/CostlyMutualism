library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(xlsx)
library(vegan)

#############     Invert Data in Benthic Samples    ##############
#bring in enclosure and treatment data, trait/taxa list, and length weight regressions
Treat<-read.csv("./FEn17_data/FEn17OKTreatments.csv", sep=",", stringsAsFactors = F) 
BiomassReg<-read.xlsx("./FEn17_data/Macroinv Power Law Coeffs TBP.xlsx", sheetIndex = 1, stringsAsFactors=F)
#THIS IS THE INSECT DATA

#clean the data frame
Inv$Treatment<-Treat[match(Inv$Enc, Treat$Enclosure2), "TreatA"]
sort(unique(Inv$Taxa)) #check to make sure no misspellings
Inv$Treatment<-factor(Inv$Treatment, c("CTRL","ACTL","ACTS","AMBL","AMBS"))
###table is wrong
treattype<-data.frame(Treatment=na.omit(unique(Inv$Treatment)),
                      Type=factor(c("Live","Sham","Ctrl","Sham","Live"),levels=c("Ctrl","Sham","Live")), 
                      Spp=factor(c("ACT","ACT","Ctrl","AMB","AMB"), levels=c("Ctrl","AMB","ACT")))
Inv$Type<-treattype[match(Inv$Treatment, treattype$Treatment),"Type"]
Inv$Spp<-treattype[match(Inv$Treatment,treattype$Treatment),"Spp"]
Inv$Family<-as.character(TaxaList$Family[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Order<-as.character(TaxaList$Order[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Length<-Inv$Length.cm*10

#how many individuals of each taxa in each enclosure/time; long format
Counts<-ddply(Inv, .variables = c("TEid","Taxa"), .fun=function(x) {
  data.frame(TEid=x[1,5],
             Enc=x[1,6],
             Week=x[1,7],
             Treatment=x[1,8],
             n=count(x,x[1,4])[,2])
})
Counts$Family<-as.character(TaxaList$Family[match(Counts$Taxa,TaxaList$Taxa)])
Counts$Order<-as.character(TaxaList$Order[match(Counts$Taxa,TaxaList$Taxa)])

#converting to density to compensate for different sampling effort
head(SlurryData) #found in Slurry Analysis sheet
Counts$Density.npm<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5]*.03315)
Counts$Density.npb<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5])

#############     Field Invert Biomass Calculation     #############
#removing taxa that give me trouble 
InvA<-Inv[!Inv$Order=="misc",]
InvB<-na.omit(InvA)
#apply appropriate biomass regressions to each length
InvBM<-ddply(.data=InvB, .var=c("Taxa"), .fun=function(x) {
  idx<-x[1,c("Family","Order")] #what family/order are we on
  if(idx$Family=="misc"){ #if not ID'd to family, use order level regressions
    plcoeffs<-BiomassReg[BiomassReg$Order == idx$Order &
                           !is.na(BiomassReg$Order == idx$Order),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==idx$Family&
                           BiomassReg$Order == idx$Order&
                           !is.na(BiomassReg$Family==idx$Family), ]  
  }
  #which regressions were actually built for insects this size
  ldply(lapply(1:dim(x)[[1]],FUN=function(i) { 
    idx2<- c(x$Length[i]>=plcoeffs[,19] & x$Length[i]<=plcoeffs[,20]) 
    idx2[is.na(idx2)]<-T #if no size range listed, use the regression anyways
    d1<-plcoeffs[idx2,]
    indmassest<-d1$a*(x$Length[i]^d1$b) #power law to determine biomass
    data.frame(
      TEid=x$TEid[i],
      Enc=x$Enc[i],
      Week=x$Week[i],
      Treatment=x$Treatment[i],
      length=x$Length[i],
      neq=length(idx2), #number of possible equations used
      ninR=sum(idx2), #number of equations used
      meanBM.mg=mean(indmassest),
      median=median(indmassest),
      stDev=sd(indmassest))}), 
    data.frame)
})
#get mean biomass and sum of each taxa for each enclosure/time
InvTotalBM<-ddply(InvBM, .var=c("TEid", "Taxa"),
                  .fun=function(x){data.frame(mean.mg=mean(x$meanBM.mg, na.rm=T), 
                                              median.mg=median(x$meanBM.mg, na.rm=T), 
                                              SD.mg=sd(x$meanBM.mg, na.rm=T),
                                              Sum.mg=sum(x$meanBM.mg,na.rm=T),
                                              Treatment=x$Treatment[1],
                                              Enc=x$Enc[1],
                                              Week=x$Week[1],
                                              mean.length=mean(x$length, na.rm=T)
                  )}) 
#converting mean biomass to biomass/meter
InvTotalBM$BMDensity<-InvTotalBM$Sum.mg/(SlurryData[match(InvTotalBM$TEid, SlurryData$TEid),"Basket."]*.03315)
###would use density because sampling was not constant (not always full basket recovery)

#get each taxa and teid with counts and biomass, and density of both
InvGraph<-merge(InvTotalBM[,-c(4,5)],Counts, by=c("TEid","Taxa","Treatment","Enc","Week"))
#get traits into the graph
InvGraph<-merge(InvGraph, TaxaList[,c(1,2,5:9)], by=c("Taxa", "Family"))

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


#### graph city ####

fungraph<-theme(axis.text.x=element_text(angle = 00,size=12,color="black", hjust=1),
                axis.text.y = element_text(size=12,color="black"),
                axis.title.y=element_text(size=20),
                plot.background = element_blank(),
                panel.border=element_blank(),
                panel.grid.major= element_line(colour=NA), 
                panel.grid.minor=element_line(colour=NA),
                title=element_text(size=20),
                panel.background = element_rect(fill = "white"),
                legend.key=element_rect(colour=NA), 
                axis.line.x=element_line(colour="black"),
                axis.line.y=element_line(colour="black"),
                strip.background=element_rect(fill="white", color="black"),
                strip.text=element_text(size=15))
library(colorspace)
CP<-diverge_hcl(5, h=c(180,70), c = 100, l = c(50, 90), power = 1)
CP[3]<-"black"
CP2<-data.frame(colorss=CP[c(3,1,2,5,4)], Treat=unique(Traitplot$Treatment))