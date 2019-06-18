
library(brms)
#bring in data from other scripts; already cleaned
source("Producers.R")
head(ChlSummary)
ChlSumMODEL<-ChlSummary %>% 
  mutate(TreF=factor(NewTreat, levels=c("Control","Live","Dead")),
         TankF=factor(Tank),
         DayFR=factor(Day, levels=c(-20,-15,-3,4,11,25,39)))
ChlSumMODEL[ChlSumMODEL$BenthicChlA.ug.cm<=0,"BenthicChlA.ug.cm"]<-0.0001
ChlSumMODEL[ChlSumMODEL$WaterColChlA.ug.L<=0,"WaterColChlA.ug.L"]<-1e-3

source("waterchem.R")
head(WaterNutrients)
WaterNutMODEL<-WaterNutrients %>%
  mutate(TreF=factor(NewTreat, levels=c("Control","Live","Dead")),
       DayF=factor(Day, levels=c(-20,-15,-3,4,18,25,39)))

set.seed(6363) #set random number seed for replicable analysis

### testing priors 
hist(rstudent_t(1000, 3,0,10)) #default prior on sd(TankIntercept)
hist(rgamma(1000,0.01,0.01), xlim=c(0,70)) # prior on shape of Treat:Day relationship
hist(rnorm(1000, 0,5)) # prior on beta (Treat:Day interaction)

#### Ammonium analysis ####
summary(WaterNutMODEL$FilterdNH3ugL)
nhBmodel<-brm(FilterdNH3ugL~TreF:DayF+(1|TankF)-1, data=WaterNutMODEL,
     family=Gamma(link='log'), 
     prior=c(prior(normal(0,5), class="b"),
             prior(gamma(0.01,0.01), class="shape")),
     chains=4, iter=2000)
print(nhBmodel, prior=T)
plot(marginal_effects(nhBmodel)) #check it is reproducing data well
#pull out posterior samples for each parameter
nh_post_model<-posterior_samples(nhBmodel) 
nhRRraw<-as_tibble(nh_post_model) %>%
  dplyr::select(-starts_with("r_Tank")) %>% #don't care about tank intercepts
  mutate( #response ratios between dead and control treatments
         rrDCd.20=`b_TreFDead:DayFM20`/`b_TreFControl:DayFM20`,
         rrDCd.15=`b_TreFDead:DayFM15`/`b_TreFControl:DayFM15`,
         rrDCd.3=`b_TreFDead:DayFM3`/`b_TreFControl:DayFM3`,
         rrDCd4=`b_TreFDead:DayF4`/`b_TreFControl:DayF4`,
         rrDCd18=`b_TreFDead:DayF18`/`b_TreFControl:DayF18`,
         rrDCd25=`b_TreFDead:DayF25`/`b_TreFControl:DayF25`,
         rrDCd39=`b_TreFDead:DayF39`/`b_TreFControl:DayF39`,
         #response ratios between live and control treatments
         rrLCd.20=`b_TreFLive:DayFM20`/`b_TreFControl:DayFM20`,
         rrLCd.15=`b_TreFLive:DayFM15`/`b_TreFControl:DayFM15`,
         rrLCd.3=`b_TreFLive:DayFM3`/`b_TreFControl:DayFM3`,
         rrLCd4=`b_TreFLive:DayF4`/`b_TreFControl:DayF4`,
         rrLCd18=`b_TreFLive:DayF18`/`b_TreFControl:DayF18`,
         rrLCd25=`b_TreFLive:DayF25`/`b_TreFControl:DayF25`,
         rrLCd39=`b_TreFLive:DayF39`/`b_TreFControl:DayF39`,
         #response ratios between dead and live treatments
         rrDLd.20=`b_TreFDead:DayFM20`/`b_TreFLive:DayFM20`,
         rrDLd.15=`b_TreFDead:DayFM15`/`b_TreFLive:DayFM15`,
         rrDLd.3=`b_TreFDead:DayFM3`/`b_TreFLive:DayFM3`,
         rrDLd4=`b_TreFDead:DayF4`/`b_TreFLive:DayF4`,
         rrDLd18=`b_TreFDead:DayF18`/`b_TreFLive:DayF18`,
         rrDLd25=`b_TreFDead:DayF25`/`b_TreFLive:DayF25`,
         rrDLd39=`b_TreFDead:DayF39`/`b_TreFLive:DayF39`,
         #total BACI response ratios for each comparison
         BACIdc=10^((rrDCd39+rrDCd25+rrDCd18+rrDCd4)/4)/10^((rrDCd.20+rrDCd.15+rrDCd.3)/3),
         BACIdl=10^((rrDLd39+rrDLd25+rrDLd18+rrDLd4)/4)/10^((rrDLd.20+rrDLd.15+rrDLd.3)/3),
         BACIlc=10^((rrLCd39+rrLCd25+rrLCd18+rrLCd4)/4)/10^((rrLCd.20+rrLCd.15+rrLCd.3)/3))
# graphing weekly response ratios, check with plotted values
nhRRweekgraph <-nhRRraw%>% 
  dplyr::select(starts_with("rr")) %>%
  gather(ratio,value) %>% 
  mutate(ratio.Type=case_when(substr(ratio, 1,4)=="rrDC"~"Dead / Ctrl",
                              substr(ratio, 1,4)=="rrLC"~"Live / Ctrl",
                              substr(ratio, 1,4)=="rrDL"~"Dead / Live")) %>%
  mutate(ratio.F=factor(ratio, levels=c("rrDCd.20", "rrDCd.15", "rrDCd.3","rrDCd4",
                                        "rrDCd18","rrDCd25", "rrDCd39",
                                        "rrLCd.20", "rrLCd.15", "rrLCd.3","rrLCd4",
                                        "rrLCd18","rrLCd25", "rrLCd39",
                                        "rrDLd.20", "rrDLd.15", "rrDLd.3","rrDLd4",
                                        "rrDLd18","rrDLd25", "rrDLd39")))
ggplot(nhRRweekgraph, 
       aes(x=ratio.F, y=value))+
  geom_hline(yintercept=1, color="red", size=2)+
  geom_boxplot()+
  facet_wrap(~ratio.Type, scales="free_x")+
  scale_x_discrete(label=rep(c(-20,-15,-3,4,18,25,39),3), name="Sampling Day")+
  scale_y_continuous(name="Ammonium log(Response Ratio)", limits = c(0.4,3.1),
                     breaks=c(0.5,1,1.5,2,2.5,3))
nhRRraw %>% dplyr::select(rrDLd4, rrDCd4) %>%
  dplyr::summarize(meanRR4dl=mean(rrDLd4), meanRR4dc=mean(rrDCd4))
# graphing BACI ratios
# these are % probability that response increased (if >1) after the impact
nhBACIgraph<-nhRRraw %>% dplyr::select(starts_with("BACI")) %>%
  gather(ratio,value)
ggplot(nhBACIgraph, aes(x=ratio, y=value))+geom_boxplot()+
  scale_x_discrete(labels=c("BACI dead.ctrl","BACI dead.live","BACI live.ctrl"))

# quantifying the % probability response increased after impact
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1 ) 
nhBACIgraph %>% group_by(ratio) %>% dplyr::summarize(meanBACI=mean(value))
1-ecdf(as.matrix(nhBACIgraph[nhBACIgraph$ratio=="BACIdc",2]))( 1.5 ) 
marginal_effects(nhBmodel)$`TreF:DayF` %>% 
  select(TreF, DayF, estimate__, lower__,upper__) %>%
  filter(DayF==4)

#### SRP Analysis ####
summary(WaterNutMODEL$FilterdSRPugL)
srpBmodel<-brm(FilterdSRPugL~TreF:DayF+(1|TankF)-1, data=WaterNutMODEL,
              family=Gamma(link='log'), 
              prior=c(prior(normal(0,5), class="b"),
                      prior(gamma(0.01,0.01), class="shape")),
              chains=4, iter=2000, inits="0")
print(srpBmodel, prior=T)
plot(marginal_effects(srpBmodel)) #check to see if it is replicating data well
#very large bands around the data, sd prior?
srp_post_model<-posterior_samples(srpBmodel)
srpRRraw<-as_tibble(srp_post_model) %>%
  dplyr::select(-starts_with("r_Tank")) %>%
  mutate( #weekly response ratio between dead and control treatments
         rrDCd.20=`b_TreFDead:DayFM20`/`b_TreFControl:DayFM20`,
         rrDCd.15=`b_TreFDead:DayFM15`/`b_TreFControl:DayFM15`,
         rrDCd.3=`b_TreFDead:DayFM3`/`b_TreFControl:DayFM3`,
         rrDCd4=`b_TreFDead:DayF4`/`b_TreFControl:DayF4`,
         rrDCd18=`b_TreFDead:DayF18`/`b_TreFControl:DayF18`,
         rrDCd25=`b_TreFDead:DayF25`/`b_TreFControl:DayF25`,
         rrDCd39=`b_TreFDead:DayF39`/`b_TreFControl:DayF39`,
         # weekly response ratio between live and control treatments
         rrLCd.20=`b_TreFLive:DayFM20`/`b_TreFControl:DayFM20`,
         rrLCd.15=`b_TreFLive:DayFM15`/`b_TreFControl:DayFM15`,
         rrLCd.3=`b_TreFLive:DayFM3`/`b_TreFControl:DayFM3`,
         rrLCd4=`b_TreFLive:DayF4`/`b_TreFControl:DayF4`,
         rrLCd18=`b_TreFLive:DayF18`/`b_TreFControl:DayF18`,
         rrLCd25=`b_TreFLive:DayF25`/`b_TreFControl:DayF25`,
         rrLCd39=`b_TreFLive:DayF39`/`b_TreFControl:DayF39`,
         # weekly response ratio between dead and live treatments
         rrDLd.20=`b_TreFDead:DayFM20`/`b_TreFLive:DayFM20`,
         rrDLd.15=`b_TreFDead:DayFM15`/`b_TreFLive:DayFM15`,
         rrDLd.3=`b_TreFDead:DayFM3`/`b_TreFLive:DayFM3`,
         rrDLd4=`b_TreFDead:DayF4`/`b_TreFLive:DayF4`,
         rrDLd18=`b_TreFDead:DayF18`/`b_TreFLive:DayF18`,
         rrDLd25=`b_TreFDead:DayF25`/`b_TreFLive:DayF25`,
         rrDLd39=`b_TreFDead:DayF39`/`b_TreFLive:DayF39`,
         # total BACI response ratios for each comparison
         BACIdc=10^((rrDCd39+rrDCd25+rrDCd18+rrDCd4)/4)/10^((rrDCd.20+rrDCd.15+rrDCd.3)/3),
         BACIdl=10^((rrDLd39+rrDLd25+rrDLd18+rrDLd4)/4)/10^((rrDLd.20+rrDLd.15+rrDLd.3)/3),
         BACIlc=10^((rrLCd39+rrLCd25+rrLCd18+rrLCd4)/4)/10^((rrLCd.20+rrLCd.15+rrLCd.3)/3)) 
# graphing weekly response ratios to investigate weekly response
srpRRweekgraph <-srpRRraw %>% 
  dplyr::select(starts_with("rr")) %>%
  gather(ratio,value) %>% 
  #easy graphing categories for fact_wrap
  mutate(ratio.Type=case_when(substr(ratio, 1,4)=="rrDC"~"Dead / Ctrl",
                              substr(ratio, 1,4)=="rrLC"~"Live / Ctrl",
                              substr(ratio, 1,4)=="rrDL"~"Dead / Live")) %>%
  mutate(ratio.F=factor(ratio, levels=c("rrDCd.20", "rrDCd.15", "rrDCd.3","rrDCd4",
                                        "rrDCd18","rrDCd25", "rrDCd39",
                                        "rrLCd.20", "rrLCd.15", "rrLCd.3","rrLCd4",
                                        "rrLCd18","rrLCd25", "rrLCd39",
                                        "rrDLd.20", "rrDLd.15", "rrDLd.3","rrDLd4",
                                        "rrDLd18","rrDLd25", "rrDLd39")))
# graphing results of BACI ratios
# evaluates if response increase (when rr >1) after impact
ggplot(srpRRweekgraph, 
       aes(x=ratio.F, y=value))+
  geom_hline(yintercept=1, color="red", size=2)+
  geom_boxplot()+
  facet_wrap(~ratio.Type, scales="free_x")+
  scale_x_discrete(label=rep(c(-20,-15,-3,4,18,25,39),3), name="Sampling Day")+
  scale_y_continuous(name="SRP Response Ratio", #limits = c(0.4,3.1),
                     breaks=c(0.5,1,1.5,2,2.5,3))
srpRRraw %>%  
  dplyr::summarize(meanrrDCd4=mean(rrDCd4), meanrrDLd4=mean(rrDLd4),
                   meanrrDC39=mean(rrDCd39), meanrrDLd39=mean(rrDLd39))
srpBACIgraph<-srpRRraw %>% dplyr::select(starts_with("BACI")) %>%
  gather(ratio,value)
ggplot(srpBACIgraph, aes(x=ratio, y=value))+geom_boxplot()+
  scale_x_discrete(labels=c("BACI dead.ctrl","BACI dead.live","BACI live.ctrl"))

# quantifying BACI ratio - % prob of an increase, 50% increase 
1-ecdf(as.matrix(srpBACIgraph[srpBACIgraph$ratio=="BACIdc",2]))( 1 ) 
srpBACIgraph %>% group_by(ratio) %>% dplyr::summarize(mean(value))
1-ecdf(as.matrix(srpBACIgraph[srpBACIgraph$ratio=="BACIdl",2]))( 1 )
marginal_effects(srpBmodel)$`TreF:DayF` %>% 
  select(TreF, DayF, estimate__, lower__,upper__) %>%
  filter(DayF==39)

#### Water Column Chlorophyll a Analysis ####
summary(ChlSumMODEL$WaterColChlA.ug.L)
wcchlaBmodel<-brm((WaterColChlA.ug.L+10)~TreF:DayF+(1|TankF)-1, data=ChlSumMODEL,
                  family=Gamma(link='log'), 
                  prior=c(prior(normal(0,6), class="b"),
                          prior(gamma(0.01,0.01), class="shape")),
               chains=4, iter=2000)
print(wcchlaBmodel, prior=T) 
plot(marginal_effects(wcchlaBmodel)) 
ggplot(ChlSumMODEL, aes(y=WaterColChlA.ug.L+10, x=TreF, color=DayF))+
  stat_summary(position=position_dodge(width=.4))
wcchla_post_model<-posterior_samples(wcchlaBmodel)
wcchlaRRraw<-as_tibble(wcchla_post_model) %>%
  dplyr::select(-starts_with("r_Tank")) %>%
  mutate( #weekly response ratio between dead and control treatments
    rrDCd.20=`b_TreFDead:DayFM20`/`b_TreFControl:DayFM20`,
    rrDCd.15=`b_TreFDead:DayFM15`/`b_TreFControl:DayFM15`,
    rrDCd.3=`b_TreFDead:DayFM3`/`b_TreFControl:DayFM3`,
    rrDCd4=`b_TreFDead:DayF4`/`b_TreFControl:DayF4`,
    rrDCd11=`b_TreFDead:DayF11`/`b_TreFControl:DayF11`,
    rrDCd25=`b_TreFDead:DayF25`/`b_TreFControl:DayF25`,
    rrDCd39=`b_TreFDead:DayF39`/`b_TreFControl:DayF39`,
    # weekly response ratio between live and control treatments
    rrLCd.20=`b_TreFLive:DayFM20`/`b_TreFControl:DayFM20`,
    rrLCd.15=`b_TreFLive:DayFM15`/`b_TreFControl:DayFM15`,
    rrLCd.3=`b_TreFLive:DayFM3`/`b_TreFControl:DayFM3`,
    rrLCd4=`b_TreFLive:DayF4`/`b_TreFControl:DayF4`,
    rrLCd11=`b_TreFLive:DayF11`/`b_TreFControl:DayF11`,
    rrLCd25=`b_TreFLive:DayF25`/`b_TreFControl:DayF25`,
    rrLCd39=`b_TreFLive:DayF39`/`b_TreFControl:DayF39`,
    # weekly response ratio between dead and live treatments
    rrDLd.20=`b_TreFDead:DayFM20`/`b_TreFLive:DayFM20`,
    rrDLd.15=`b_TreFDead:DayFM15`/`b_TreFLive:DayFM15`,
    rrDLd.3=`b_TreFDead:DayFM3`/`b_TreFLive:DayFM3`,
    rrDLd4=`b_TreFDead:DayF4`/`b_TreFLive:DayF4`,
    rrDLd11=`b_TreFDead:DayF11`/`b_TreFLive:DayF11`,
    rrDLd25=`b_TreFDead:DayF25`/`b_TreFLive:DayF25`,
    rrDLd39=`b_TreFDead:DayF39`/`b_TreFLive:DayF39`,
    # total BACI response ratios for each comparison
    BACIdc=10^((rrDCd39+rrDCd25+rrDCd11+rrDCd4)/4)/10^((rrDCd.20+rrDCd.15+rrDCd.3)/3),
    BACIdl=10^((rrDLd39+rrDLd25+rrDLd11+rrDLd4)/4)/10^((rrDLd.20+rrDLd.15+rrDLd.3)/3),
    BACIlc=10^((rrLCd39+rrLCd25+rrLCd11+rrLCd4)/4)/10^((rrLCd.20+rrLCd.15+rrLCd.3)/3))
    #back transforming RR4
    #RR4dc = 10^`b_TreFDead:DayF4`/10^`b_TreFControl:DayF4`,
    #RR4dl = 10^`b_TreFDead:DayF4`/10^`b_TreFLive:DayF4`) 
# graphing weekly response ratios to investigate weekly response
wcchlaRRweekgraph <-wcchlaRRraw%>% 
  dplyr::select(starts_with("rr")) %>%
  gather(ratio,value) %>% 
  #easy graphing categories for fact_wrap
  mutate(ratio.Type=case_when(substr(ratio, 1,4)=="rrDC"~"Dead / Ctrl",
                              substr(ratio, 1,4)=="rrLC"~"Live / Ctrl",
                              substr(ratio, 1,4)=="rrDL"~"Dead / Live")) %>%
  mutate(ratio.F=factor(ratio, levels=c("rrDCd.20", "rrDCd.15", "rrDCd.3","rrDCd4",
                                        "rrDCd11","rrDCd25", "rrDCd39",
                                        "rrLCd.20", "rrLCd.15", "rrLCd.3","rrLCd4",
                                        "rrLCd11","rrLCd25", "rrLCd39",
                                        "rrDLd.20", "rrDLd.15", "rrDLd.3","rrDLd4",
                                        "rrDLd11","rrDLd25", "rrDLd39")))
# graphing results of BACI ratios
# evaluates if response increase (when rr >1) after impact
ggplot(wcchlaRRweekgraph, 
       aes(x=ratio.F, y=value))+
  geom_hline(yintercept=1, color="red", size=2)+
  geom_boxplot()+
  facet_wrap(~ratio.Type, scales="free_x")+
  scale_x_discrete(label=rep(c(-20,-15,-3,4,11,25,39),3), name="Sampling Day")+
  scale_y_continuous(name="WaterChlA Response Ratio")
wcchlaBACIgraph<-wcchlaRRraw %>% dplyr::select(starts_with("BACI")) %>%
  gather(ratio,value) %>% filter(!is.infinite(value), value < 10000)
ggplot(wcchlaBACIgraph, aes(x=ratio, y=value))+
  geom_boxplot()+
  scale_x_discrete(labels=c("BACI dead.ctrl","BACI dead.live","BACI live.ctrl"))

# quantifying BACI ratio - % prob of an increase, 50% increase 
1-ecdf(as.matrix(wcchlaBACIgraph[wcchlaBACIgraph$ratio=="BACIdc",2]))( 1 ) 
1-ecdf(as.matrix(wcchlaBACIgraph[wcchlaBACIgraph$ratio=="BACIdc",2]))( 1.5 )
qplot(as.matrix(wcchlaBACIgraph[wcchlaBACIgraph$ratio=="BACIdc",2]), 
      geom="histogram", xlim=c(0,5), bins=60,ylab="frequency",xlab="BACI d/c", 
      fill= as.matrix(wcchlaBACIgraph[wcchlaBACIgraph$ratio=="BACIdc",2]) >= 1)+
  scale_fill_manual(values=c("grey","black"), guide=F)

#### Benthic Chlorophyll a Analysis ####
summary(ChlSumMODEL$BenthicChlA.ug.cm)
bchlaBmodel<-brm((BenthicChlA.ug.cm+10)~TreF:DayFR+(1|TankF)-1, data=ChlSumMODEL,
                  family=Gamma(link='log'), 
                  prior=c(prior(normal(0,5), class="b"),
                          prior(gamma(0.01,0.01), class="shape")),
                  chains=4, iter=2000)
print(bchlaBmodel, prior=T) #negative beta predictions
plot(marginal_effects(bchlaBmodel))
ggplot(ChlSumMODEL, aes(y=BenthicChlA.ug.cm+10, x=TreF, color=DayFR))+
  stat_summary(position=position_dodge(width=.4))
bchla_post_model<-posterior_samples(bchlaBmodel)
bchlaRRraw<-as_tibble(bchla_post_model) %>%
  dplyr::select(-starts_with("r_Tank")) %>%
  mutate( #weekly response ratio between dead and control treatments
    rrDCd.20=`b_TreFDead:DayFRM20`/`b_TreFControl:DayFRM20`,
    rrDCd.15=`b_TreFDead:DayFRM15`/`b_TreFControl:DayFRM15`,
    rrDCd.3=`b_TreFDead:DayFRM3`/`b_TreFControl:DayFRM3`,
    rrDCd4=`b_TreFDead:DayFR4`/`b_TreFControl:DayFR4`,
    rrDCd11=`b_TreFDead:DayFR11`/`b_TreFControl:DayFR11`,
    rrDCd25=`b_TreFDead:DayFR25`/`b_TreFControl:DayFR25`,
    rrDCd39=`b_TreFDead:DayFR39`/`b_TreFControl:DayFR39`,
    # weekly response ratio between live and control treatments
    rrLCd.20=`b_TreFLive:DayFRM20`/`b_TreFControl:DayFRM20`,
    rrLCd.15=`b_TreFLive:DayFRM15`/`b_TreFControl:DayFRM15`,
    rrLCd.3=`b_TreFLive:DayFRM3`/`b_TreFControl:DayFRM3`,
    rrLCd4=`b_TreFLive:DayFR4`/`b_TreFControl:DayFR4`,
    rrLCd11=`b_TreFLive:DayFR11`/`b_TreFControl:DayFR11`,
    rrLCd25=`b_TreFLive:DayFR25`/`b_TreFControl:DayFR25`,
    rrLCd39=`b_TreFLive:DayFR39`/`b_TreFControl:DayFR39`,
    # weekly response ratio between dead and live treatments
    rrDLd.20=`b_TreFDead:DayFRM20`/`b_TreFLive:DayFRM20`,
    rrDLd.15=`b_TreFDead:DayFRM15`/`b_TreFLive:DayFRM15`,
    rrDLd.3=`b_TreFDead:DayFRM3`/`b_TreFLive:DayFRM3`,
    rrDLd4=`b_TreFDead:DayFR4`/`b_TreFLive:DayFR4`,
    rrDLd11=`b_TreFDead:DayFR11`/`b_TreFLive:DayFR11`,
    rrDLd25=`b_TreFDead:DayFR25`/`b_TreFLive:DayFR25`,
    rrDLd39=`b_TreFDead:DayFR39`/`b_TreFLive:DayFR39`,
    # total BACI response ratios for each comparison
    BACIdc=10^((rrDCd39+rrDCd25+rrDCd11+rrDCd4)/4)/10^((rrDCd.20+rrDCd.15+rrDCd.3)/3),
    BACIdl=10^((rrDLd39+rrDLd25+rrDLd11+rrDLd4)/4)/10^((rrDLd.20+rrDLd.15+rrDLd.3)/3),
    BACIlc=10^((rrLCd39+rrLCd25+rrLCd11+rrLCd4)/4)/10^((rrLCd.20+rrLCd.15+rrLCd.3)/3)) 
# graphing weekly response ratios to investigate weekly response
bchlaRRweekgraph <-bchlaRRraw%>% 
  dplyr::select(starts_with("rr")) %>%
  gather(ratio,value) %>% 
  #easy graphing categories for fact_wrap
  mutate(ratio.Type=case_when(substr(ratio, 1,4)=="rrDC"~"Dead / Ctrl",
                              substr(ratio, 1,4)=="rrLC"~"Live / Ctrl",
                              substr(ratio, 1,4)=="rrDL"~"Dead / Live")) %>%
  mutate(ratio.F=factor(ratio, levels=c("rrDCd.20", "rrDCd.15", "rrDCd.3","rrDCd4",
                                        "rrDCd11","rrDCd25", "rrDCd39",
                                        "rrLCd.20", "rrLCd.15", "rrLCd.3","rrLCd4",
                                        "rrLCd11","rrLCd25", "rrLCd39",
                                        "rrDLd.20", "rrDLd.15", "rrDLd.3","rrDLd4",
                                        "rrDLd11","rrDLd25", "rrDLd39")))
# graphing results of BACI ratios
# evaluates if response increase (when rr >1) after impact
ggplot(bchlaRRweekgraph, 
       aes(x=ratio.F, y=value))+
  geom_hline(yintercept=1, color="red", size=2)+
  geom_boxplot()+
  facet_wrap(~ratio.Type, scales="free_x")+
  scale_x_discrete(label=rep(c(-20,-15,-3,4,11,25,39),3), name="Sampling Day")+
  scale_y_continuous(name="Benthic ChlA Response Ratio")
bchlaBACIgraph<-bchlaRRraw %>% dplyr::select(starts_with("BACI")) %>%
  gather(ratio,value)
ggplot(bchlaBACIgraph, aes(x=ratio, y=value))+geom_boxplot()+
  scale_x_discrete(labels=c("BACI dead.ctrl","BACI dead.live","BACI live.ctrl"))

# quantifying BACI ratio - % prob of an increase, 50% increase 
1-ecdf(as.matrix(bchlaBACIgraph[bchlaBACIgraph$ratio=="BACIdc",2]))( 1 ) 
bchlaBACIgraph %>% group_by(ratio) %>% dplyr::summarize(meanBACI=mean(value))
1-ecdf(as.matrix(bchlaBACIgraph[bchlaBACIgraph$ratio=="BACIdc",2]))( 1.16 )
qplot(as.matrix(bchlaBACIgraph[bchlaBACIgraph$ratio=="BACIdc",2]), 
      geom="histogram", xlim=c(0,5), bins=60,ylab="frequency",xlab="BACI d/c", 
      fill= as.matrix(bchlaBACIgraph[bchlaBACIgraph$ratio=="BACIdc",2]) >= 1)+
  scale_fill_manual(values=c("grey","black"), guide=F)

#### metabolism & decomposition ####
head(Metstats)
Metstats1 <- Metstats %>%
  mutate(DayF=factor(Day, levels=c(4,11,25,39)),
         TreF=factor(NewTreat, levels=c("Control","Live","Dead")),
         TankF=factor(Tank))
MetModel<-brm(meanGPP~TreF*Day+(1|TankF)-1, data=Metstats1,
      #family=Gamma(link='log'), 
      #prior=c(prior(normal(0,5), class="b"),
      #        prior(gamma(0.01,0.01), class="shape")),
      chains=4, iter=2000)
MetModel.tankonly<-brm(meanGPP~(1|TankF), data=Metstats1,
              #family=Gamma(link='log'), 
              #prior=c(prior(normal(0,5), class="b"),
              #        prior(gamma(0.01,0.01), class="shape")),
              chains=4, iter=2000)
print(MetModel, prior=T) 
plot(marginal_effects(MetModel))
ggplot(Metstats1, aes(x=TreF, y=meanGPP, color=DayF)) +
  stat_summary(position=position_dodge(1))
bayes_factor(MetModel)
