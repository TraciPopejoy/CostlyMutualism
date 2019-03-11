library(maps);library(maptools);library(calibrate);library(sp);
library(rgdal);library(hydroMap);library(dataRetrieval); library(raster)
library(rgeos);library(ggplot2);library(ggmap);library(mapdata)

OKsites<-c("07339000", "07338500","07335790","07336200","07336200")
GAsites<- c("02357000","02353500","02354500","02350900",
            "02352500","02356000","02353000")
test<-c("02358000")
basinsOK <- getBasin(OKsites)
basinsOK@bbox
basinsGA <- getBasin(GAsites)
basinsGA@bbox
bigbasin <- crop(getBasin(test), extent(-85, -83.40847, 30.69960, 31.41))
bottomFlint <- erase(bigbasin,basinsGA)

#from StackExchange: 
#https://gis.stackexchange.com/questions/32732/proper-way-to-rbind-spatialpolygonsdataframes-with-identical-polygon-ids
makeUniform<-function(SPDF){
  pref<-substitute(SPDF)  #just putting the file name in front.
  newSPDF<-spChFIDs(SPDF,as.character(paste(pref,rownames(as(SPDF,"data.frame")),sep="_")))
  return(newSPDF)
}
basinsOK<-makeUniform(basinsOK)
basinsGA<-makeUniform(basinsGA)

## Plot the output
basins1<-spRbind(basinsOK, basinsGA)
basins<-spRbind(basins1,bottomFlint)

plot(basins)
map('usa')
map('state', add=T)
map('state', region=c('oklahoma','georgia'),col='lightgrey', fill=T,add=T)
map(basins, fill=T, col='black', add=T)
map('state',
    region = c('texas','oklahoma','arkansas','louisiana',
                        'mississippi','alabama','tennessee', 'georgia',
                       'north carolina','south carolina',
                        'florida'), 
    lwd=2.5, col="darkgray") 
map('state',
    xlim=c((basins@bbox[1]-2),(basins@bbox[3]+2)),
    ylim=c(basins@bbox[2]-2,basins@bbox[4]+2),
    lwd=2.5, col="darkgray") 
plot(basins,add=T, col='black')
plot(OKriver, add=T)
plot(GArivers, add=T)

GAflowLines <- getFlowLines(c(basinsGA@bbox[1], basinsGA@bbox[3],
                              30.69960,basinsGA@bbox[4]),
                            streamorder = 4)
SC<-GAflowLines[GAflowLines@data$gnis_name=="Spring Creek",]
IC<-GAflowLines[GAflowLines@data$gnis_name=="Ichawaynochaway Creek",]
CC<-GAflowLines[GAflowLines@data$gnis_name=="Chickasawhatchee Creek",]
KC<-GAflowLines[GAflowLines@data$gnis_name=="Kinchafoonee Creek",]
MC<-GAflowLines[GAflowLines@data$gnis_name=="Muckalee Creek",]
Flint<-GAflowLines[GAflowLines@data$gnis_name=="Flint River",]
GAr1<-spRbind(SC,IC)
GAr2<-spRbind(GAr1,CC)
GAr3<-spRbind(GAr2,KC)
GAr4<-spRbind(GAr3,MC)
GArivers<-spRbind(GAr4,Flint)
galakes<-crop(shape_hydropoly, c(basinsGA@bbox[1], basinsGA@bbox[3],
                        30.69960,31.4))
plot(GArivers, ylim=c(30.69960,32.5))
plot(galakes, add=T, col='lightblue')
plot(getBasin(test), add=T, border='darkgrey')
GAflowLines <- getFlowLines(c(basinsGA@bbox[1], basinsGA@bbox[3],
                              30.69960,basinsGA@bbox[4]),
                            streamorder = 4)

USflowlinese <- getFlowLines(c(-124.6813, -67.007,
                             25.1299,49.3832),
                              streamorder = 7)
USflow8 <- getFlowLines(c(-124.6813, -67.007,
                               25.1299,49.3832),
                             streamorder = 8)
plot(USflowlinese)
Kiam<-OKflowLines[OKflowLines@data$gnis_name=="Kiamichi River",]
LR<-OKflowLines[OKflowLines@data$gnis_name=="Little River",]
Mt<-OKflowLines[OKflowLines@data$gnis_name=="Mountain Fork",]
OKr1<-spRbind(Kiam,LR)
OKriver<-spRbind(OKr1,Mt)
oklakes<-crop(shape_hydropoly, c(basinsOK@bbox[1], basinsOK@bbox[3],
                                 basinsOK@bbox[2],basinsOK@bbox[4]))

plot(OKriver)
plot(oklakes, add=T, col='lightblue')
plot(basinsOK, add=T, border='darkgrey')
plot(OKflowLines[OKflowLines@data$gnis_name=="Red River",], add=T, col='lightgrey')

require(devtools)
devtools::install_github("dkahle/ggmap", ref = "tidyup")

usa <- map_data("usa")
states <- map_data("state")
gatidy<-broom::tidy(GArivers)
oktidy<-broom::tidy(OKriver)
oklakestidy<-broom::tidy(oklakes)
galakestidy<-broom::tidy(galakes)
okbasins<-broom::tidy(basinsOK)
gabasins<-broom::tidy(basinsGA)
canada<-map_data("world","canada")
mexico<-map_data("world",'mexico')
usstreams<-broom::tidy(USflow8)

bigmap<-ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="white", color="darkgrey")+
  geom_polygon(data=canada, aes(x=long, y=lat, group=group), fill="darkgreen", color="darkgrey")+
  geom_polygon(data=mexico, aes(x=long, y=lat, group=group), fill="darkgreen", color="darkgrey")+
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill=NA, color="black") +
  geom_polygon(data = okbasins, aes(x=long, y = lat, group = group), fill="goldenrod2", color="black") +
  geom_polygon(data = gabasins, aes(x=long, y = lat, group = group), fill="goldenrod2", color="black") +
  geom_path(data = usstreams, aes(x=long, y = lat, group = group), color="lightblue", 
            size=1.3) +
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim= c(-125.6813, -66.007), ylim=c(24.1299,50.3832))+theme_bw()
okmap<-ggplot()+
  geom_path(data = oktidy, aes(x = long, y = lat, group=group), size=1.3)+
  geom_polygon(data=oklakestidy, aes(x=long,y=lat,group=group), fill="lightblue")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="darkgrey")+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim=c(-95.75, -94.25),ylim=c(33.75,34.75))+theme_bw()
gamap<-ggplot()+
  geom_path(data = gatidy, aes(x = long, y = lat, group=group), size=1.25)+
  geom_polygon(data=galakestidy, aes(x=long,y=lat,group=group), fill="lightblue")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="darkgrey")+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim=c(-85.5, -83.5),ylim=c(30.5,33.5))+theme_bw()

#add smaller rivers to map
#add labels
#add cities?
p1<-arrangeGrob(bigmap,okmap,gamap, layout_matrix = rbind(c(1,1),c(2,3)), 
            respect=TRUE)
ggsave("testmap.tiff", p1)

library(LifeTables)
data(MLTobs) 
test.mx.m <- mlt.mx[,1]
# build the life table 
lt.mx(nmx=test.mx.m, sex="male")
