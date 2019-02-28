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

OKflowLines <- getFlowLines(c(basinsOK@bbox[1], basinsOK@bbox[3],
                              basinsOK@bbox[2],basinsOK@bbox[4]),
                              streamorder = 4)
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

usa <- map_data("usa")
states <- map_data("state")
oktidy<-broom::tidy(OKriver)
oklakestidy<-broom::tidy(oklakes)
myMap <- get_map(location=OKriver@bbox,source="stamen",
                 maptype="terrain")
ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="darkgrey")+
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill=NA, color="black") +
  coord_fixed(1.3)+theme_bw()
ggmap(myMap)+
  geom_path(data = oktidy, aes(x = long, y = lat, group=group))+
  geom_polygon(data=oklakestidy, aes(x=long,y=lat,group=group), fill="black")
