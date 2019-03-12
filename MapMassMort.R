# idk what libraries are actually in use anymore
library(maps);library(maptools);library(calibrate);library(sp);
library(rgdal);library(hydroMap);library(dataRetrieval); library(raster)
library(rgeos);library(ggplot2);library(ggmap);library(mapdata)
#usgs sites (used to grab water basins)
OKsites<-c("07339000", "07338500","07335790","07336200","07336200")
GAsites<- c("02357000","02353500","02354500","02350900",
            "02352500","02356000","02353000")
test<-c("02358000") #bottom of flint river - needs to be cropped
basinsOK <- getBasin(OKsites) #grab water basins in oklahoma
basinsOK@bbox #bounding box of water basins
basinsGA <- getBasin(GAsites) #grab water basins in georgia
basinsGA@bbox #bounding box of water basins
#cropping flint water basin (too large)
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
#binding basins into one spatial polygon
basins1<-spRbind(basinsOK, basinsGA) 
basins<-spRbind(basins1,bottomFlint)

#grab all streams 4 order and above within the georgia bounding box
GAflowLines <- getFlowLines(c(-85.2, -83.6,
                              30.69960,basinsGA@bbox[4]),
                            streamorder = 4)
#pull the specific rivers mentioned in Golladay et al 2004
SC<-GAflowLines[GAflowLines@data$gnis_name=="Spring Creek",]
IC<-GAflowLines[GAflowLines@data$gnis_name=="Ichawaynochaway Creek",]
CC<-GAflowLines[GAflowLines@data$gnis_name=="Chickasawhatchee Creek",]
KC<-GAflowLines[GAflowLines@data$gnis_name=="Kinchafoonee Creek",]
MC<-GAflowLines[GAflowLines@data$gnis_name=="Muckalee Creek",]
PC<-GAflowLines[GAflowLines@data$gnis_name=="Pachitla Creek",]
Flint<-GAflowLines[GAflowLines@data$gnis_name=="Flint River",]
#bind them together to make one Spatial Lines of the Flint Rivers 
GAr1<-spRbind(SC,IC)
GAr2<-spRbind(GAr1,CC)
GAr3<-spRbind(GAr2,KC)
GAr4<-spRbind(GAr3,MC)
GAr5<-spRbind(GAr4,PC)
GArivers<-spRbind(GAr5,Flint) 
galakes<-crop(shape_hydropoly, c(basinsGA@bbox[1], basinsGA@bbox[3],
                        30.69960,31.4)) #pull out that big lake
#plot to check if it matches Golladay et al 2004
plot(GArivers, ylim=c(30.69960,32.5))
plot(galakes, add=T, col='lightblue')
plot(getBasin(test), add=T, border='darkgrey')

#pull stream orders 8 and above for continental US
#takes about 8 minutes
USflow8 <- getFlowLines(c(-124.6813, -67.007, 
                               25.1299,49.3832),
                             streamorder = 8)
#plot(USflowlinese)

#get the Oklahoma rivers stream order 4 and above
OKflowLines <- getFlowLines(c(basinsOK@bbox[1], basinsOK@bbox[3],
                              33.75,34.75),
                            streamorder = 4)
#pull out rivers of interest
Kiam<-OKflowLines[OKflowLines@data$gnis_name=="Kiamichi River",]
LR<-OKflowLines[OKflowLines@data$gnis_name=="Little River",]
Mt<-OKflowLines[OKflowLines@data$gnis_name=="Mountain Fork",]
OKr1<-spRbind(Kiam,LR)
OKriver<-spRbind(OKr1,Mt) #bind them into one Spatial Lines file
#pull out lakes (these look ugly :( )
oklakes<-crop(shape_hydropoly, c(basinsOK@bbox[1], basinsOK@bbox[3],
                                 basinsOK@bbox[2],basinsOK@bbox[4]))
#check I have what I need
plot(OKriver)
plot(oklakes, add=T, col='lightblue')
plot(basinsOK, add=T, border='darkgrey')
plot(OKflowLines[OKflowLines@data$gnis_name=="Red River",], add=T, col='lightgrey')

library(ggmap)
#these calls 'tidy' the data so I can plot it with ggplot
usa <- map_data("usa")
states <- map_data("state")
gatidy<-broom::tidy(GArivers) #lots of warnings about character strings
GAflow<-broom::tidy(GAflowLines) #lots of warnings about character strings
oktidy<-broom::tidy(OKriver) #lots of warnings about character strings
OKflow<-broom::tidy(OKflowLines) #lots of warnings about character strings
oklakestidy<-broom::tidy(oklakes)
galakestidy<-broom::tidy(galakes)
okbasins<-broom::tidy(basinsOK) #lots of warnings about character strings
gabasins<-broom::tidy(basinsGA)
canada<-map_data("world","canada")
mexico<-map_data("world",'mexico')
usstreams<-broom::tidy(USflow8) #lots of warnings about character strings

#make the continental map to show the distance between our two study sites
bigmap<-ggplot() + 
  geom_polygon(data=states, aes(x=long, y=lat, group=group),
               fill="white", color="darkgrey")+
  geom_polygon(data=canada, aes(x=long, y=lat, group=group),
               fill="darkgreen", color="darkgrey")+
  geom_polygon(data=mexico, aes(x=long, y=lat, group=group), 
               fill="darkgreen", color="darkgrey")+
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
               fill=NA, color="black") +
  geom_polygon(data = okbasins, aes(x=long, y = lat, group = group), 
               fill="black", color="black") +
  geom_polygon(data = gabasins, aes(x=long, y = lat, group = group),
               fill="black", color="black") +
  geom_path(data = usstreams, aes(x=long, y = lat, group = group), 
            color="lightblue", size=1.3) +
  geom_text(aes(x=c(-97,-82.5), y=c(34.75,32.55)), 
            label=c("A","B"))+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim= c(-125.6813, -66.007), ylim=c(24.1299,50.3832))+theme_bw()
okmap<-ggplot()+
  geom_path(data=OKflow, aes(x=long, y=lat, group=group), size=1.3, 
            color="lightblue")+
  geom_path(data = oktidy, aes(x = long, y = lat, group=group), size=1.3)+
  geom_polygon(data=oklakestidy, aes(x=long,y=lat,group=group), 
               fill="lightblue")+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="darkgrey")+
  geom_point(aes(x=c(-95.354508,-95.062153), y=c(34.574689,34.6569991)), 
             color="black", fill="white",size=3, shape=24)+
  geom_text(aes(x=c(-95.33,-95.01), y=c(34.522,34.611)), 
          label=c("K3","K2"))+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim=c(-95.7, -94.3),ylim=c(33.8,34.75))+theme_bw()
gamap<-ggplot()+
  geom_path(data=GAflow, aes(x=long, y=lat, group=group), size=1.3, 
            color="lightblue")+
  geom_path(data = gatidy, aes(x = long, y = lat, group=group), size=1.25)+
  geom_polygon(data=galakestidy, aes(x=long,y=lat,group=group), fill="lightblue", 
               size=1.2)+
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill=NA, color="darkgrey")+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim=c(-85, -83.65),ylim=c(30.70,32.3))+theme_bw()
BigMapF<-grid.arrange(bigmap, top="Continental USA")
sites<-plot_grid(okmap,gamap, labels=c("A","B"), nrow=1)
mapmatrix<-plot_grid(BigMapF,sites, ncol=1)
ggsave("test.tiff",mapmatrix)

#add labels
#add cities?
library(gridExtra)
p1<-arrangeGrob(, okmap, gamap, layout_matrix = rbind(c(1,1),c(2,3)),
                respect=TRUE)
ggsave("testmap.tiff",p1)

library(LifeTables)
data(MLTobs) 
test.mx.m <- mlt.mx[,1]
# build the life table 
lt.mx(nmx=test.mx.m, sex="male")
