###########################
#This script labels latlon tows as outer, middle, or inner domain in EBS


#Inshore of the 50 m isobath, the coastal domain is vertically homogeneous and 
#separated from the adjacent middle domain by a narrow ('"10 km) front. Between 
#the 50 m and 100 m isobaths, the middle domain tends toward a strongly 
#stratified two-layered structure, and is separated from the adjacent outer 
#domain by a weak front. Between the 100 m isobath and the shelf break ('" 1 70 m depth),
#the outer domain has surface and bottom mixed layers above and below a stratified interior.
#Kinder and Schumacher "Hydrographic Structure Over the Continental Shelf of the Southeastern Bering Sea"

#######################
##VERSIONS##
#######################
#R 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
#macOS Big Sur 11.7

#######################
##PACKAGES##
#######################
library(data.table)
library(sf)
library(raster)
library(ggplot2)

#######################
##DATA##
#######################

#Pull in cleaned EBS data from 1_PullEBS_data.R

EBS.r <- fread(file.path("Data","EBS.r.csv"))

#######################
##DATA PROCESSING
#######################

#visualize depth distribution
hist(EBS.r$depth)
abline(v = 50, col = "orange")
abline(v = 100, col = "orange")

#unique points full region
EBS.latlon <- unique(EBS.r[,.(longitude, latitude, depth)])

#new column labeling depending on domain
EBS.domains.r <- EBS.latlon[,domain := ifelse(depth <= 50,"Inner",ifelse(depth > 100,"Outer", "Middle"))]

#create single polygon feature
EBS.domains.sf <- st_as_sf(EBS.domains.r, coords=c("longitude","latitude"), crs = 4326)

# Aggregate - which unions the points into MULTIPOINTS (default do_union=TRUE)
# and uses FUN to assign the first value in other fields
EBS.domains.ag  <- st_sf(
  aggregate(
    EBS.domains.sf,
    by=list(ID=EBS.domains.sf$domain),
    FUN=function(vals){vals[1]}))

#transform projection
EBS.domains.t <- st_transform(EBS.domains.ag,
                         crs=3467)

#convex hull
EBS.domains.chull <- EBS.domains.ag
st_geometry(EBS.domains.chull) <- st_convex_hull(EBS.domains.ag$geometry)

plot(EBS.domains.chull['domain'])
plot(EBS.domains.ag['domain'], pch = 16)

#make map of Alaska plus points colored by domain

#create data frame of east coast CA and USA from rnaturalearth
state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")

alaska <- st_crop(state_prov, xmin = -185, ymin = 50, xmax = -140, ymax = 70)

#change projections
alaska.t <- st_transform(alaska,
                             crs=3467)
  
#(Warning) arises in recent versions of sp due to an update from PROJ4 to PROJ6.
#In the meantime, the warning is just a nuisance and has no implications (Jim Thorson, 2020)

#basemap
Alaska_domains <- ggplot() + 
  geom_sf(data = alaska.t, color="white", fill="gainsboro", linewidth = 0.1) + 
  geom_sf(data = EBS.domains.sf, aes(color = domain), alpha = 0.5, size = 0.1) +
  labs(color = "Domain") +
  scale_x_continuous(breaks=seq(-180,-150, 10), labels = c("180˚W","170˚W","160˚W","150˚W")) +
  scale_y_continuous(breaks=seq(52,65,5)) +
  coord_sf(datum = sf::st_crs(4326), expand = F) +
  scale_color_manual(values = c("#AA4499","#44AA99","#999933"), labels = c("Inner\n(to 50m)","Middle\n(to 100m)","Outer")) +
  theme_classic() +
  #make legend visible
  guides(color = guide_legend(override.aes = list(size=4, alpha = 0.8)))

Alaska_domains

ggsave(Alaska_domains, path = file.path("Figures"),
       filename = "Alaska_domains_map.jpg", height = 4, width = 5, units = "in")




