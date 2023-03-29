###########################
#This script loads data from FishGlob specifically for EBS
#The FishGlob cleaned and standardized trawl data is currently under 
#However, these data are already available online at this link: https:, , www.fisheries.noaa.gov, foss, f?p=215%3A28

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
library(sp)
library(concaveman)
library(raster)
library(ggplot2)

#######################
##DATA##
#######################

#Pull in compiled and cleaned data from FishGlob, downloaded on November 28, 2022 (V 1.5) at 
#(https:, , drive.google.com, file, d, 1MgXKhmIufUtjE_Y_mWpDeldJKfd7nDEM, view?usp=share_link). 
#This is typically compiled by Dr. Aurore Maureaud.

FishGlob_1.5 <- fread("/Users/zoekitchel/Documents/grad_school/Rutgers/Repositories/trawl_spatial_turnover_git/data/FISHGLOB_v1.5_clean.csv")

#######################
##DATA PROCESSING
#######################

##Restrict to EBS only
EBS <- FishGlob_1.5[survey_unit == "EBS",]

rm(FishGlob_1.5)

##Clean and Standardize to maintain sample footprint through time
###Because EBS using stations, it makes the most sense to do data standardization and Q/C using stations instead of using arbitrary grid cells that make more sense when you're doing multiple surveys at once

#Sampling years prior to 1984 (data begin in 1982) were excluded from analysis due
#to large apparent increases in the number of species recorded in the first two years. 
#(Batt et al. 2017)

EBS <-  EBS[year>=1984,]

#Months 6,7,8 sampled the most, exclude other months

EBS <-  EBS[month %in% c(6,7,8),]

#unique lat lon combos
EBS.latlon <- unique(EBS[,.(longitude,latitude,haul_id,year)])

#unique year, #stations  sampled
year_stations_sampled <- unique(EBS[,.(year,station)])
year_stations_sampled <- year_stations_sampled[,yearly_station_count := .N,year]
year_stations_sampled <- unique(year_stations_sampled[,.(year,yearly_station_count)])

#eliminate years with less than 75% of stations sampled
year_benchmark <- 0.75
benchmark_value <- year_benchmark*max(year_stations_sampled[,yearly_station_count])

#only keep years where over 75% of stations are sampled
year_stations_sampled[,benchmark := yearly_station_count >= benchmark_value]

years_deleted <- unique(year_stations_sampled[benchmark == F,year]) #which years are left out?

years_kept <-unique(year_stations_sampled[benchmark ==T,year]) #which years to keep

#identify any stations that are not sampled in 75% of years
EBS[,year_station_count := length(unique(haul_id)),.(year,station)] # unique haul ids per station per year

station_by_year <- unique(EBS[, .(station,year)])

station_by_year[,years_per_station := .N,station]

#stations to remove and keep
#in any year, which stations are sampled in less than 90% of years
station_benchmark <- 0.90
benchmark_value_year_count <- station_benchmark*max(station_by_year[,years_per_station])

station_remove <- unique(station_by_year[years_per_station < benchmark_value_year_count,station])

#where are these stations?
removed_stations <- EBS[station %in% station_remove,.(latitude, longitude, year, station)]
plot(EBS.latlon$longitude, EBS.latlon$latitude, pch = ".")
points(removed_stations$longitude, removed_stations$latitude, col = "pink")

#make map of Alaska plus stations, highlighting excluded stations

#create simple feature all stations
EBS.latlon.sf <- st_as_sf(EBS.latlon, coords=c("longitude","latitude"), crs = 4326)

#create simple feature excluded stations
removed_stations.sf <- st_as_sf(removed_stations, coords=c("longitude","latitude"), crs = 4326)

#create data frame of east coast CA and USA from rnaturalearth
state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")

alaska <- st_crop(state_prov, xmin = -185, ymin = 50, xmax = -140, ymax = 70)

#change projections
alaska.t <- st_transform(alaska,
                         crs=4326)

#basemap
EBS_stations <- ggplot() + 
  geom_sf(data = alaska.t, color="white", fill="gainsboro", linewidth = 0.1) + 
  geom_sf(data = removed_stations.sf, color = "pink", size = 1) +
  geom_sf(data = EBS.latlon.sf, size = 0.01) +
   coord_sf(expand = T, xlim = c(-180,-140), datum = sf::st_crs(4326)) +
   theme_classic()

EBS_stations

#save map
ggsave(EBS_stations, path = file.path("Figures", "Supplement"),
       filename = "EBS_stations.jpg", height = 4, width = 4, units = "in")

#what percent stations deleted

stations_deleted_percent <- round(length(station_remove)/length(unique(EBS[,station]))*100,1)

#removes 4.5% of stations

#reduce to stations that are well sampled
EBS.r <- EBS[!(station %in% station_remove),]

#What percent of hauls does this remove?
hauls_removed_yearstation <- round((length(unique(EBS[,haul_id]))-length(unique(EBS.r[,haul_id])))/length(unique(EBS[,haul_id]))*100,1) 


#What percent of observations does this remove?
obs_removed_yearstation <- round((nrow(EBS)-nrow(EBS.r))/nrow(EBS)*100,1)

#removes 4.1% of observations

#######################
##Save finalized dataset
fwrite(EBS.r, file.path("Data","EBS.r.csv"))


