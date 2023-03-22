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

#Sampling years prior to 1984 (data begin in 1982) were excluded from analysis due
#to large apparent increases in the number of species recorded in the first two years. 
#(Batt et al. 2017)

EBS <-  EBS[year>=1984,]

#Months 6,7,8 sampled the most, exclude other months

EBS <-  EBS[month %in% c(6,7,8),]

#Match lat/lon sampling points to hexagonal cells, so that we can see how many cells to keep to maintain a lot of observation points

EBS.cells <- data.table()

#unique lat lon combos
EBS.latlon <- unique(EBS[,.(longitude,latitude,haul_id,year)])
  
#coordinates to simple feature points
EBS.sf <- st_as_sf(EBS.latlon, coords = c("longitude","latitude"), crs = 4326)
  
#sf to sp (spatial points)
EBS.sp <- as(EBS.sf, "Spatial")

#change projection  
proj4string(EBS.sp) <- CRS("+proj=longlat")
  
#projection
proj <- "+proj=robin +lon_0=-140 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
  
#work with units
EBS.sp <- spTransform(EBS.sp, CRS(proj)) #note km^2 units
  
#use Concaveman to convert points to polygon
EBS.polygon  <- concaveman(EBS.sf, 2, 3)
  
EBS.polygon_spapol <- as(EBS.polygon, "Spatial") #convert simple polygon feature to spatial polygon
  
proj4string(EBS.polygon_spapol) <- CRS("+proj=longlat")
  
EBS.polygon_spapol.p <- spTransform(EBS.polygon_spapol, CRS(proj)) #note km^2 units
  
#create grid 

#set cell size
cell_area  <- 7774.2 #km2 (8 from dggrdr)

#calculate cell_diameter of hexagons from cell_areas
cell_diameter_km <- sqrt(2 * cell_area / sqrt(3)) # in meters
  
  
ext <- as(extent(EBS.polygon_spapol.p)
            + 2*cell_diameter_km #add a buffer to make sure all observations are assigned a cell
            , "SpatialPolygons")
  # plot(ext)
  # plot(sp.p, add = T, pch = ".")
  
projection(ext) <- projection(EBS.polygon_spapol.p) #match projection
  
# generate array of hexagon centers
g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter_km, offset = c(0.5, 0.5))
  
# convert center points to hexagons
g <- HexPoints2SpatialPolygons(g, dx = cell_diameter_km)
  
#visualize
plot(g)
plot(EBS.polygon_spapol.p, add = T, pch = ".")

#link lat lon to cell#
#where do they overlap
#over(x=location of queries, y = layer from which geometries are queried)
EBS.sp$cell_ID <- sp::over(EBS.sp, g) 
  
#link lat long to cell #s
EBS.latlon[,cell_ID := EBS.sp$cell_ID][,cell_year_count := .N, .(cell_ID, year)]
  
#link back to subsetted database of observations
EBS <- EBS[EBS.latlon, on = c("longitude", "latitude","year","haul_id")]

#unique year, #cells  sampled
year_cells_sampled <- unique(EBS[,.(year,cell_ID)])
year_cells_sampled <- year_cells_sampled[,yearly_cell_count := .N,year]
year_cells_sampled <- unique(year_cells_sampled[,.(year,yearly_cell_count)])

#eliminate years with less than 70% of cells sampled
year_benchmark <- 0.70
benchmark_value <- year_benchmark*max(year_cells_sampled[,yearly_cell_count])

#only keep years where over 70% of cells are sampled
year_cells_sampled[,benchmark := yearly_cell_count >= benchmark_value]

years_deleted <- unique(year_cells_sampled[benchmark == F,year]) #which years are left out?

years_kept <-unique(year_cells_sampled[benchmark ==T,year]) #which years to keep

#identify any cells that are not sampled in 70% of years
EBS[,year_cell_count := length(unique(haul_id)),.(year,cell_ID)] # unique haul ids per cell per year

cell_by_year <- unique(EBS[, .(cell_ID,year)])

cell_by_year[,years_per_cell := .N,cell_ID]

#cell ids to remove and keep
#in any year, which cells are sampled in less than 70% of years
#we'll make benchmark 70% just for now
cell_benchmark <- 0.70
benchmark_value_year_count <- cell_benchmark*max(cell_by_year[,years_per_cell])

cell_id_remove <- unique(cell_by_year[years_per_cell<benchmark_value_year_count,cell_ID])

cells_deleted_percent <- round(length(cell_id_remove)/length(unique(EBS[,cell_ID]))*100,1)

#reduce to cells that are well sampled
EBS.r <- EBS[!(cell_ID %in% cell_id_remove),]

#What percent of hauls does this remove?
hauls_removed_yearcell <- round((length(unique(EBS[,haul_id]))-length(unique(EBS.r[,haul_id])))/length(unique(EBS[,haul_id]))*100,1) 

#removes 0.1% of cells

#What percent of observations does this remove?
obs_removed_yearcell <- round((nrow(EBS)-nrow(EBS.r))/nrow(EBS)*100,1)

#removes 0.1% of observations

#######################
##Save finalized dataset
fwrite(EBS.r, file.path("Data","EBS.r.csv"))


