library(dplyr)
library(data.table)
library(raster)
library(sp)
library(lubridate)
library(ggplot2)

## Script loads and extracts CDEC Snow Water equivalents nearest
## each lake in sierra lakes dataset for years 2004-2019 on April 1 

# Load in lake data and CDEC data
lake_data = fread("../data/archival/sites_nov2020.csv")
cdec_stations = fread("../data/environmental_covariates/snowpack_data/cdec_data_from_dozier/LocationStation.csv")
swe_dat = fread("../data/environmental_covariates/snowpack_data/cdec_data_from_dozier/SnowPillow_SWE.csv")
swe_dat[, MeasDate:=as.POSIXct(MeasDate, format="%Y-%m-%d")][, c("year", "month", "day"):=.(year(MeasDate), month(MeasDate), day(MeasDate))]
cdec_w_dat = unique(swe_dat$CDEC)
cdec_stations = cdec_stations[CDEC %in% cdec_w_dat]

# Set projection strings
crs_latlon = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
epsg_11N = "EPSG:32611"

# Transform to latlong
lake_trans = SpatialPoints(lake_data[, .(utme, utmn)], proj4string=CRS(epsg_11N)) %>%
	   		 spTransform(CRS(crs_latlon))
lake_coords = lake_trans@coords
cdec_spatial = SpatialPoints(cdec_stations[, .(Longitude, Latitude)], proj4string=CRS(epsg_11N))

# Check overlay...Looks good
plot(cdec_spatial, col='red')
points(lake_trans)

#### Find which CDEC stations are closest to which lakes ####

nearest_cdec = array(NA, dim=nrow(lake_data))
dists_array = array(NA, dim=nrow(lake_data))
for(i in 1:nrow(lake_data)){
	dists = spDistsN1(as.matrix(cdec_stations[, .(Longitude, Latitude)]), lake_coords[i, ], longlat=TRUE)
	min_ind = which.min(dists)
	nearest_cdec[i] = cdec_stations[min_ind, ]$CDEC
	dists_array[i] = min(dists)
}

lake_data$CDEC = nearest_cdec
lake_data$min_dist_km = dists_array

# Add elevations
lake_data = merge(lake_data, cdec_stations[, .(CDEC, Elevation)], by="CDEC")
colnames(lake_data)[length(colnames(lake_data))] = "elevation_cdec"
fwrite(lake_data, "../data/environmental_covariates/lake_cdec_distances.csv")


##### Extract the CDEC SWE data on April 1 ####

swe_dat = fread("../data/environmental_covariates/snowpack_data/cdec_data_from_dozier/SnowPillow_SWE.csv")
swe_dat[, MeasDate:=as.POSIXct(MeasDate, format="%Y-%m-%d")][, c("year", "month", "day"):=.(year(MeasDate), month(MeasDate), day(MeasDate))]

# Get a sense of what the april 1 date is looking like
test_dat = swe_dat[CDEC == "ADM"]
april_test_dat = test_dat[month == 4 & day == 1]
ggplot(test_dat) + geom_line(aes(x=MeasDate, y=SWEmm)) +
				   geom_point(data=april_test_dat, aes(x=MeasDate, y=SWEmm), color='red')

# Get all the data on April 1 of each year
april_dat = swe_dat[month == 4 & day == 1 & year >= 2000]
temp = merge(april_dat, cdec_stations[, .(Latitude, Longitude, Elevation, CDEC)], by="CDEC")

plot(temp$Elevation, temp$SWEmm)

# Look at the year effect 
# The variability in year is greater than the variability in space
ggplot(temp) + geom_line(aes(x=year, y=SWEmm, group=CDEC))

# Make a full lake data set for all years 
all_lakes_years = list()
for(y in 2000:2020){
	tlake = lake_data[, .(id, CDEC)]
	tlake$year = y
	all_lakes_years[[y]] = tlake
}

all_lakes_years_dt = do.call(rbind, all_lakes_years)

# Join CDEC data and lake years data 
snowpack = merge(all_lakes_years_dt, april_dat[, .(CDEC, SWEmm, year)], by=c("CDEC", "year"), all.x=TRUE)
colnames(snowpack)[3] = "lake_id"

# Save CDEC snowpack data
fwrite(snowpack, "../data/environmental_covariates/cdec_snowpack_april1.csv")

# Spot check
# snowpack[50000, ]
# april_dat[CDEC == "KIB" & year == 2009]

