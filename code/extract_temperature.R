library(data.table)
library(raster)
library(sp)
library(rgdal)
library(lubridate)
library(ncdf4)

## This script extracts temperature covariates for all Sierra lakes for 20 years from 2000-2020
## Month minimum and maximum summer temperature are extracted. Summer is defined as
## months c(6, 7, 8, 9)
##
## Temperature data is obtained from http://www.climatologylab.org/gridmet.html

### Preliminaries ###

# Load data
lake_dat = fread("../data/archival/sites_nov2020.csv")

# Set projection strings
crs_latlon = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
epsg_11N = "EPSG:32611"

# Transform to latlong
tdat = SpatialPoints(lake_dat[, .(utme, utmn)], proj4string=CRS(epsg_11N))
tdat_trans = spTransform(tdat, CRS(crs_latlon))
plot(tdat_trans)

### Check outliers. One point in way too far south it seems.
shp = shapefile("../data/environmental_covariates/ca-state-boundary/CA_State_TIGER2016.shp")
shp = spTransform(shp, CRS(crs_latlon))
plot(tdat_trans, pch=19)
lines(shp)
lake_dat[utmn == min(lake_dat$utmn)]
lake_dat[utme == max(lake_dat$utme)]
lake_dat[county == "fresno" & topo == "slide_bluffs"]

#### Extract temperature data ####

covpath = "/Volumes/SeagateBackup/serology_covariates"
folders = c("minimum_temperature", "maximum_temperature")
prefixes = c("tmmn", "tmmx")
years = 2000:2020
summer_months = c(6, 7, 8, 9)
keys = c("lake_id", "year")

year_temps = list()
for(i in 1:length(years)){

    yr = years[i]
    cat("Working on", yr, "\n")

    minmax_temps = list()
    for(j in 1:length(folders)){

        fl = folders[j]
        pre = prefixes[j]
        cat("Loading", fl, "for", yr, "\n")

        # Load temperature data
        tras = brick(file.path(covpath,
                              fl,
                              paste0(pre, "_", yr, ".nc")))
        cat("Loaded", "\n")

        crs(tras) = crs_latlon
        rasnames = names(tras)
        dayssince = unlist(lapply(rasnames, function(x) as.numeric(strsplit(x, "X")[[1]][2])))
        basetime = as.POSIXct("1900-01-01 00:00:00", tz="GMT")
        yeardates = basetime + as.difftime(dayssince, units="days")


        ttemp = extract(tras, tdat_trans) - 273.15 # Convert K to C
        colnames(ttemp) = yeardates
        tlake_temp = cbind(lake_dat$id, ttemp)
        colnames(tlake_temp)[1] = "lake_id"

        melt_lake = melt(data.table(as.data.frame(tlake_temp)), id.vars=c("lake_id"),
                         value.name = fl,
                         variable.name = "unix_time")

        # Convert to datetime
        melt_lake$unix_time = as.numeric(levels(melt_lake$unix_time))[melt_lake$unix_time]
        melt_lake$unix_time = as.POSIXct(melt_lake$unix_time, origin="1970-01-01", tz="GMT")

        # Group by summer months in year
        melt_lake$month = month(melt_lake$unix_time)
        melt_lake$year = year(melt_lake$unix_time)
        colnames(melt_lake)[colnames(melt_lake) == fl] = "temp_val"
        melt_lake_month = melt_lake[month %in% summer_months][, .(temp_val=mean(temp_val)),
                                                         by=.(lake_id, year)]
        colnames(melt_lake_month)[colnames(melt_lake_month) == "temp_val"] = fl
        minmax_temps[[j]] = melt_lake_month

    }

    full_minmax = merge(minmax_temps[[1]], minmax_temps[[2]], by=keys)
    year_temps[[i]] = full_minmax

}

# Merge and save temperatures
full_temps = do.call(rbind, year_temps)
colnames(full_temps) = c("lake_id", "year", "min_summer_temp_C", "max_summer_temp_C")
fwrite(full_temps, "../data/environmental_covariates/summer_lake_temperatures.csv")


































