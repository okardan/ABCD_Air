# Map of US air pollution
# "Neighborhood air pollution is negatively associated with neurocognitive maturation in early adolescence"
# Omid Kardan, Chacriya Sereeyothin, Kathryn E. Schertz, Mike Angstadt, Alexander S. Weigard, Marc G. Berman, Monica D. Rosenberg

# script by Kathryn E Schertz 
# contact omidk@med.umich.edu
# Produces US county-level air pollution maps in Fig 1 and Supp Fig S1 (only for visualization purpose, not used in analysis)
# Requires air pollutant data for US counties from CDC and EPS

###########
library(tidyverse)
library(sf)
library(data.table)
library(RColorBrewer)

# Download county shapefile from https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
county_map <- st_read("cb_2018_us_county_20m/cb_2018_us_county_20m.shp")
county_map_mainland <- county_map %>% filter(!STATEFP %in%c("02","15","72"))
county_map_mainland[(county_map_mainland$GEOID=="46102"),"GEOID"] <- "46113"

# Download Daily PM2.5 Concentrations All County, 2001-2016 from https://data.cdc.gov/browse/select_dataset?tags=pm2.5
pm25_all <- fread("Daily_PM2.5_Concentrations_All_County__2001-2016.csv")
pm25_2016 <- pm25_all %>% filter(year == 2016)
pm25_2016_avg <- pm25_2016 %>% group_by(statefips, countyfips) %>% summarise(PM25_mean = mean(PM25_mean_pred, na.rm=TRUE))

pm25_2016_avg <- pm25_2016_avg %>% mutate(STATEFP = str_pad(as.character(statefips), 2, side = "left", pad="0"), COUNTYFP = str_pad(as.character(countyfips), 3, side = "left", pad = "0"))
pm25_2016_avg <- pm25_2016_avg %>% mutate(GEOID = paste0(STATEFP, COUNTYFP))
pm25_2016_avg <- pm25_2016_avg %>% ungroup()

county_map_mainland_pm25 <- county_map_mainland %>% left_join(pm25_2016_avg, by = "GEOID")

mid = mean(county_map_mainland_pm25$PM25_mean, na.rm=TRUE)

#With legend
ggplot(county_map_mainland_pm25) + aes(fill=PM25_mean) +
  geom_sf(size = .1, color = "white") + 
  scale_fill_gradient2(low = "seagreen3", high = "firebrick4", mid = "yellow", midpoint = mid, name = "Average PM 2.5", guide = "colourbar") + 
  theme_void() +
  theme(legend.direction = 'horizontal', legend.position = c(.25,.15)) + 
  guides(colourbar = guide_legend(title.position = "top", title.hjust = 0.5))

#Without legend
ggplot(county_map_mainland_pm25) +
  geom_sf(aes(fill=PM25_mean), size = .1, color = "white") + 
  scale_fill_gradient2(low = "seagreen3", high = "firebrick4", mid = "yellow", midpoint = mid) + 
  geom_sf(data = abcd_sites_sf) + 
  theme_void() +
  theme(legend.title = element_blank())

## NO2
# download NO2 data from from https://aqs.epa.gov/aqsweb/airdata/download_files.html#Annual

no2_2016 <- read_csv("annual_conc_by_monitor_2016.csv")
no2_2017 <- read_csv("annual_conc_by_monitor_2017.csv")
no2_2018 <- read_csv("annual_conc_by_monitor_2018.csv")
no2_2021 <- read_csv("annual_conc_by_monitor_2021.csv")

no2_allyears <- rbind(no2_2016, no2_2017, no2_2018, no2_2021)

no2_allyears2 <- no2_allyears %>% filter(str_detect(`Parameter Name`,"NO2"))
no2_hourlystandard <- no2_allyears2 %>% filter(`Pollutant Standard` == "NO2 1-hour 2010")
no2_hourlystandard <- no2_hourlystandard %>% mutate(GEOID = paste0(`State Code`, `County Code`))
no2_bycounty <- no2_hourlystandard %>% group_by(GEOID) %>% summarise(HourlyStandardPPB = mean(`Arithmetic Mean`, na.rm=TRUE))
no2_bycounty <- no2_bycounty %>% mutate(NO2_Z = scale(HourlyStandardPPB))

county_map_mainland_no2 <- county_map_mainland %>% left_join(no2_bycounty)
mid2= mean(no2_bycounty$HourlyStandardPPB, na.rm=TRUE)

#With legend
ggplot(county_map_mainland_no2) +
  geom_sf(aes(fill=HourlyStandardPPB), size = .1, color = "white") + 
  scale_fill_gradient2(low = "seagreen3", high = "firebrick4", mid = "yellow", midpoint = mid2) + 
  geom_sf(data = abcd_sites_sf) + 
  theme_void() +
  theme(legend.direction = 'horizontal', legend.position = c(.25,.15), legend.title = element_blank())

#Without legend
ggplot(county_map_mainland_no2) +
  geom_sf(aes(fill=HourlyStandardPPB), size = .1, color = "white") + 
  scale_fill_gradient2(low = "seagreen3", high = "firebrick4", mid = "yellow", midpoint = mid2) + 
  geom_sf(data = abcd_sites_sf, size = .2, color = "dodgerblue1") + 
  theme_void() +
  theme(legend.title = element_blank())


### Together
pm25_2016_avg <- pm25_2016_avg %>% mutate(PM25_Z = scale(PM25_mean))
no2_bycounty <- no2_bycounty %>% mutate(NO2_Z = scale(HourlyStandardPPB))
county_map_mainland_combo <- county_map_mainland %>% left_join(pm25_2016_avg, by="GEOID") %>% left_join(no2_bycounty)
county_map_mainland_combo <- county_map_mainland_combo %>% mutate(AQ = rowMeans(cbind(PM25_Z, NO2_Z), na.rm=TRUE))

mid3 = mean(county_map_mainland_combo$AQ, na.rm=TRUE)

ggplot(county_map_mainland_combo) +
  geom_sf(aes(fill=AQ), size = .1, color = "white") + 
  scale_fill_gradient2(low = "seagreen3", high = "firebrick4", mid = "yellow", midpoint = mid3) + 
  geom_sf(data = abcd_sites_sf, size = .2, color = "dodgerblue1") + 
  theme_void() +
  theme(legend.title = element_blank())
