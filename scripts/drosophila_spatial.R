#
# Area distribution based on species temperature tolerance from WorldClim v2.1
# Drosophila Project - Hjelman & Curnow
# Code - Olson 10/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # #

library(terra);library(dplyr);library(geodata)
library(sf);library(tidyr)

# WorldClim v2.1 data
outpath = "./worldclim"
bioclim19 <- worldclim_global('bio', res=10, path=outpath) #10 mins of a degree

# subset variables of interest
tmax <- subset(bioclim19, 5);names(tmax) <- "BIO5 - Max Temperature of Warmest Month"
tmin <- subset(bioclim19, 6);names(tmin) <- "BIO6 - Min Temperature of Coldest Month"
tmeanWq <- subset(bioclim19, 10);names(tmin) <- "BI10 - Mean Temperature of Warmest Quarter"
plot(c(tmin,tmax), col=rev(heat.colors(15)))

# read in data - omit NA values
df <- read.csv("./climate_data_drosophila_Sept17.csv")
df2 <- df %>% dplyr::select(Species,Subgenus,TMax,TMin) %>% na.omit()

# replicate grid for all species
tmin.cpy <- rast(replicate(nrow(df2),tmin))
tmax.cpy <- rast(replicate(nrow(df2),tmax))
names(tmin.cpy) <- df2$Species

# place threshold over lower bounds
tmin.dist <- rast(sapply(1:nrow(df2), 
                         function(x) clamp(subset(tmin.cpy, x), lower=df2[,'TMin'][x], value=FALSE)))

# place threshold over upper bounds
tmax.dist <- rast(sapply(1:nrow(df2), 
                         function(x) clamp(subset(tmax.cpy, x), upper=df2[,'TMax'][x], value=FALSE)))

# conditions for all species
dros_dist <- !is.na(tmin.dist) & !is.na(tmax.dist) & !is.na(rast(replicate(nlyr(tmin.dist),tmin)))

# global area correction layers
cz <- cellSize(tmin[[1]], mask=TRUE, unit="km")
naras = (!is.na(cz));naras[naras==0]=NA 

# total counts and heatmap
dros_heat = app(dros_dist, fun="sum")
dros_heat[dros_heat==0] = NA

# function to summarize by area
lyr_summary_area <- function(lyr){
  xdf_summary = cbind(crds(lyr*cz, df=TRUE), as.data.frame(lyr*cz))
  names(xdf_summary) = c("x","y", "area_km2")
  xdf_summary = xdf_summary %>% 
    filter(na.omit(values(lyr*naras))==1) %>% 
    summarize(area_mkm2=sum(area_km2)/1e6, latmean=mean(y), latmax=max(y), latmin=min(y) ) %>% 
    mutate(latrange = abs(latmin - latmax))
  return(xdf_summary)
}

# summarize area and location stats by species
xd_summarize_01 <- sapply(1:nlyr(dros_dist), function(x) lyr_summary_area(dros_dist[[x]])  ) 
df_area <- as.data.frame(t(xd_summarize_01)) %>% mutate(Species = df2$Species, .before=area_mkm2) %>% full_join(df2, by='Species') %>%
  mutate(area_mkm2 = unlist(area_mkm2), latmean = unlist(latmean),
         latmax = unlist(latmax),latmin = unlist(latmin), latrange = unlist(latrange))
head(df_area)

# combine with original
df_comb = left_join(df_area, df %>% select(Species, MbDNA_Male, MbDNA_Female), by='Species')

# aggregate layer
wsh = vect("./data_1024/world-simple.geojson")
wsh_cont = aggregate(wsh, "continent")

# zonal statistics by continent
df_cont = df2 %>% select(Species)
df_cont[wsh_cont$continent] = NA
for (i in 1:nrow(df_cont)){
  print(paste(i,"of", nrow(df_cont)))
  df_cont[i,2:7] = t(zonal(dros_dist[[df_cont$Species[i]]]*cz, wsh_cont, fun = "sum", na.rm=T))
}

# nature conservancy ecoregions
eco_realm2 <- vect("./data_1024/eco-realm2.geojson")

# zonal stats by ecoregion
df_eco = df2 %>% select(Species)
df_eco[eco_realm2$WWF_REALM2] = NA
for (i in 1:nrow(df_eco)){
  print(paste(i,"of", nrow(df_eco)))
  df_eco[i,2:ncol(df_eco)] = t(zonal(dros_dist[[df_eco$Species[i]]]*cz, eco_realm2, fun = "sum", na.rm=T))
}

# full area of eco-zones (for normalization)
eco_tot = zonal(cz, eco_realm2, fun = "sum", na.rm=T)
df_eco_tot = data.frame(realm = eco_realm2$WWF_REALM2, total_area = eco_tot)

# distribution 
df_eco %>% pivot_longer(cols = !Species,
                        names_to = "Biogeographic",
                        values_to = "values") %>% 
  mutate(values = values/1e6) %>% 
  filter(Biogeographic != "Antarctic", Biogeographic != "Oceania" ) %>% 
  head()

### FIGURES
# drosophila_figures.R
