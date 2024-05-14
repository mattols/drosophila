#
# Area distribution based on temperature tolerance
# Drosophila Species - Hjelman & Curnow
#
library(terra);library(dplyr)

# 30s resolution
# download data (1/2 min or 30s resolution)
outpath = "~/data/drosophila/worldclim_05"
library(geodata)
bioclim19 <- worldclim_global('bio', res=10, path=outpath) #0.5 mins of degree
# ls_bio <- list.files(list.files(outpath, full.names = T), full.names = T)
# ls_bio <- list.files("~/data/drosophila/worldclim_05/wc2.1_10m/", full.names = T)

# 10m res
# outpath = "archive/data/wc2.1_10m/"
# ls_bio <- list.files(outpath, full.names = T)
# bioclim19 <- rast(ls_bio)

# info
nlyr(bioclim19) # dimensions of data (layers)
ext(bioclim19)
names(bioclim19)

# subset variables of interest
tmax <- subset(bioclim19, 5);names(tmax) <- "BIO5 - Max Temperature of Warmest Month"
tmin <- subset(bioclim19, 6);names(tmin) <- "BIO6 - Min Temperature of Coldest Month"
tmeanWq <- subset(bioclim19, 10);names(tmin) <- "BI10 - Mean Temperature of Warmest Quarter"
plot(c(tmin,tmax), col=rev(heat.colors(15)))

# scaling?
# https://worldclim.org/data/v1.4/formats.html


# read in data - omit NA values
df <- read.csv("archive/climate_data_drosophila_June.csv")
df2 <- df %>% dplyr::select(Species,Subgenus,Tmax,Tmin) %>% na.omit()
head(df2)

df2 <- df %>% dplyr::select(Species,Subgenus,Tmax,Tmin) %>% na.omit()


# tmin for first two species
tmin.cpy <- rast(replicate(nrow(df2),tmin))
tmax.cpy <- rast(replicate(nrow(df2),tmax))
nlyr(tmin.cpy)
names(tmin.cpy) <- df2$Species
# place threshold over lower bounds
tmin.dist <- rast(sapply(1:nrow(df2), 
                         function(x) clamp(subset(tmin.cpy, x), lower=df2[,'Tmin'][x], value=FALSE)))

# place threshold over upper bounds
tmax.dist <- rast(sapply(1:nrow(df2), 
                         function(x) clamp(subset(tmax.cpy, x), upper=df2[,'Tmax'][x], value=FALSE)))


# distribution of first four species in the list
plot(tmin.dist[[1:4]], col=rev(heat.colors(15)))
plot(tmax.dist[[1:4]], col=rev(heat.colors(15)))

# where both conditions are met
plot(!is.na(tmin.dist[[1]]) & !is.na(tmin))

# stats
xxcomb <- !is.na(tmin.dist) & !is.na(tmax.dist) & !is.na(rast(replicate(nlyr(tmin.dist),tmin)))

# ## OLD PIXEL COUNT
# #
# lyr_summary <- function(lyr){
#   rowinfo <- crds(lyr, df=TRUE) %>% filter(values(lyr)==1) %>% 
#     summarize(n=n(), latmean=mean(y), latmax=max(y), latmin=min(y) )
#   return(rowinfo)
# }
# 
# xd2 <- sapply(1:nlyr(xxcomb), function(x) lyr_summary(xxcomb[[x]])  ) 
# 
# df3 <- as.data.frame(t(xd2)) %>% mutate(Species = df2$Species, .before=n) %>% full_join(df2, by='Species') %>%
#   mutate(n = unlist(n), latmean = unlist(latmean),latmax = unlist(latmax),latmin = unlist(latmin))
# 
# # Top Area Code
# df3 %>% arrange(desc(n)) %>% head()
# top_area <- df3 %>% slice_max(n, n=4) %>% select(Species)
# 
# op <- par(cex = 0.5)
# tstk = xxcomb[[top_area[,1]]]
# plot(tstk*cz, col=c(NA,rev(colorspace::terrain_hcl(7))[2:6]), cex = 1,
#      plg = list(title=bquote("area"~km^2), cex=0.8), mar=c(0,1.5,0,3.5))
# 
# # plot all four
# library(tidyterra)
# library(ggtext)
# 
# test01 = tstk*cz
# test01[test01==0] = NA
# 
# ggplot() +
#   geom_spatraster(data = test01) +
#   facet_wrap(~lyr, ncol = 2) +
#   scale_fill_whitebox_c(
#     palette = "muted",
#     # labels = scales::label_number(suffix = " km^2"),
#     n.breaks = 10,
#     guide = guide_legend(reverse = TRUE)
#   ) +
#   labs(
#     fill = bquote("area"~km^2),
#     title = "Drosophila species with largest distribution",
#     subtitle = "italic(Based on Tmin + Tmax)"
#   ) +
#   theme_minimal()
#
#
# ggplot() +
#   geom_spatraster(data = test01) +
#   facet_wrap(~lyr, ncol = 2) +
#   scale_fill_whitebox_c(
#     palette = "arid",
#     # labels = scales::label_number(suffix = " km^2"),
#     n.breaks = 10,
#     guide = guide_legend(reverse = TRUE,
#                          title.position="left")
#   ) +
#   labs(
#     fill = bquote("area"~km^2),
#     title = "Drosophila species with largest distribution",
#     subtitle = "italic(Based on Tmin + Tmax)"
#   ) +
#   geom_spatvector(fill = NA) +
#   theme_light()

# map and projection tangent
world_map = map_data("world") %>% 
  filter(! long > 180)

countries = world_map %>% 
  distinct(region) %>% 
  tibble::rowid_to_column()

countries %>% 
  ggplot(aes(fill = rowid, map_id = region)) +
  geom_map(map = world_map) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("moll") +
  ggthemes::theme_map()
##

# Correct for area
cz <- cellSize(tmin[[1]], mask=TRUE, unit="km")
cz
plot(cz)

# PLOT
library(sf)
# Transform data to sf object
d <- st_as_sf(as.data.frame(cz, xy = TRUE), coords = c("x", "y"))
# Assign CRS
st_crs(d) <- 4326
# Plot
library(ggplot2)
ggplot(d) + geom_sf() +
  geom_raster(data = as.data.frame(cz, xy = TRUE),
              aes(x = x, y = y, fill = area))

# plot
plot(xxcomb[[c(1,3,5,7)]]*cz, col=c(NA,rev(colorspace::terrain_hcl(10))),
     plg = list(title="Area_km2", cex=0.75), mar=c(2,2,1,3))


## AREA CORRECTION
###
# calculate correct area
naras = (!is.na(cz));naras[naras==0]=NA
xdf1 = cbind(crds(xxcomb[[1]]*cz, df=TRUE), as.data.frame(xxcomb[[1]]*cz))
names(xdf1) = c("x","y", "area_km2")
xdf1 %>% 
  filter(na.omit(values(xxcomb[[1]]*naras))==1) %>% 
  summarize(area_mkm2=sum(area_km2)/1e6, latmean=mean(y), latmax=max(y), latmin=min(y) ) %>% 
  mutate(latrange = abs(latmin - latmax))

# LOOP summarize area
lyr_summary_area <- function(lyr){
  xdf1 = cbind(crds(lyr*cz, df=TRUE), as.data.frame(lyr*cz))
  names(xdf1) = c("x","y", "area_km2")
  xdf1 = xdf1 %>% 
    filter(na.omit(values(lyr*naras))==1) %>% 
    summarize(area_mkm2=sum(area_km2)/1e6, latmean=mean(y), latmax=max(y), latmin=min(y) ) %>% 
    mutate(latrange = abs(latmin - latmax))
  return(xdf1)
}

# new method
xd2 <- sapply(1:nlyr(xxcomb), function(x) lyr_summary_area(xxcomb[[x]])  ) 
df4 <- as.data.frame(t(xd2)) %>% mutate(Species = df2$Species, .before=area_mkm2) %>% full_join(df2, by='Species') %>%
  mutate(area_mkm2 = unlist(area_mkm2), latmean = unlist(latmean),
         latmax = unlist(latmax),latmin = unlist(latmin), latrange = unlist(latrange))
head(df4)

df5 = left_join(df4, df %>% select(Species, MbDNA_Male, MbDNA_Female), by='Species')
write.csv(df5,"drosophila_area.csv",row.names=FALSE)

df5 %>% ggplot(aes(x = area_mkm2, y = MbDNA_Male)) +
  geom_point()

# # # # # # #
## MAP PLOTS
# add world outline
# w <- world(path=tempdir())
w %>% ggplot() + geom_sf()

# all together
tot = app(xxcomb, fun="sum")
xarea = !is.na(cz);xarea[xarea==0]=NA;tot1 = tot*xarea
tot2 = tot;tot2[tot2==0] = NA
tot3 = project(tot2, "+proj=eck4 +lon_0=0 +datum=WGS84 +units=m +no_defs")
#
tot4 = app(xxcomb, fun="sum", na.rm=T)
xarea = !is.na(tmax[[1]]);xarea[xarea==0]=NA;tot1 = tot4*xarea


w %>% ggplot() + 
  # geom_spatvector(fill = alpha("grey75", 0.5) ) +
  #geom_sf(colour = alpha("khaki", 0.5) ) +
  geom_spatraster(data = tot1) +
  # facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(
    palette = "muted",
    # labels = scales::label_number(suffix = " km^2"),
    n.breaks = 10,
    na.value=alpha("lightblue",0.5),
    guide = guide_legend(reverse = TRUE) 
    # legend.location = "bottom",
    # label.position = "bottom")
  ) +
  labs(
    fill = "Count",
    title = "Drosophila Species Distribution (n = 61)",
    subtitle = "Based on temperature range"
  ) +
  # coord_map("moll") +
  geom_spatvector(fill = NA ) +
  # theme(legend.position = "bottom") +
  theme_minimal()



# test 2
w %>% ggplot() + 
  # geom_spatvector(fill = alpha("grey75", 0.5) ) +
  #geom_sf(colour = alpha("khaki", 0.5) ) +
  geom_spatraster(data = tot1) +
  # facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(
    palette = "viridi",
    # labels = scales::label_number(suffix = " km^2"),
    n.breaks = 10,
    # na.value=alpha("lightblue",0.5),
    guide = guide_legend(reverse = TRUE) 
    # legend.location = "bottom",
    # label.position = "bottom")
  ) +
  labs(
    fill = "Count",
    title = "Drosophila Species Distribution (n = 61)",
  ) +
  # coord_map("moll") +
  # geom_spatvector(fill = NA ) +
  # theme(legend.position = "bottom") +
  theme_minimal()

#
w %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  #geom_sf(colour = alpha("khaki", 0.5) ) +
  geom_spatraster(data = tot2) +
  # facet_wrap(~lyr, ncol = 2) +
  scale_fill_whitebox_c(
    palette = "muted",
    # labels = scales::label_number(suffix = " km^2"),
    n.breaks = 10,
    na.value=alpha("lightblue",0.3),
    guide = "colorbar"
    # guide = guide_legend(reverse = TRUE) 
    # legend.location = "bottom",
    # label.position = "bottom")
  ) +
  labs(
    fill = "",
    title = "Drosophila Species Distribution Counts (n = 61)"
    # subtitle = "Temperature tolerance and WorldCtim 2.1 temperature range"
  ) +
  # coord_map("moll") +
  # geom_spatvector(fill = alpha("grey75", 0.5) ) +
  # theme(legend.position = "bottom") +
  # theme(plot.subtitle = element_text(face="italic", size=6)) +
  theme_minimal() + theme(legend.position = "bottom", legend.direction="horizontal",
                          legend.key.width = unit(1, "in"), legend.title.position = "left")



w %>% ggplot() + 
  # geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_spatraster_contour_filled(data = tot1) +
  scale_fill_whitebox_d(
    palette = "viridi",
    # labels = scales::label_number(suffix = " km^2"),
    na.value=alpha("lightblue",0.3),
    guide = "colorbar"
    # guide = guide_legend(reverse = TRUE) 
    # legend.location = "bottom",
    # label.position = "bottom")
  ) +
  labs(
    fill = "",
    title = "Drosophila Distribution Contours (n = 61)"
    # subtitle = "Temperature tolerance and WorldCtim 2.1 temperature range"
  ) +
  theme_minimal() + theme(legend.position = "bottom", legend.direction="horizontal",
                          legend.key.width = unit(1, "in"), legend.title.position = "left")

## ????
# different theme
countries %>% 
  ggplot() +
  geom_map(map = world_map, aes(map_id = region), fill = "grey") +
  geom_spatraster(data = tot2) +
  scale_fill_whitebox_c(
    palette = "viridi",
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE) ) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  ggthemes::theme_map()



# # # # # # # # # #
##########################
# ZONAL stats for continents
wsh = vect("../../data/drosophila/world_shape/world-administrative-boundaries.geojson")
head(wsh)
plot(wsh, "continent")
wshd = aggregate(wsh, "continent")
head(wshd)
dim(wsh);dim(wshd)
plot(wshd, "continent")

# zonal stats
i = 1
df_cont = df2 %>% select(Species)
df_cont[wshd$continent] = NA
for (i in 1:nrow(df_cont)){
  print(paste(i,"of", nrow(df_cont)))
  df_cont[i,2:7] = t(zonal(xxcomb[[df_cont$Species[i]]]*cz, wshd, fun = "sum", na.rm=T))
}

# distribution 
library(tidyr)
df_cont %>% pivot_longer(cols = !Species,
                         names_to = "Continent",
                         values_to = "values") %>% 
  # mutate(Coninent = as.factor(Continent)) %>% 
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = Continent, y = values, fill = Species)) +
  labs(
    x = "", y = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila Continental Distribution"
  ) +
  geom_bar(stat="identity") + 
  scale_fill_discrete(guide="none") +
  # guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  # guides(color = guide_legend(override.aes = list(size = 0.5)))+
  # theme(legend.title = element_text(size = 3), 
  #       legend.text = element_text(size = 3))+
  # guides(fill = guide_legend(show=F)) +
  theme_classic()

df_cont %>% pivot_longer(cols = !Species,
                         names_to = "Continent",
                         values_to = "values") %>% 
  # mutate(Coninent = as.factor(Continent)) %>% 
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = Continent, y = values)) +
  labs(
    x = "", y = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila Continental Distribution"
  ) +
  geom_bar(stat="identity") + 
  scale_fill_discrete(guide="none") +
  # guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  # guides(color = guide_legend(override.aes = list(size = 0.5)))+
  # theme(legend.title = element_text(size = 3), 
  #       legend.text = element_text(size = 3))+
  # guides(fill = guide_legend(show=F)) +
  theme_classic()


df5 %>% ggplot(aes(x = area_mkm2, y = MbDNA_Male)) +
  geom_point()

df67 = left_join(df_cont, df %>% select(Species, Subgenus,MbDNA_Male, MbDNA_Female), by='Species')

df67 %>% pivot_longer(cols = !c(Species,MbDNA_Male, MbDNA_Female),
                      names_to = "Continent",
                      values_to = "values") %>%  
  ggplot(aes(x = values, y = MbDNA_Male, colour=Continent)) +
  geom_point() +facet_wrap(.~Continent)

left_join(df_cont, df %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>%
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Continent",
               values_to = "values") %>%  
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = values, y = MbDNA_Male, colour=Subgenus)) +
  labs(title = "Drosophila Continental Distribution and Genome Size")
# guides(fill = guide_legend(show=F)) +
geom_point() + facet_wrap(.~Continent) +
  theme_minimal()





# # # # # # # # # #
##########################
#### NATURE CONSERVANCY ecoregions
# https://geospatial.tnc.org/datasets/b1636d640ede4d6ca8f5e369f2dc368b/about
v <- vect("./../../data/drosophila/world_shape/terr-ecoregions-TNC/tnc_terr_ecoregions.shp")
head(v)
plot(v, "WWF_REALM")
glimpse(v)

vrealm <- aggregate(v, "WWF_REALM")
plot(vrealm, "WWF_REALM")

vrealm2 <- aggregate(v, "WWF_REALM2")
plot(vrealm2, "WWF_REALM2")

# omit islands
dim(v)
crs(v)
exp_v <- expanse(v, unit="km")
v$area <- exp_v
as.data.frame(v) %>% arrange(area) %>% select(area) %>% plot()
abline(v=1.4e6)
hist(v$area)

dim(v[v$area > 1.4e6,])
vrealm2n <- aggregate(v[v$area > 1.4e6,], "WWF_REALM2")
plot(vrealm2n, "WWF_REALM2")

# check out
# ZONAL STATS by ecoregion

# zonal stats
i = 1
df_eco = df5 %>% select(Species)
df_eco[vrealm2$WWF_REALM2] = NA
for (i in 1:nrow(df_eco)){
  print(paste(i,"of", nrow(df_eco)))
  df_eco[i,2:ncol(df_eco)] = t(zonal(xxcomb[[df_eco$Species[i]]]*cz, vrealm2, fun = "sum", na.rm=T))
}

# distribution 
library(tidyr)
library(ggplot2)
df_eco %>% pivot_longer(cols = !Species,
                        names_to = "Biogeographic",
                        values_to = "values") %>% 
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = Biogeographic, y = values, fill = Species)) +
  labs(
    x = "", y = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila Biogeographic Distributions"
  ) +
  geom_bar(stat="identity") + 
  scale_fill_discrete(guide="none") +
  theme_classic()

### MAP 1
m1 = left_join(df_eco, df5 %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  filter(Biogeographic != "Antarctic", Biogeographic != "Oceania" ) %>% 
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = values, y = MbDNA_Male, colour=Biogeographic)) +
  labs(title = "Drosophila Distribution and Genome Size",
       subtitle = "by Biogeographic Region",
       x = bquote("Area by temperature tolerance (million"~km^2~")")) +
  scale_colour_brewer(palette = "Dark2", guide = 'none') +
  # guides(color=guide_legend(title="Bio Region")) +
  # scale_color_discrete(guide = 'none') +  
  geom_point() + facet_wrap(.~Biogeographic) +
  theme(legend.text=element_text(size=rel(0.5))) +
  theme_minimal()

# MAP 2
# vrsf = st_as_sf(vrealm2) 
vrsf$species_tot_area_mkm = colSums(df_eco[,2:ncol(df_eco)])/1e6
dmax = t(zonal(tot2, vrealm2, fun = "max", na.rm=T))
vrsf$species_max = as.numeric(dmax)
m2 = vrsf %>%
  filter(WWF_REALM2 != "Antarctic", WWF_REALM2 != "Oceania" ) %>% 
  st_transform(6933) %>% 
  ggplot() + 
  labs(#title = "Drosophila Species Count by Biogeographic Region (n=61)",
    x = "",y="") +
  geom_sf(aes(fill = WWF_REALM2)) + 
  guides(fill=guide_legend(title="Species count (n=61)")) +
  scale_fill_brewer(palette = "Dark2") +
  geom_sf_label(aes(label = species_max, alpha = 0.7))  +
  scale_alpha(guide = 'none') +
  theme_minimal() 
#+ coord_sf("moll")
# coord_equal()


## PLOT BOTH
library(gridExtra)
# png("figures/drosophila_biogegraphy.png", width = 9, height = 4, units = "in", res = 300)
grid.arrange(m1,m2)
# dev.off()


vrsf %>% 
  ggplot(aes(fill = rowid, map_id = WWF_REALM2)) +
  geom_map(map = vrealm2) +
  # geom_spatraster(data = tot2) +
  # expand_limits(x = vrealm2$long, y = world_map$lat) +
  # coord_sf("moll") +
  theme_map()

# left_join(df_eco, df5 %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>%
#   pivot_longer(cols = !c(Species, Subgenus, MbDNA_Male, MbDNA_Female),
#                names_to = "Biogeographic",
#                values_to = "values") %>%  
#   mutate(values = values/1e6) %>% 
#   ggplot(aes(x = values, y = MbDNA_Male, colour=Subgenus)) +
#   labs(title = "Drosophila Continental Distribution and Genome Size") +
#   # guides(fill = guide_legend(show=F)) +
#   geom_point() + facet_wrap(.~Biogeographic) +
#   theme_minimal()


left_join(df_eco, df5 %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>%
  pivot_longer(cols = !c(Species, Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>%  
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = values, y = MbDNA_Male, colour=Biogeographic)) +
  labs(title = "Drosophila Biogeographic Distribution and Genome Size",
       x = bquote("Area by temperature min/max threshold (million"~km^2~")"),
       legend = "Region") +
  scale_colour_brewer(palette = "Dark2") +
  # guides(fill = guide_legend(show=F)) +
  geom_point() + 
  # facet_wrap(.~Biogeographic) +
  theme_minimal()







# SENT TO CARL
# simple plot
library(terra)
dros_heat = rast("/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_1.tif")
maps::map("world", col="grey40")
plot(dros_heat, col=heat.colors(12), main = "Drosophila temperature tolerance")
maps::map("world", add=T)

# load world political data
library(sf);library(ggplot2);library(dplyr)
world = read_sf("../../data/drosophila/world_shape/world-simple.geojson")
# ggplot w/ tidyterra
library(tidyterra)
world %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_spatraster(data = dros_heat) +
  scale_fill_whitebox_c(
    palette = "muted", # can use "viridi"
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE)) +
  labs(
    fill = "Counts",
    title = "Drosophila Distribution (n = 61)",
    subtitle = "based on temperature tolerance") +
  ylim(c(-60,85)) + theme_minimal()


# or
m = maps::map("world", plot = F)
m %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_spatraster(data = dros_heat) +
  scale_fill_whitebox_c(
    palette = "muted", # can use "viridi"
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE)) +
  labs(
    fill = "Counts",
    title = "Drosophila Distribution (n = 61)",
    subtitle = "based on temperature tolerance") +
  theme_minimal()


world %>% ggplot() + 
  geom_spatvector(fill = alpha("grey75", 0.5) ) +
  geom_spatraster(data = dros_heat) +
  scale_fill_whitebox_c(
    palette = "muted", # can use "viridi
    n.breaks = 12,
    na.value=alpha("lightblue",0.2),
    guide = guide_legend(reverse = TRUE)) +
  ylim(c(-60,85)) +
  labs(
    fill = "Count",
    title = "Drosophila Distribution (n = 61)",
    subtitle = "Based on temperature tolerance"
  ) +
  theme_minimal()

#
# horizontal bar on bottom
# theme(
#     legend.position = "bottom", 
#     legend.direction="horizontal",
#     legend.key.width = unit(1, "in"), 
#     legend.title.position = "top")

#
###### SIMPLIFY
world_o = read_sf("../../data/drosophila/world_shape/world-administrative-boundaries.geojson")
# st_write(w, dsn = "../../data/drosophila/world_shape/world-simple.geojson")
# layer = "world-simple.geojson")
w = world %>% 
  st_transform(6933) %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 1e3 ) %>% 
  select(color_code, name, continent, region) %>% 
  st_transform(4326)

plot(w)



#### Write/READ files

# writeRaster(xxcomb, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/dros_Tdist_n61.grd")
# writeRaster(cz, "/Users/mattolson/data/drosophila/results/rast0/global_area_10m.grd")
# writeRaster(tot1, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_0.grd")
# writeRaster(tot2, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_1.grd")

xxcomb = rast("/Users/mattolson/data/drosophila/results/Tdistribution_n61/dros_Tdist_n61.grd")
cz = rast("/Users/mattolson/data/drosophila/results/rast0/global_area_10m.grd")
tot1 = rast(, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_0.grd")
tot2 = rast(, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_1.grd")

# write.csv(df5, "/Users/mattolson/data/drosophila/results/tables/drosophila_area_table.csv", row.names = T)
# write.csv(df67, "/Users/mattolson/data/drosophila/results/tables/drosophila_continent_areakm_table.csv", row.names = T)

df5 = read.csv("/Users/mattolson/data/drosophila/results/tables/drosophila_area_table.csv")
df67 = read.csv("/Users/mattolson/data/drosophila/results/tables/drosophila_continent_areakm_table.csv")

dros_heat = rast(, "/Users/mattolson/data/drosophila/results/Tdistribution_n61/drosophila_heatmap_n61_10m_1.grd")


