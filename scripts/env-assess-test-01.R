
# Threshold - relevance drops off at a certain point

library(spatstat.core)

# Preprocess: get full set of candidate u points
all_u <- vect(dfr, geom = c("long", "lat"), crs = "EPSG:4326")

# Number of top-ranked anchor points to compare
anchor_n <- 500
anchor_n <- 1000
anchor_n <- 1500

# Subset top-ranked u points
uv_top <- vect(dfr[dfr$Rank <= anchor_n, ], 
               geom = c("long", "lat"), crs = "EPSG:4326")

# Observed point pattern (ppp object with "univ" and "obs" marks)
ppp_obs <- prep_k_data(uv_top, gv, sp2)

# Define simulation function: sample random u points, return pcfcross
simulate_fn <- function(i) {
  uv_rand <- all_u[sample(1:nrow(dfr), anchor_n), ]
  ppp_rand <- prep_k_data(uv_rand, gv, sp2)
  return(ppp_rand)
}

# Envelope computation
env <- envelope(
  Y = ppp_obs,
  fun = pcfcross.inhom,
  simulate = simulate_fn,
  i = "univ", j = "obs",
  nsim = 199,
  correction = "trans",
  savefuns = TRUE,
  verbose = TRUE
)

# Plot the envelope
plot(env, main = paste("Envelope: Top", anchor_n, "u vs Random"))



# # # # # # # # #
# subsets to assess?
dfr2 = dfr; dfr2$cnt = 1

z = zonal(vect(dfr2[dfr2$Rank<=1000, ], geom=c("long", "lat"), crs="EPSG:4326"), eco_simple2, fun='sum', as.polygons=T)
plot(z, "cnt")
text(z, "cnt", col='white')


# crazy MAUP in Palearctic

zonal(vect(dfr2[dfr2$Rank<=100, ], geom=c("long", "lat"), crs="EPSG:4326"), 
      eco_simple2, fun='sum', as.polygons=T) %>% 
  plot("cnt")

# 10 - only 1-Neotropic; 3-Nearctic, 6-Palearctic
# 

wrld <- vect("~/data/drosophila/world_shape/world-administrative-boundaries.geojson")
head(wrld)
plot(wrld, "name", col=rainbow(length(wrld$name)))

us = vect('/Users/mattolson/data/drosophila/world_shape/s_18mr25/s_18mr25.shp')
us = vect(ms_simplify(st_as_sf(us)))

# 4762 universities

zonal(vect(dfr2[dfr2$Rank<=1500, ], geom=c("long", "lat"), crs="EPSG:4326"), 
      wrld, fun='sum', as.polygons=T) %>% 
  as.data.frame() %>% 
  arrange(desc(cnt)) %>% 
  select(name, cnt) %>% 
  head()

zonal(vect(dfr2, geom=c("long", "lat"), crs="EPSG:4326"), 
      wrld, fun='sum', as.polygons=T) %>% 
  as.data.frame() %>% 
  arrange(desc(cnt)) %>% 
  select(name, cnt) %>% 
  head(12)

# smaller GBIF
gv2 = gv;gv2$cnt=1

# ALL GBIF
gv_all_csv = read.csv('/Users/mattolson/src/drosophila-bioclim/data/output/all.drosophila.gbif.csv')
dim(gv_all_csv)
gv_all <- vect(gv_all_csv, geom=c("long", "lat"), crs="EPSG:4326")

gv2 = gv_all;gv2$cnt=1


zonal(gv2,
      wrld, fun='sum', as.polygons=T) %>% 
  as.data.frame() %>% 
  arrange(desc(cnt)) %>% 
  select(name, cnt) %>% 
  head(12)


# normalization
wrld$area = terra::expanse(wrld, unit='km')

zonal(gv2, wrld, fun='sum', as.polygons=T) %>% 
  as.data.frame() %>% 
  mutate(cnt_area=cnt/(area/(10*10))) %>% #100sqkm 
  arrange(desc(cnt_area)) %>% 
  select(name, area, cnt_area, cnt) %>%
  mutate(area = round(area)) %>% 
  head(12)


zonal(gv2, wrld, fun='sum', as.polygons=T) %>% 
  as.data.frame() %>% 
  full_join(., as.data.frame(zonal(vect(dfr2, geom=c("long", "lat"), crs="EPSG:4326"),
                wrld, fun='sum', as.polygons=T) ), by='name') %>% 
  as.data.frame() %>% 
  mutate(cnt_spec_uv=cnt.x/cnt.y) %>% 
    arrange(desc(cnt_spec_uv)) %>%
  rename(cnt_obs = cnt.x, cnt_univ = cnt.y) %>% 
  select(name, cnt_spec_uv, cnt_obs, cnt_univ) %>%
  head(12)


# to assess better we need to subset for areas with most dens observations
# normalized


# to do:
# - split - create new windows for top 10 countries
# create subset - run for all 


