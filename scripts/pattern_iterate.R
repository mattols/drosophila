#
# L/K at multiple levels
#
#

library(terra);library(sf)
library(spatstat)
library(dplyr)

# pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.drosophila.gbif.csv"
pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.locale.gbif.csv"
pth2 = "/Users/mattolson/src/drosophila-bioclim/data/output/university_loc_rankings.csv"

# GBIF DATA
dfg = read.csv(pth)
gv <- vect(dfg, geom=c("long", "lat"), crs="EPSG:4326")

# UNIVERSITY DATA
dfr = read.csv(pth2)
head(dfr);dim(dfr)
uv <- vect(dfr[dfr$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")
dim(uv)

# classify points in eco
gv$eco = NA
for(ue in unique(eco_realm2$WWF_REALM2)){
  # gv[eco_realm2[eco_realm2$WWF_REALM2=="Australasia"],]
  pts_nam = gv[eco_realm2[eco_realm2$WWF_REALM2==ue],]
  gv$eco[pts_nam$X] = ue
}
plot(gv, "eco")

# remove oceania?
sum(length(gv[gv$eco=="Oceania",1]))
# 198 pts???
# - should still remove due to sparseness

prep_k_data <- function(uv, gv, sp2){
  #
  # required
  #   uv - university vector object
  #   gv - gbif spatial object
  #   sp2 - window object
  # 
  
  # project and crop to extent
  univ = project(uv, sp2)[sp2,]
  gbif = project(gv, sp2)[sp2,]
  
  # create dataframes
  obs_coords <- as.data.frame(geom(gbif)[, c("x", "y")])
  school_coords <- as.data.frame(geom(univ)[, c("x", "y")])
  school_coords = school_coords %>% na.omit()
  
  # remove duplicates
  if(sum(duplicated(obs_coords))>100){
    print(paste("removing", sum(duplicated(obs_coords)), "duplicates..."))
    obs_coords = obs_coords[!duplicated(obs_coords[, c("x","y")]),]
  }
  
  # create ppp
  obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = my_window)
  school_ppp <- ppp(school_coords$x, school_coords$y, window = my_window)
  
  # rescale
  ppp_obs <- rescale(obs_ppp, 1000, "km")
  ppp_univ <- rescale(school_ppp, 1000, "km")
  
  # combine
  # Ensure both ppp objects share the same window
  common_window <- union.owin(ppp_obs$window, ppp_univ$window)
  # Set marks as you create the objects (not afterward)
  ppp_obs_marked <- ppp(x = ppp_obs$x, y = ppp_obs$y, window = common_window,
                        marks = factor(rep("obs", ppp_obs$n)))
  ppp_univ_marked <- ppp(x = ppp_univ$x, y = ppp_univ$y, window = common_window,
                         marks = factor(rep("univ", ppp_univ$n)))
  # Combine them using superimpose()
  ppp_all <- superimpose(ppp_obs_marked, ppp_univ_marked)
}





