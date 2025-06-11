
# spatial patterns
# https://cran.r-project.org/web/packages/spatialTIME/vignettes/intro.html
# http://www.integraess.com/GISPortfolio/Assign13/SpatialStatistics.html
# https://jo-wilkin.github.io/GEOG0030/coursebook/analysing-spatial-patterns-iii-point-pattern-analysis.html

# # # # #
library(terra);library(sf)
library(spatstat)
library(dplyr)

# GBIF records
pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.drosophila.gbif.csv"
pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.locale.gbif.csv"
dfg = read.csv(pth)
head(dfg);dim(dfg)
gv <- vect(dfg, geom=c("long", "lat"), crs="EPSG:4326")


# read university data
pth2 = "/Users/mattolson/src/drosophila-bioclim/data/output/university_loc_rankings.csv"
dfr = read.csv(pth2)
# dfr = dfr[dfr$Rank<=1000, ]
head(dfr);dim(dfr)
uv <- vect(dfr[dfr$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")
dim(uv)

# WINDOW
# results     <- "~/data/drosophila/results/data_1024"
# world       <- vect(file.path(results, "world-simple.geojson"))
# world       <- read_sf(file.path(results, "world-simple.geojson"))
# https://stackoverflow.com/questions/75927165/error-in-wk-handle-wk-wkbwkb-s2-geography-writeroriented-oriented-loop-0
# sf_use_s2(FALSE)
# area1 = st_area(world)
# world <- world[as.numeric(area1) > 1e10, 1]
# my_boundary <- aggregate(world)

# https://gis.stackexchange.com/questions/421257/plot-filtered-patch-sizes-in-r-terra
r = rast("~/data/drosophila/worldclim_05/wc2.1_10m/wc2.1_10m_bio_1.tif")
ry = patches(r)
# z = cellSize(ry,unit="ha") |> zonal(ry, sum)
# head(z);dim(z)
# remove greenland and antarctica
plot(ry == 1634 )
plot(ry == 1 )
ry[ry==1 | ry>=1634] = NA

rz = zonal(cellSize(ry,unit="km"), ry, sum, as.raster=TRUE)
rz
plot(rz)

# size of island to eliminate?
quantile(unique(values(rz)), na.rm=T)
hist(values(rz)) 
plot(rz<120000) # 1.2e05 keep iceland and SE ASIA
plot(rz<=1e06) # aggressive (eliminate Iceland, NZ, Indonesia, Madagascar)
plot(gv,add=T,cex=0.6)

# remove small patches
s = ifel(rz < 1.2e5, NA, 1)
plot(s)
sp1 = as.polygons(s)
sp2 = project(sp1, "EPSG:6933")
plot(sp2, col='orange3')
plot(project(uv, sp2), col='red', add=T, cex=0.5)
# or, if you want to keep the patch numbers
# s = ifel(rz < 100, NA, y)

my_window <- as.owin(st_as_sf(sp2))  


# create datafram
df_obs <- rbind(cbind(dfr[, c("long", "lat")], type = "university"),
                cbind(dfg[, c("long", "lat")], type = "observation") ) %>% 
  na.omit() %>% mutate(type = as.factor(type))
dfov <- vect(df_obs, geom=c("long", "lat"), crs = "EPSG:4326")
dfov2 <- project(dfov, sp2)
df_obs2 <- cbind(geom(dfov2)[, c("x", "y")], as.data.frame(dfov2)) %>% mutate(type = as.factor(type))
# plot(dfov, "type")
dim(df_obs)
univ_obs <- ppp(x = df_obs2$x, y = df_obs2$y, marks = df_obs2$type, window = my_window)
univ_obs <- ppp(x = df_obs$long, y = df_obs$lat, marks = df_obs$type, window = owin(c(-180, 180), c(-90, 90)) )
univ_obs <- ppp(x = df_obs2$x, y = df_obs2$y, marks = df_obs2$type, 
                window = owin(c(-17367530, 17367530), c(-7342230, 7297021)))

  # Warning messages:
  # 1: 1961 points were rejected as lying outside the specified window 
  # 2: data contain duplicated points 
plot(univ_obs)

# https://jo-wilkin.github.io/GEOG0030/coursebook/analysing-spatial-patterns-iii-point-pattern-analysis.html
# duplicated points 
anyDuplicated(univ_obs)
# Count the number of duplicated points and sum this
# sum(multiplicity(univ_obs) > 1) # does not finish

# Multiple kcross
uobs <- Kcross(univ_obs, i = "university", j = "observation")
plot(uobs, main = "", ylim=c(0, 18500))

uobs2 <- Kcross(univ_obs, i = "university", j = "observation")
plot(uobs2, main = "")

uobs3 <- Kcross(univ_obs, i = "university", j = "observation")
plot(uobs3, main = "")


stop()



# my_window <- as.owin(my_boundary)  
# birds <- readr::read_csv("https://drive.google.com/uc?export=download&id=1uaQu9LTLqnIjiIRQZLNlraRcLe6nsqkr")
#> Warning: Missing column names filled in: 'X1' [1]
# species <- factor(birds$species_name)
birds <- ppp(x = birds$X, y = birds$Y, marks = species, window = my_window)
#> Warning: data contain duplicated points
Kbirds <- Kcross(birds, i = "Species A", j = "Species B")
plot(Kbirds)





# # # # # # # # # # # # # # # #
# GPT - using Kest (univariate) - this actually ignores the second argument
# Install and load necessary packages
# install.packages("terra")
# install.packages("spatstat")

# Step 1: Load spatial data for observations1 and schools1
# Example: loading data using terra package
# observations1 <- vect("path_to_observations1.shp")  # Spatial points of observations
observations1 <- vect(dfg[,2:4], geom = c("long", "lat"), crs="EPSG:4326")
# schools1 <- vect("path_to_schools1.shp")  # Spatial points of schools
schools1 <- vect(dfr[,2:5], geom = c("long", "lat"), crs="EPSG:4326")

# Step 2: Convert the data into point patterns using spatstat package
# First, extract coordinates from the terra vectors (observations1, schools1)
# obs_coords <- as.data.frame(observations1)[, c("x", "y")]
obs_coords <- as.data.frame(geom(observations1)[, c("x", "y")])
# school_coords <- as.data.frame(schools1)[, c("x", "y")]
school_coords <- as.data.frame(geom(schools1)[, c("x", "y")])
school_coords = school_coords %>% na.omit()

# Create point pattern objects using spatstat
obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = owin(c(min(obs_coords$x), max(obs_coords$x)), 
                                                         c(min(obs_coords$y), max(obs_coords$y))))
school_ppp <- ppp(school_coords$x, school_coords$y, window = owin(c(min(school_coords$x), max(school_coords$x)), 
                                                                  c(min(school_coords$y), max(school_coords$y))))

# Step 3: Compute the Bivariate Ripley's K function
# This tests the spatial relationship between the two point patterns (observations1 and schools1)
cross_k <- cross.kppm(obs_ppp, school_ppp)
cross_k <- Kest(obs_ppp, school_ppp)

# Step 4: Plot the results
plot(cross_k, main = "Bivariate Ripley's K Function")

# Step 5: Interpretation of the results
# The plot shows the observed K-function and the expected K-function (under complete spatial randomness)
# If the observed K-function is above the expected K-function, it indicates positive spatial correlation (clustering)
# If the observed K-function is below the expected K-function, it indicates negative spatial correlation (dispersion)