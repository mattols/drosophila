#
# Point pattern analysis
# https://mgimond.github.io/Spatial/chp11_0.html
# 06/2025
# 

# IMPORTANT
# projection:
 # - important to project data
 # - K, L, PCF measure distances Euclidean distance
 # - use Mollweide or Eckert
# duplicates:
 # - Duplicate points can cause spurious clustering at small scales
 # - close but not exact points are usually fine if they are real observations
 # - too many zero-distance pairs can inflate statistics

# Methods: (bivariate)
#   K: Test for 
#   L:
#   PCF:


# Carl Feedback:
#   Test: gbif vs universities (top 100, 500, 1000, 1500, all)
#   Assess gbif with population density
#       points for pop over 10,000, 100,000 or continuous
#   Window: include a window
#
# redo figures and stats

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



r = rast("~/data/drosophila/worldclim_05/wc2.1_10m/wc2.1_10m_bio_1.tif")
ry = patches(r)
# z = cellSize(ry,unit="ha") |> zonal(ry, sum)
# head(z);dim(z)
# remove greenland and antarctica
ry[ry==1 | ry>=1634] = NA

rz = zonal(cellSize(ry,unit="km"), ry, sum, as.raster=TRUE)
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


univ = project(uv, sp2)[sp2,]
plot(sp2, col='orange3')
plot(univ, add=T, cex=0.5)
gbif = project(gv, sp2)[sp2,]
plot(gbif, add=T, pch=3, col='red')

obs_coords <- as.data.frame(geom(gbif)[, c("x", "y")])
school_coords <- as.data.frame(geom(univ)[, c("x", "y")])
school_coords = school_coords %>% na.omit()

### issues - data contain duplicate points - for both
# jitter or use marks or weights?

obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = my_window)
school_ppp <- ppp(school_coords$x, school_coords$y, window = my_window)

# duplicates
dim(obs_coords)
sum(duplicated(obs_coords))
sum(duplicated(school_coords))

obs_coords = obs_coords[!duplicated(obs_coords[, c("x","y")]),]
obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = my_window)

# obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = owin(c(min(obs_coords$x), max(obs_coords$x)), 
#                                                          c(min(obs_coords$y), max(obs_coords$y))))
# school_ppp <- ppp(school_coords$x, school_coords$y, window = owin(c(min(school_coords$x), max(school_coords$x)), 
#                                                                   c(min(school_coords$y), max(school_coords$y))))

# add jitter - DID NOT WORK (only 4 and 7 repeats so relatively small)
# o_coords <-  coords(obs_ppp)
# o_coords_jitter <- o_coords + runif(length(o_coords), -0.000005, 0.000005)
# o_ppp_jitt <- ppp(o_coords_jitter[,1], o_coords_jitter[,2], window = obs_ppp$window)
# none
# s_coords <-  coords(school_ppp)
# s_coords_jitter <- s_coords + runif(length(s_coords), -0.00005, 0.00005)
# s_ppp_jitt <- ppp(s_coords_jitter[,1], s_coords_jitter[,2], window = school_ppp$window)

sum(duplicated(o_coords))
sum(duplicated(o_coords_jitter))



gbif.km <- rescale(obs_ppp, 1000, "km")
univ.km <- rescale(school_ppp, 1000, "km")




# ppp_obs = gbif.km
# ppp_univ = univ.km
# # Add marks
# ppp_obs$marks <- factor("obs")
# ppp_univ$marks <- factor("univ")
# # Combine into a single multitype pattern
# ppp_all <- superimpose(ppp_obs, ppp_univ, W = union.owin(ppp_obs$window, ppp_univ$window))


# Autocorrelation - clustering
# clustering of single - second argument is ignored!
cross_k <- Kest(gbif.km, univ.km)
cross_L <- Lest(gbif.km, univ.km)
# cross_pcf <- pcf(gbif.km, univ.km) # DID NOT WORK

par(mar=c(5,5,1,6)) # 5.1 4.1 4.1 2.1
plot(cross_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))

par(mar=c(5,5,1,8)) # 5.1 4.1 4.1 2.1
plot(cross_L, . -r ~ r, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))


### CROSS CORRELATION
ppp_obs = gbif.km
ppp_univ = univ.km

# Ensure both ppp objects share the same window
common_window <- union.owin(ppp_obs$window, ppp_univ$window)

# Set marks as you create the objects (not afterward)
ppp_obs_marked <- ppp(x = ppp_obs$x, y = ppp_obs$y, window = common_window,
                      marks = factor(rep("obs", ppp_obs$n)))

ppp_univ_marked <- ppp(x = ppp_univ$x, y = ppp_univ$y, window = common_window,
                       marks = factor(rep("univ", ppp_univ$n)))

# Combine them using superimpose()
ppp_all <- superimpose(ppp_obs_marked, ppp_univ_marked)

# Now you can use Kcross
K_obs_univ <- Kcross(ppp_all, i = "univ", j = "obs")
plot(K_obs_univ)

# L_obs_univ <- Lcross(ppp_all, i = "obs", j = "univ")
# univ is i - obs (j) compute expected number of points around universities
L_obs_univ <- Lcross(ppp_all, i = "univ", j = "obs")

plot(L_obs_univ, . -r ~ r, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))



# K_obs_univ <- Kcross(ppp_all, i = "obs", j = "univ")
# plot(K_obs_univ)
env_cross <- envelope(ppp_all, Lcross, i = "univ", j = "obs", nsim = 99)
plot(env_cross,. -r ~ r)




# PCF
pcf_univ_obs <- pcfcross(ppp_all, i = "univ", j = "obs")

svpth = "~/src/drosophila/figs_0625"
png(file.path(svpth, "pcf_global_univ_obs_2.png"), width = 8, height=6, unit="in", res=300)
plot(pcf_univ_obs, main=NULL, lwd=2, axes=F, col=c("black", "firebrick", "grey"))
axis(1, cex=1.3)
axis(2, cex=1.3)
box()
dev.off()

head(as.data.frame(pcf_univ_obs), 20)
# r  6.471837    - 108.00851


pcf_univ_obs_inhom <- pcfcross.inhom(ppp_all, i = "univ", j = "obs")

svpth = "~/src/drosophila/figs_0625"
png(file.path(svpth, "pcf_inhom_global_univ_obs_2.png"), width = 8, height=6, unit="in", res=300)
plot(pcf_univ_obs_inhom, main=NULL, lwd=2, axes=F, col=c("black", "firebrick", "grey"))
axis(1, cex=1.3)
axis(2, cex=1.3)
box()
dev.off()
head(as.data.frame(pcf_univ_obs_inhom), 20)




png(file.path(svpth, "pcf_global_univ_obs_comb.png"), width = 10, height=5, unit="in", res=300)
par(mfrow=c(1,2))
plot(pcf_univ_obs, main=NULL, lwd=2, axes=F, col=c("black", "firebrick", "grey"), legendpos="topright")
# plot(pcf_univ_obs, legend.only=T, legend.position="topright", plot=F)
axis(1, cex=1.3):axis(2, cex=1.3);box()
plot(pcf_univ_obs_inhom, main=NULL, lwd=2, axes=F, col=c("black", "firebrick", "grey"), legend=F)
axis(1, cex=1.3):axis(2, cex=1.3);box()
dev.off()

#### EXTRAS - in case



sp3 = project(sp1, crs="+proj=moll +datum=WGS84 +units=km +no_defs")


# windows of different dimensions

eco_realm2 = vect("/Users/mattolson/src/drosophila-bioclim/data/dros_spatial/eco-realm2.geojson")
plot(eco_realm2)


gv2 = gv
gv2$realm = extract(eco_realm2[,"WWF_REALM2"], gv2)[,2]
plot(gv2, "realm")
