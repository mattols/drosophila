

library(terra);library(dplyr)

library(geodata)
if (!file.exists("data")){dir.create("data")}
outpath = "data"
tmin <- worldclim_global('tmin',res=10,path=outpath)
plot(tmin) # monthly temperatures

bioclim19 <- worldclim_global('bio', res=10, path=outpath)
nlyr(bioclim19) # dimensions of data (layers)
ext(bioclim19)
names(bioclim19)

# subset variables of interest
tmax <- subset(bioclim19, 5);names(tmax) <- "BIO5 - Max Temperature of Warmest Month"
tmin <- subset(bioclim19, 6);names(tmin) <- "BIO6 - Min Temperature of Warmest Month"
plot(c(tmin,tmax), col=rev(heat.colors(15)))

# read in data - omit NA values
df <- read.csv("climate_data_drosophila_June.csv")
df2 <- df %>% dplyr::select(Species,Subgenus,Tmax,Tmin) %>% na.omit()
head(df2)

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
# summarize for one
df.stats <- data.frame(matrix(nrow=nrow(df2),ncol=4), row.names = c(""))
crds(xxcomb[[1]], df=TRUE) %>% filter(values(xxcomb[[1]])==1) %>% summarize(n=n(), latmean=mean(y), latmax=max(y), latmin=min(y) )
# xd <- sapply(1:nlyr(xxcomb), function(x) (crds(xxcomb[[x]], df=TRUE) %>% filter(values(xxcomb[[x]])==1) %>% summarize(n=n(), latmean=mean(y), latmax=max(y), latmin=min(y) )) )
lyr_summary <- function(lyr){
  rowinfo <- crds(lyr, df=TRUE) %>% filter(values(lyr)==1) %>% summarize(n=n(), latmean=mean(y), latmax=max(y), latmin=min(y) )
  return(rowinfo)
}

xd <- sapp(xxcomb, function(x) lyr_summary(x)  ) 
xd2 <- sapply(1:nlyr(xxcomb), function(x) lyr_summary(xxcomb[[x]])  ) 

df3 <- as.data.frame(t(xd2)) %>% mutate(Species = df2$Species, .before=n) %>% full_join(df2, by='Species') %>%
  mutate(n = unlist(n), latmean = unlist(latmean),latmax = unlist(latmax),latmin = unlist(latmin))



# plot first species with world map
worldmap <- world(resolution=5,level=0,outpath,version="latest")
plot(subset(t.dist,1), col=rev(heat.colors(15)), main=paste("Temperature threshold for",df2$Species[1]) )
plot(worldmap,add=TRUE, lwd=0.5, border=alpha(rgb(0,0,0), 0.5))

# you can check out some other biotic factors
# BIO12 = Annual Precipitation

if (bounds!="lower" | bounds!="upper"){stop("bounds should be either: 'lower' or 'upper'")}
