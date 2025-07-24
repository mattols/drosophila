#
# Obs & University Maps
#
#
#

# load world political data
library(sf);library(terra);library(dplyr)
library(ggplot2);library(gridExtra);library(tidyterra)

pth = "~/Downloads/all.drosophila.gbif.csv"
pth = "~/Downloads/all.locale.gbif.csv"
dfg = read.csv(pth)
head(dfg);dim(dfg)
gv <- vect(dfg, geom=c("long", "lat"), crs="EPSG:4326")


# load output data for plotting
results     <- "../../data/drosophila/results/data/"
dfgs        <- read.csv(file.path(results, "dros_gs_occur_june.csv")) # OLD DATA
world       <- read_sf(file.path(results, "world-simple.geojson"))
dros_heat   <- rast(file.path(results, "drosophila_heatmap_n67_10m_1.grd"))


# read university data
pth2 = "~/Downloads/university_loc_rankings.csv"
dfr = read.csv(pth2)
# dfr = dfr[dfr$Rank<=1000, ]
head(dfr);dim(dfr)
uv <- vect(dfr[dfr$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")
univ <- st_as_sf(uv)





# Figure 1.1 global species w/ heatmap heatmap
# point map
p1 <- world %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5)) +
  geom_point(data = dfgs, aes(x = long, y = lat, colour = GS), size = 0.5) +
  labs(
    tag = "A",
    fill = "Species Occurences",
    title = "GBIF observations (subset)") +
  scale_color_viridis_c() + 
  ylim(c(-60,85)) + theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))
# universities
p2 <- world %>% ggplot() +
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_point(data = dfr, aes(x = long, y = lat, colour = Rank), size = 0.5) +
  labs(
    tag = "B",
    fill = "Rank",
    title = "Universities Ranked (All)") +
  scale_color_viridis_c() + 
  ylim(c(-60,85)) + theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))
p3 <- world %>% ggplot() +
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_point(data = dfr[dfr$Rank<=500, ], aes(x = long, y = lat, colour = Rank), size = 0.5) +
  labs(
    tag = "C",
    fill = "Rank",
    title = "Universities Ranked (Top 500)") +
  scale_color_viridis_c() + 
  ylim(c(-60,85)) + theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))
# adjust margins to align
p1_adjusted <- p1 + theme(plot.margin = margin(t = 0, r = 18, b = 10, l = 0))  
p2_adjusted <- p2 + theme(plot.margin = margin(t = 0, r = 0, b = 10, l = 5)) 
p3_adjusted <- p3 + theme(plot.margin = margin(t = 0, r = 0, b = 10, l = 5)) 
## plot combined with layout adjustment

png("~/Downloads/map-univ-obs-v01.png", width = 8, height = 10, unit = "in", res = 300)
# grid.arrange(p1_adjusted, p2_adjusted, ncol = 1, heights = c(1, 1) )
grid.arrange(p1_adjusted, p2_adjusted, p3_adjusted, ncol = 1, heights = c(1, 1, 1) )
dev.off()
