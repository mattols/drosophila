#
# Vary pcf by population & university
# Do we see a stronger pattern with universities vs population?
#

library(terra);library(sf)
library(spatstat)
library(dplyr)

# # # # # # # # # SET UP DATA # # # # # # # # #

# pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.drosophila.gbif.csv"
eco_realm2 = vect("/Users/mattolson/src/drosophila-bioclim/data/dros_spatial/eco-realm2.geojson")
pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.locale.gbif.csv"
pth2 = "/Users/mattolson/src/drosophila-bioclim/data/output/university_loc_rankings.csv"

# not big enough (n=2540)
# mcty <- vect('~/data/drosophila/world_pop/World_Cities.geojson')

# population data
popr <- rast('~/data/drosophila/world_pop/global_ppp_2022_1km_UNadj_constrained.tif')
pop <- as.points(popr)
pop2 = pop[pop$global_ppp_2022_1km_UNadj_constrained>10000,]
pop20 = pop[pop$global_ppp_2022_1km_UNadj_constrained>20000,]
#
# plot(sp2, col=adjustcolor('navajowhite', 0.3), main="Population")
# plot(project(pop3, sp2),"global_ppp_2022_1km_UNadj_constrained", add=T)

eco_simple2 = vect(ms_simplify(st_as_sf(eco_realm2)))

# GBIF DATA
dfg = read.csv(pth)
gv <- vect(dfg, geom=c("long", "lat"), crs="EPSG:4326")

# UNIVERSITY DATA
dfr = read.csv(pth2)
head(dfr);dim(dfr)
uv <- vect(dfr[dfr$Rank<=1000, ], geom=c("long", "lat"), crs="EPSG:4326")
dim(uv)

# reprojected global window
create_window_global <- function(global=TRUE){
  r = rast("~/data/drosophila/worldclim_05/wc2.1_10m/wc2.1_10m_bio_1.tif")
  ry = patches(r)
  ry[ry==1 | ry>=1634] = NA
  rz = zonal(cellSize(ry,unit="km"), ry, sum, as.raster=TRUE)
  # remove small patches
  s = ifel(rz < 1.2e5, NA, 1)
  sp2 = project(as.polygons(s), "EPSG:6933")
  return(sp2)
}

# SHOW ALL (point comparison)
sp2 <- create_window_global()
par(mfrow=c(3,1))
plot(sp2, col=adjustcolor('navajowhite', 0.5), main="Population > 10,000")
plot(project(pop2[,1], sp2), add=T, cex=0.5, pch=18, col='deepskyblue3')
plot(sp2, col=adjustcolor('navajowhite', 0.5), main="Universities (Top 1000)")
plot(project(uv[,1], sp2),add=T, cex=0.5, pch=17, col='firebrick')
plot(sp2, col=adjustcolor('navajowhite', 0.5), main="GBIF Observations")
plot(project(gv[,1], sp2),add=T, cex=0.7, pch=1)
par(mfrow=c(1,1))




# # # # # # # # # TEST RUNS # # # # # # # # #


# TEST for 1000
uv0 <- vect(dfr[dfr$Rank<=1000, ], 
            geom=c("long", "lat"), crs="EPSG:4326")
# prep data
cross_all = prep_k_data(uv0, gv, sp2)
cross_all_pop = prep_k_data(pop3, gv, sp2)

# CLUSTERING FOR UNV \ DISPERSION FOR POP
# run cross
L_univ_obs <- Lcross(cross_all, i = "univ", j = "obs")
strt = Sys.time()
L_pop_obs <- Lcross(cross_all_pop, i = "univ", j = "obs")
Sys.time() - strt
# 20000 pop- 7.32 mins 

strt = Sys.time()
cross_all_pop2 = prep_k_data(pop2, gv, sp2)
L_pop_obs2 <- Lcross(cross_all_pop2, i = "univ", j = "obs")
Sys.time() - strt 
# 10000 pop - 37.76545 mins (for single)

plot(L_pop_obs)
plot(L_univ_obs)
plot(L_pop_obs2)

# CHECK OUT PCF
pcf_univ_obs <- pcfcross(ppp_all, i = "univ", j = "obs")
pcf_univ_obs_inhom <- pcfcross.inhom(ppp_all, i = "univ", j = "obs")


strt = Sys.time()
cross_univ0 = prep_k_data(uv0, gv, sp2)
pcf_univ <- pcfcross(cross_univ0, i = "univ", j = "obs")
Sys.time() - strt # 1.706967 mins
cross_pop0 = prep_k_data(pop3, gv, sp2)
pcf_pop <- pcfcross(cross_pop0, i = "univ", j = "obs")
Sys.time() - strt # 7.6 mins

strt = Sys.time()
pcf_popi <- pcfcross.inhom(cross_pop0, i = "univ", j = "obs")
Sys.time() - strt # 1.13

plot(pcf_pop)



# # # # # # # # # FUNCTIONS # # # # # # # # #

prep_k_data <- function(univ_vect, gbif_vect, sp2){
  #
  # required
  #   uv - university vector object
  #   gv - gbif spatial object
  #   sp2 - window object
  # 
  
  # project and crop to extent
  univ = project(univ_vect, sp2)[sp2,]
  gbif = project(gbif_vect, sp2)[sp2,]
  
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
  my_window <- as.owin(st_as_sf(sp2))  
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
  
  return(ppp_all)
}



# # # # # # # # # Rank Runs # # # # # # # # #

# BY DISTANCE TO UNIVERSITY
univ_lists = c(10, 50, 100, 500, 1000, 1500)  
#
# LOOP
for (rank_level in univ_lists){
  # rank_level = 10
  print(paste("...running rank:", rank_level,"or", match(rank_level, univ_lists), "of", length(univ_lists)))
  uv0 <- vect(dfr[dfr$Rank<=rank_level, ], 
              geom=c("long", "lat"), crs="EPSG:4326")
  print(dim(uv0))
  # prep data
  cross_all = prep_k_data(uv0, gv, sp2)
  
  # run cross
  # L_obs_cross <- Lcross(cross_all, i = "univ", j = "obs")
  L_obs_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
  
  L_obs_cross = as.data.frame(L_obs_cross)
  L_obs_cross$univ_rank = rank_level
  
  if(match(rank_level, univ_lists)==1){
    df_cross = L_obs_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_cross = rbind(df_cross, L_obs_cross)
  }
}
# names(L_obs_cross)
head(df_cross);dim(df_cross)
plot(df_cross)

# saveRDS(df_cross, "/Users/mattolson/data/drosophila/results/k/dfcross_pop20_pcf.rds")

# saveRDS(df_cross, "/Users/mattolson/data/drosophila/results/k/dfcross_univ_pcf.rds")
df_cross = readRDS("/Users/mattolson/data/drosophila/results/k/dfcross_univ_pcf.rds")

# 1. Compute relative values to theo

# 2. Pivot to long format
df_long1 <- df_cross %>%
  pivot_longer(cols = !c(r, univ_rank),
               names_to = "variable",
               values_to = "value")

# 3. Line type mapping
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "iso" = "dashed"
)

df_long1 <- df_long1 %>%
  mutate(univ_rank = as.factor(univ_rank))

# Plot
df_long1 %>% 
  filter(variable=='trans') %>%
  filter(univ_rank%in%c(10, 100, 1000)) %>%
  # filter(univ_rank==10) %>%
  ggplot( aes(x = r, y = value,
                          color = univ_rank)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Option 1 (automatic greyscale)
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)",
       y = expression(hat(g)[obs/univ]),
       color = "Rank") +
  ylim(c(0,1e5)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),  # X-axis tick labels
    axis.text.y = element_text(size = 14),   # Y-axis tick labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    # legend.position = "inside",
    legend.position = c(0.8,0.85),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))
  # facet_wrap(.~univ_rank, ncol=1)


svpth = "~/src/drosophila/figs_0725/"
ggsave(file.path(svpth, "g_iso_rank.png"), p1, width = 5, height = 6, dpi = 300, bg = "white")



# # # # # # # # # Rank Runs # # # # # # # # #

# BY DISTANCE TO POPULATION
pop_lists = c(10, 20, 30, 40)  *1e3
#
# LOOP
for (pop_level in pop_lists){
  print(paste("...running rank:", pop_level,"or", match(pop_level, pop_lists), "of", length(pop_lists)))
  pop0 <- pop[pop$global_ppp_2022_1km_UNadj_constrained>pop_level,]
  print(dim(uv0))
  # prep data
  # cross_all = prep_k_data(uv0, gv, sp2)
  cross_all = prep_k_data(pop0, gv, sp2)
  
  # run cross
  # L_obs_cross <- Lcross(cross_all, i = "univ", j = "obs")
  L_obs_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
  
  L_obs_cross = as.data.frame(L_obs_cross)
  L_obs_cross$univ_rank = pop_level
  
  if(match(pop_level, pop_lists)==1){
    df_cross_pop = L_obs_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_cross_pop = rbind(df_cross_pop, L_obs_cross)
  }
}
# names(L_obs_cross)
head(df_cross_pop);dim(df_cross_pop)

# saveRDS(df_cross_pop, "/Users/mattolson/data/drosophila/results/k/dfcross_pop_pcf.rds")


# 1. Compute relative values to theo

# 2. Pivot to long format
df_long1 <- df_cross_pop %>%
  pivot_longer(cols = !c(r, univ_rank),
               names_to = "variable",
               values_to = "value")

# 3. Line type mapping
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "iso" = "dashed"
)

df_long1 <- df_long1 %>%
  mutate(univ_rank = as.factor(univ_rank))

# Plot
df_long1 %>% 
  filter(variable=='trans') %>%
  filter(univ_rank%in%c(10000, 20000)) %>%
  filter(univ_rank==10000) %>%
  ggplot( aes(x = r, y = value,
              color = univ_rank)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  # Option 1 (automatic greyscale)
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)",
       y = expression(hat(g)[obs/univ]),
       color = "Population") +
  # ylim(c(0,1e5)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),  # X-axis tick labels
    axis.text.y = element_text(size = 14),   # Y-axis tick labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    # legend.position = "inside",
    legend.position = c(0.8,0.85),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))
# facet_wrap(.~univ_rank, ncol=1)





# # # # # # # # # ECO Runs # # # # # # # # #
# classify points in eco
gv$eco = NA
for(ue in unique(eco_realm2$WWF_REALM2)){
  # gv[eco_realm2[eco_realm2$WWF_REALM2=="Australasia"],]
  pts_nam = gv[eco_realm2[eco_realm2$WWF_REALM2==ue],]
  gv$eco[pts_nam$X] = ue
}

# FIXED at 1000 Universities (10, 100, 1000)
uv <- vect(dfr[dfr$Rank<=1000, ], geom=c("long", "lat"), crs="EPSG:4326")
ue_names = setdiff(unique(eco_realm2$WWF_REALM2), c("Antarctic", "Oceania"))
#
# LOOP
for(ue in ue_names){
  
  print(paste("...running realm:", ue,"or", match(ue, unique(eco_realm2$WWF_REALM2)), "of", length(unique(eco_realm2$WWF_REALM2))))
  
  # subsets
  uv_sub = uv[eco_realm2[eco_realm2$WWF_REALM2==ue],]
  gv_sub = gv[eco_realm2[eco_realm2$WWF_REALM2==ue],]
  
  # prep data
  cross_eco_all = prep_k_data(uv_sub, gv_sub, sp2)
  
  # run cross
  L_eco_cross <- Lcross(cross_eco_all, i = "univ", j = "obs")
  # cross_list[match(rank_level, univ_lists)] <- L_obs_cross
  
  L_eco_cross = as.data.frame(L_eco_cross)
  L_eco_cross$realm = ue
  
  if(match(ue, unique(eco_realm2$WWF_REALM2))==1){
    df_eco_cross = L_eco_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_eco_cross = rbind(df_eco_cross, L_eco_cross)
  }
  
}
head(df_eco_cross)

# saveRDS(df_eco_cross, "/Users/mattolson/data/drosophila/results/k/df_eco_cross.rds")
df_cross = readRDS("/Users/mattolson/data/drosophila/results/k/df_eco_cross.rds")


# PLOTS
df_eco_cross
