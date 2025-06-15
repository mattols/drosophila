#
# L/K at multiple levels
#
#

library(terra);library(sf)
library(spatstat)
library(dplyr)

# pth = "/Users/mattolson/src/drosophila-bioclim/data/output/all.drosophila.gbif.csv"
eco_realm2 = vect("/Users/mattolson/src/drosophila-bioclim/data/dros_spatial/eco-realm2.geojson")
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
# gv = gv[gv$eco=="Oceania",1]

# # # # # # # # # # # # # # # # # #

## BY ECO REALM

# select one ecorealm
uv <- vect(dfr[dfr$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")
ue_names = setdiff(unique(eco_realm2$WWF_REALM2), c("Antarctic", "Oceania"))
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


# BY DISTANCE TO UNIVERSITY
# 100, 500, 1000, 1500
# univ_lists = c(10, 25, 50, 100, 500, 1000, 1500, dim(dfr)[1]+1) # Inf) 
univ_lists = c(10, 25, 50, 100, 500, 1000, 1500)  
# cross_list = list()
# df_cross = data.frame()
for (rank_level in univ_lists){
  # rank_level = 10
  print(paste("...running rank:", rank_level,"or", match(rank_level, univ_lists), "of", length(univ_lists)))
  uv0 <- vect(dfr[dfr$Rank<=rank_level, ], 
              geom=c("long", "lat"), crs="EPSG:4326")
  print(dim(uv0))
  # prep data
  cross_all = prep_k_data(uv0, gv, sp2)

  # run cross
  L_obs_cross <- Lcross(cross_all, i = "univ", j = "obs")
  # cross_list[match(rank_level, univ_lists)] <- L_obs_cross
  
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


# saveRDS(df_cross, "/Users/mattolson/data/drosophila/results/k/df_cross_1500.rds")
df_cross = readRDS("/Users/mattolson/data/drosophila/results/k/df_cross_1500.rds")

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






library(ggplot2)
library(tidyr)
library(dplyr)

# Example: assume df is your dataframe
# df <- read.csv("your_file.csv")

# Reshape data from wide to long format
df_long <- df_cross %>%
  pivot_longer(cols = c(theo, trans, border),
               names_to = "variable",
               values_to = "value")

# Set line types
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "border" = "dashed"
)

# Plot
ggplot(df_long, aes(x = r, y = value, color = univ_rank, linetype = variable)) +
  geom_line() +
  scale_linetype_manual(values = line_types) +
  labs(x = "r", y = "Value", color = "University Rank", linetype = "Variable") +
  theme_minimal()



# Example: assume df is your dataframe
# df <- read.csv("your_file.csv")

# Reshape data from wide to long format
df_long <- df_cross %>%
  pivot_longer(cols = c(theo, trans, border),
               names_to = "variable",
               values_to = "value")

# Set line types for each variable
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "border" = "dashed"
)

# Plot with grouping by both variable and univ_rank
ggplot(df_long[complete.cases(df_long$variable),], aes(x = r, y = value,
                    color = as.factor(univ_rank),
                    linetype = variable,
                    group = interaction(univ_rank, variable))) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = line_types) +
  labs(x = "r", y = "Value", color = "University Rank", linetype = "Variable") +
  theme_minimal() #+
  # facet_wrap(~univ_rank)


  
  
  
  
  
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Reshape to long format
df_long <- df_cross %>%
  # pivot_longer(cols = c(theo, trans, border),
  # pivot_longer(cols = c(theo, trans),
  pivot_longer(cols = c(theo, border),
               names_to = "variable",
               values_to = "value")

# 2. Normalize: subtract the value at r = 0 within each group
df_long_norm <- df_long %>%
  group_by(univ_rank, variable) %>%
  mutate(value_norm = value - value[r == 0]) %>%
  ungroup()

# 3. Define linetypes
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "border" = "dashed"
)

# 4. Plot
ggplot(df_long_norm, aes(x = r, y = value_norm,
                         color = univ_rank,
                         linetype = variable,
                         group = interaction(univ_rank, variable))) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = line_types) +
  labs(x = "r", y = "Normalized Value (centered at r = 0)",
       color = "University Rank", linetype = "Model") +
  theme_minimal()










# # # # # # # # # # #
library(tidyr);library(ggplot2)


######### HERERERE

# 1. Compute relative values to theo
df_relative1 <- df_cross %>%
  mutate(trans = trans - theo,
         border = border - theo,
         iso = iso - theo,
         theo = 0)  # Set theo as the baseline (0)

# 2. Pivot to long format
df_long1 <- df_relative1 %>%
  # pivot_longer(cols = c(theo, trans, border),
  # pivot_longer(cols = border,
  pivot_longer(cols = iso,
               names_to = "variable",
               values_to = "value")

# 3. Line type mapping
line_types <- c(
  "theo" = "dotted",
  "trans" = "solid",
  "border" = "dashed"
)

df_long1 <- df_long1 %>%
  mutate(univ_rank = as.factor(univ_rank))

# Plot
p1 = ggplot(df_long1, aes(x = r, y = value,
                    color = univ_rank)) +
                    # linetype = variable,
                    # group = interaction(univ_rank, variable))) +
  # geom_line(linewidth = 1) +
  geom_line() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Option 1 (automatic greyscale)
  scale_color_grey(start = 0.1, end = 0.7) +
  # Option 2 (comment out one of these)
  # scale_color_manual(values = c("#000000", "#4D4D4D", "#7F7F7F", "#B2B2B2", "#D9D9D9")) +
  # scale_linetype_manual(values = line_types) +
  labs(x = "r (km)",
       # y = expression(hat(L)),
       y = expression(hat(L)[obs/univ]),
       # y = "L̂[observations near universities]",
       # y = "L̂ (univ,obs)",
       # y = expression(hat(L)[observations * " near " * universities]),
       color = "Rank") +
  # xlim(0, 1750) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),  # X-axis tick labels
    axis.text.y = element_text(size = 14),   # Y-axis tick labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    # legend.position = "inside",
    legend.position = c(0.1,0.75),
    legend.background = element_rect(fill = "white", color = "black"))

p1

svpth = "~/src/drosophila/figs_0625"
ggsave(file.path(svpth, "L_iso_rank.png"), p1, width = 5, height = 6, dpi = 300, bg = "white")


  


# # # # # # # # # # 
### ECO FIGS



# # # # # # # # # # #

# 1. Compute relative values to theo
df_relative <- df_eco_cross %>%
  mutate(trans = trans - theo,
         border = border - theo,
         iso = iso - theo,
         theo = 0)  # Set theo as the baseline (0)

# 2. Pivot to long format
df_long <- df_relative %>%
  # pivot_longer(cols = c(theo, trans, border),
  # pivot_longer(cols = border,
  # pivot_longer(cols = trans,
  pivot_longer(cols = iso,
               names_to = "variable",
               values_to = "value")


df_long <- df_long %>%
  mutate(realm = as.factor(realm))

# Plot
p2 = ggplot(df_long, aes(x = r, y = value,
                         color = realm)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "r (km)",
       y = expression(hat(L)[obs/univ]),
       color = "Region") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),  # X-axis tick labels
    axis.text.y = element_text(size = 14),   # Y-axis tick labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    # legend.position = "inside",
    legend.position = c(0.15,0.80),
    legend.background = element_rect(fill = "white", color = "black"))

p2

svpth = "~/src/drosophila/figs_0625"
ggsave(file.path(svpth, "L_iso_500u_bioregion.png"), p2, width = 5, height = 6, dpi = 300, bg = "white")



library(patchwork)

# Combine side-by-side
combined_plot <- p1 + p2

# Save
# ggsave("combined_plot.png", combined_plot, width = 12, height = 6, dpi = 300)
svpth = "~/src/drosophila/figs_0625"
ggsave(file.path(svpth, "L_iso_univ_bio.png"), combined_plot, width = 12, height = 6, dpi = 300, bg = "white")

#






df_long %>%  filter(realm != "Antarctic", realm != "Oceania" ) %>%
  mutate(realm = as.factor(realm)) %>% 
  ggplot(aes(x = r, y = value,
                    color = realm)) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(aes(color = realm)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

  labs(x = "r (km)",
       y = expression(hat(L)[obs/univ]),
       color = "Region") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),  # X-axis tick labels
    axis.text.y = element_text(size = 14),   # Y-axis tick labels
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    # legend.position = "inside",
    legend.position = c(0.15,0.80),
    legend.background = element_rect(fill = "white", color = "black"))
 

