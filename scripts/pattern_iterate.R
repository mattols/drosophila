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
# gv = gv[gv$eco=="Oceania",1]

# # # # # # # # # # # # # # # # # #

## BY ECO REALM

# select one ecorealm
uv <- vect(dfr[dfr$Rank<=100, ], geom=c("long", "lat"), crs="EPSG:4326")
for(ue in unique(eco_realm2$WWF_REALM2)){
  
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



# BY DISTANCE TO UNIVERSITY
# 100, 500, 1000, 1500
univ_lists = c(10, 25, 50, 100, 500, 1000, 1500, dim(dfr)[1]+1) # Inf) 
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
  facet_wrap(~univ_rank)


