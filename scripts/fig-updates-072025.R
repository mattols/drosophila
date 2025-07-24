#
# New figure changes
# July 2025
#

# load world political data
library(sf);library(terra);library(dplyr)
library(ggplot2);library(gridExtra);library(tidyterra)

# load output data for plotting
results     <- "~/data/drosophila/results/data_1024"
df_comb     <- read.csv(file.path(results, "drosophila_area_table.csv"))
df_cont     <- read.csv(file.path(results, "drosophila_continent_table.csv"))
df_eco      <- read.csv(file.path(results, "drosophila_eco_table.csv"))
eco_regions <- st_read(file.path(results, "eco-realm-summary.geojson"))


# # # # # # # # # # # #
# equations
lm_eqn <- function(df, xvar, yvar){
  formula <- reformulate(xvar, yvar)
  m <- lm(formula, data = df)
  
  # Check if the model actually has a slope term
  coefs <- coef(summary(m))
  if (nrow(coefs) < 2) return("Not enough data")
  
  a <- format(coef(m)[2], digits = 2)
  b <- format(coef(m)[1], digits = 2)
  pval <- format(coefs[2, 4], digits = 3)
  r2 <- format(summary(m)$r.squared, digits = 2)
  
  paste0("y = ", a, "x + ", b, "\n",
         "p = ", pval, "\n",
         "r² = ", r2)
}


# RECREATE FIGURE 3
# join and pivot data
df_eco2 <- left_join(df_eco, df_comb %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  filter(Biogeographic == "Afrotropic" | Biogeographic == "Neotropic") %>%
  # filter(Biogeographic !="X", Biogeographic !="Oceania", Biogeographic !="Antarctic") %>% 
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6)


# Generate data labels (already working version)
data.label <- df_eco2 %>%
  group_by(Biogeographic) %>%
  summarise(
    x = 10,
    y = 310,
    label = lm_eqn(cur_data(), "values", "MbDNA_Female"),
    .groups = "drop"
  )

# Plot
fig03 <- df_eco2 %>%
  ggplot(aes(x = values, y = MbDNA_Female)) +
  geom_point(shape = 1, size = 2) +  # open circles
  geom_smooth(method = "lm", formula = y ~ x,
              linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = data.label, aes(x = x, y = y, label = label),
            size = 5, hjust = 0) +  # increase equation text size
  facet_wrap(~Biogeographic) +
  labs(
    y = "MbDNA (Female)",
    x = bquote("Potential species distribution by temperature tolerance (million"~km^2~")")
  ) +
  theme_classic(base_size = 14) +  # black-and-white theme with larger base font
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )

fig03

# Save the last plot you made
# ggsave("./figs_0725/fig03_bw.png", fig03, width = 9.5, height = 5, dpi = 300)




##### ALL REGIONS

df_eco3 <- left_join(df_eco, df_comb %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  # filter(Biogeographic == "Afrotropic" | Biogeographic == "Neotropic") %>%
  filter(Biogeographic !="X", Biogeographic !="Oceania", Biogeographic !="Antarctic") %>%
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6) 
  # mutate(Biogeographic = ifelse(Biogeographic=='Indo.Malay', 'Indo-Malay', Biogeographic))

# repeat
data.label3 <- df_eco3 %>%
  group_by(Biogeographic) %>%
  summarise(
    x = 7.5,
    y = 285,
    label = lm_eqn(cur_data(), "values", "MbDNA_Female"),
    .groups = "drop"
  )

fig04 <- df_eco3 %>%
  ggplot(aes(x = values, y = MbDNA_Female)) +
  geom_point(shape = 1, size = 2) +  # open circles
  geom_smooth(method = "lm", formula = y ~ x,
              linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = data.label3, aes(x = x, y = y, label = label),
            size = 5, hjust = 0) +  # increase equation text size
  facet_wrap(~Biogeographic) +
  labs(
    y = "MbDNA (Female)",
    x = bquote("Potential species distribution by temperature tolerance (million"~km^2~")")
  ) +
  theme_classic(base_size = 14) +  # black-and-white theme with larger base font
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )

fig04

# ggsave("./figs_0725/fig04_bw.png", fig04, width = 9.5, height = 5, dpi = 300)



##### FIGURE WITH MAP


m1 = fig04
# m1 = df_eco3 %>%
#   ggplot(aes(x = values, y = MbDNA_Female)) +
#   geom_point(shape = 1, size = 2) +  # open circles
#   geom_smooth(method = "lm", formula = y ~ x,
#               linetype = "dashed", color = "black", size = 0.7) +
#   geom_text(data = data.label3, aes(x = x, y = y, label = label),
#             size = 5, hjust = 0) +  # increase equation text size
#   facet_wrap(~Biogeographic) +
#   labs(
#     y = "MbDNA (Female)",
#     x = bquote("Potential species distribution by temperature tolerance (million"~km^2~")")
#   ) +
#   theme_classic(base_size = 14) +  # black-and-white theme with larger base font
#   theme(
#     strip.text = element_text(size = 16, face = "bold"),
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     panel.grid = element_blank()
#   )

# simplify R
library(rmapshaper)
eco_simple = ms_simplify(eco_regions)
eco_simple %>%
  filter(WWF_REALM2 != "Antarctic", WWF_REALM2 != "Oceania" ) %>% 
  st_transform(6933) %>% 
  ggplot() + 
  labs(x = "",y="") +
  geom_sf(aes(fill = WWF_REALM2)) + 
  guides(fill=guide_legend(title="Species count (n=67)")) +
  scale_fill_brewer(palette = "Greys") +
  geom_sf_label(aes(label = species_max, alpha = 0.7))  +
  scale_alpha(guide = 'none') +
  theme_minimal() 

# CHANGES
m1 <- df_eco3 %>%
  ggplot(aes(x = values, y = MbDNA_Female)) +
  geom_point(shape = 1, size = 2) +  # open circles
  geom_smooth(method = "lm", formula = y ~ x,
              linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = data.label3, aes(x = x, y = y, label = label),
            size = 5, hjust = 0) +  # increase equation text size
  facet_wrap(~Biogeographic) +
  labs(
    y = "MbDNA (Female)",
    x = bquote("Potential species distribution by temperature tolerance (million"~km^2~")")
  ) +
  theme_classic(base_size = 14) +  # black-and-white theme with larger base font
  theme(
    strip.text = element_text(size = 16, face = "bold", color = "white"),  # white text for contrast
    strip.background = element_rect(fill = "#636363", color = NA),         # medium grey from "Greys"
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  )

m2 <- eco_simple %>%
  filter(WWF_REALM2 != "Antarctic", WWF_REALM2 != "Oceania") %>% 
  st_transform(6933) %>% 
  ggplot() + 
  labs(x = "", y = "") +
  geom_sf(aes(fill = WWF_REALM2)) + 
  guides(
    fill = guide_legend(
      title = NULL,             # Remove title
      direction = "horizontal", # Horizontal layout
      nrow = 1                  # Force single line
    )
  ) +
  scale_fill_brewer(palette = "Greys") +
  geom_sf_label(aes(label = species_max, alpha = 0.7), size = 5) +  # Increase label size
  scale_alpha(guide = 'none') +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.margin = margin(0, 0, 0, 0)  # Extra space at bottom for annotation
    # plot.caption = element_text(hjust = 0, vjust = 0.9, size = 12)
  ) +
  # labs(subtitle = "Species count (n=67)")
  annotate("text", 
           x = -Inf, y = -Inf, 
           label = "Total species \n counts (n=67)", 
           hjust = -0.41, vjust = -0.7, fontface='bold',
           size = 5)
m2
                    
## plot combined
# png("figs_2025/f6-1-2-MbFem_EcoRegions-Map.png", width = 11, height = 8, unit = "in", res = 300)

png("~/Downloads/map-dros-ecoregion.png", width = 10, height = 10, unit = "in", res = 300)
grid.arrange(m1,m2)
dev.off() 

# ggsave("./figs_0725/fig03_bw.png", fig03, width = 9.5, height = 5, dpi = 300)


# takes an incredibly long time to plot
