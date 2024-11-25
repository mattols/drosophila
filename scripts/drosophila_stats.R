#
# Area distribution based on species temperature tolerance from WorldClim v2.1
# Drosophila Project - Hjelman & Curnow
# Relations - Olson 10/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # #

# load world political data
library(sf);library(terra);library(dplyr)
library(ggplot2);library(gridExtra);library(tidyterra)

# load output data for plotting
results     <- "~/data/drosophila/results/data_1024"
df_comb     <- read.csv(file.path(results, "drosophila_area_table.csv"))
df_cont     <- read.csv(file.path(results, "drosophila_continent_table.csv"))
df_eco      <- read.csv(file.path(results, "drosophila_eco_table.csv"))


# # # # # # # # # # # # # # # # # #
# function - plot standard metrics
lm_eqn <- function(df, x, y){
  m <- lm(y ~ x, df)
  pval <- summary(m)$coefficients[2, 4]
  r2 <- summary(m)$r.squared
  eq <- substitute( italic(y) == a ~italic(x) + b * "," ~ italic(pvalue) ~ "=" ~ pval * "," ~ italic(r)^2 ~ "=" ~ r2, 
                   list(a = format(unname(coef(m)[2]), digits = 2),
                        b = format(unname(coef(m)[1]), digits = 2),
                        pval = format(pval, digits = 3),
                        r2 = format(r2, digits = 2) ) )
  as.character(as.expression(eq))
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 2.1 MbDna female vs temperature-based area distribution (WorldClim)
## plot with trend & metrics
# create label
data.label <- data.frame(
  x = 25, y = 300,
  label = c(lm_eqn(df_comb, x = df_comb$area_mkm2, y = df_comb$MbDNA_Female))
)
# plot
df_comb %>% ggplot(aes(x = area_mkm2, y = MbDNA_Female)) +
  labs(
    title = "Species MbDNA as a function of spatial distribution",
    subtitle = "based on temperature threshold and WorldClim v2.1",
    colour = "MbDNA (Male)",
    y = "MbDNA (Female)",
    x = bquote("Species distribution (million km"^2*")")
  ) + 
  geom_smooth(method = "lm", color="red", linetype = "dashed", size = 0.5, formula = y ~ x) +
  geom_text(data = data.label, aes(x = x , y = y , label = label), size=4, parse = TRUE) +
  geom_point() + theme_minimal() +
  scale_color_viridis_c()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 6.1 biogeographic distribution count (excluding Antarctic and Oceania)
## plot with trend
# join and pivot data
df_eco2 <- left_join(df_eco, df_comb %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  filter(Biogeographic == "Afrotropic" | Biogeographic == "Neotropic") %>% 
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6)
# create labels
data.label <- df_eco2 %>%
  group_by(Biogeographic) %>%
  summarise(
    x = 10,  
    y = 310,  
    label = lm_eqn(., values, MbDNA_Female)
)
# plot
df_eco2 %>% 
  ggplot(aes(x = values, y = MbDNA_Female, colour=Biogeographic)) +
  labs(title = "Drosophila distribution and genome size",
       subtitle = "by biogeographic region", y = "MbDNA (Female)",
       x = bquote("Area by temperature tolerance (million"~km^2~")")) +
  scale_colour_brewer(palette = "Dark2", guide = 'none') +
  geom_point() + facet_wrap(.~Biogeographic) +
  geom_smooth(method = "lm", color="red", linetype = "dashed", size = 0.5, formula = y ~ x) +
  geom_text(data = data.label, aes(x = x , y = y , label = label), size=3.2, parse = TRUE) +
  theme(legend.text=element_text(size=rel(0.5))) +
  theme_minimal()


# End
