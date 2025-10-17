
# determine where g is significant

# ASSESS UNIVERSITY STRENGTH
uv0 <- vect(dfr[dfr$Rank<=500,], 
            geom=c("long", "lat"), crs="EPSG:4326")

# gv - 11022
# gv_all. - 63320

# prep data
cross_all = prep_k_data(uv0, gv, sp2)
# run cross
# pcf_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
r_vals <- seq(0, 100, by = 2.5)
pcf_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs", r = r_vals)

set.seed(123)  # For reproducibility
envelope_cross <- envelope(cross_all,
                           fun = pcfcross.inhom,
                           i = "univ",
                           j = "obs",
                           nsim = 99,
                           simulate = expression(rlabel(cross_all)),  # random labeling
                           correction = "translation",
                           savefuns = TRUE,
                           verbose = FALSE,
                           r = r_vals)

plot(envelope_cross, main = "Cross-PCF with Simulation Envelopes")

plot(envelope_cross, main = "", cex.axis=1.2, cex.lab=1.4)

# plot(envelope_cross, xlim=c(0,50), main = "Cross-PCF with Simulation Envelopes")

as.data.frame(envelope_cross)[,c('r','obs','hi')]

# Universities
# Top below is not significant
# Top 100 - 0-25 km was significant (up to 300K times more likely )
# Top 500 - 0-25 km significant with x1328 times likelehood of observing
#.    again 45-100 km up to 96x
# Top 1000 & 2000 not significant


# ASSESS UNIVERSITY STRENGTH
pop0 <- pop[pop$global_ppp_2022_1km_UNadj_constrained>20000,]
pop2
pop20
# prep data
cross_all_p = prep_k_data(pop2, gv_all, sp2)
# run cross
r_vals <- seq(0, 100, by = 2.5)
pcf_cross <- pcfcross.inhom(cross_all_p, i = "univ", j = "obs", r = r_vals)

set.seed(123)  # For reproducibility
envelope_cross_p <- envelope(cross_all_p,
                           fun = pcfcross.inhom,
                           i = "univ",
                           j = "obs",
                           nsim = 99,
                           simulate = expression(rlabel(cross_all_p)),  # random labeling
                           correction = "translation",
                           savefuns = TRUE,
                           verbose = FALSE,
                           r = r_vals)

plot(envelope_cross_p, main = "Cross-PCF with Simulation Envelopes")

as.data.frame(envelope_cross_p)[,c('r','obs','hi')]











uv0_all <- vect(dfr, 
            geom=c("long", "lat"), crs="EPSG:4326")

# nn = nearby(uv0_all, gv_all)

nn = nearby(gv_all, uv0_all)


x_df <- left_join(as.data.frame(uv0_all), as.data.frame(nn), by = c("X" = "id"))

x_df %>% 
  as.data.frame() %>% 
  arrange(desc(id)) %>% 
  select(univ.Institution, id) %>% 
  head(25)

# # PRODUCE P-VALUE
# envelope_cross <- envelope(cross_all,
#                            fun = pcfcross.inhom,
#                            i = "univ",
#                            j = "obs",
#                            nsim = 99,
#                            simulate = expression(rlabel(cross_all)),  # random labeling
#                            correction = "translation",
#                            savefuns = TRUE,
#                            verbose = FALSE,
#                            r = r_vals,
#                            global=T)
# 
# # Extract the maximum absolute deviation of the observed function from the mean
# dev_obs <- max(abs(envelope_cross$obs - envelope_cross$mean))
# 
# # Get the same for each simulation
# dev_sim <- sapply(attr(envelope_cross, "simfuns")[, -1], function(sim) {
#   max(abs(sim - envelope_cross$mean))
# })
# 
# # Calculate p-value
# p_val <- (sum(dev_sim >= dev_obs) + 1) / (length(dev_sim) + 1)
# print(p_val)
