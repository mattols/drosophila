





# ASSESS UNIVERSITY STRENGTH
uv0 <- vect(dfr[dfr$Rank<=500,], 
            geom=c("long", "lat"), crs="EPSG:4326")

# prep data
cross_all = prep_k_data(uv0, gv, sp2)

# run cross
pcf_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")

set.seed(123)  # For reproducibility
envelope_cross <- envelope(cross_all,
                           fun = pcfcross.inhom,
                           i = "univ",
                           j = "obs",
                           nsim = 99,
                           simulate = expression(rlabel(cross_all)),  # random labeling
                           correction = "translation",
                           savefuns = TRUE,
                           verbose = FALSE)

plot(envelope_cross, main = "Cross-PCF with Simulation Envelopes")

pcf_cross_data <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
plot(pcf_cross$r, pcf_cross$numerator, type = "h",
     main = "Number of Point Pairs per Distance Bin",
     xlab = "r (km)", ylab = "Number of Pairs")



