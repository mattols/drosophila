

# ASSESS UNIVERSITY STRENGTH
uv0 <- vect(dfr[dfr$Rank<=500,], 
            geom=c("long", "lat"), crs="EPSG:4326")

# prep data
cross_all = prep_k_data(uv0, gv_all, sp2)
# cross_all = prep_k_data(uv0, gv, sp2)

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

plot(envelope_cross, xlim=c(0,50), main = "Cross-PCF with Simulation Envelopes")





# # # # #

# PCF within 6.5 km
barplot(df_cross[df_cross$r<10 & df_cross$r>0,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])
# within next
# # # # #

# PCF within 15 km
barplot(df_cross[df_cross$r<13 & df_cross$r>7,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])

# within next HUGE! (700,000 max)
barplot(df_cross[df_cross$r<20 & df_cross$r>13,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])

# within next (370,000 max)
barplot(df_cross[df_cross$r<26 & df_cross$r>20,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])

# 100,000 max
barplot(df_cross[df_cross$r<34 & df_cross$r>26,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])

# 1,200 max
barplot(df_cross[df_cross$r<40 & df_cross$r>34,'iso'], 
        names.arg=df_cross[df_cross$r<10 & df_cross$r>0,'univ_rank'])
