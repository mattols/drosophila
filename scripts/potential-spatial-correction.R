
#
# New stability test
#
# 

library(spatstat.core)
library(spatstat.geom)
library(terra)

# ---- SETUP PARAMETERS ----
# Choose number of anchor points for testing
anchor_n <- 100  # e.g. same size as your ranked subset

# Number of simulations
nsim <- 199

# ---- STEP 1: OBSERVED CROSS-PCF FOR RANKED U ----
# Subset top N u points
uv_top <- vect(dfr[dfr$Rank <= anchor_n, ], 
               geom = c("long", "lat"), crs = "EPSG:4326")

# Prepare spatstat ppp object (your function)
ppp_obs <- prep_k_data(uv_top, gv, sp2)

# Calculate observed inhomogeneous cross-pcf
pcf_obs <- pcfcross.inhom(ppp_obs, i = "univ", j = "obs", 
                          correction = "trans")

# ---- STEP 2: SIMULATED ENVELOPE FOR RANDOM U SUBSETS ----
# Extract all possible u locations from dfr
all_u <- vect(dfr, geom = c("long", "lat"), crs = "EPSG:4326")

# Initialize matrix to store simulated g(r)
r_vals <- pcf_obs$r
g_sim <- matrix(NA, nrow = length(r_vals), ncol = nsim)

set.seed(42)  # for reproducibility

for (i in 1:nsim) {
  # Randomly sample anchor_n u points
  uv_rand <- all_u[sample(1:nrow(dfr), anchor_n), ]
  
  # Prepare ppp object with randomly chosen u points
  ppp_rand <- prep_k_data(uv_rand, gv, sp2)
  
  # Compute cross-pcf for this random sample
  pcf_sim <- tryCatch({
    pcfcross.inhom(ppp_rand, i = "univ", j = "obs", correction = "trans")
  }, error = function(e) return(NULL))
  
  # Save the g(r) values if successful
  if (!is.null(pcf_sim)) {
    g_sim[, i] <- pcf_sim$trans  # or use $iso if you prefer
  }
}

# ---- STEP 3: CALCULATE ENVELOPE ----
# Compute pointwise envelope (e.g., 5%–95%)
g_low <- apply(g_sim, 1, quantile, probs = 0.025, na.rm = TRUE)
g_high <- apply(g_sim, 1, quantile, probs = 0.975, na.rm = TRUE)

# ---- STEP 4: PLOT OBSERVED VS ENVELOPE ----
plot(r_vals, pcf_obs$trans, type = "l", lwd = 2, col = "blue",
     ylab = "g(r)", xlab = "Distance r (km)",
     main = paste("PCF: Top", anchor_n, "u vs Random Envelope"))
lines(r_vals, g_low, col = "gray", lty = 2)
lines(r_vals, g_high, col = "gray", lty = 2)
legend("topright", legend = c("Observed", "95% envelope"),
       col = c("blue", "gray"), lty = c(1, 2), lwd = c(2, 1))
abline(h = 1, col = "black", lty = 3)


# OPTIONALLY

auc_obs <- sum(pcf_obs$trans[r_vals <= 100]) * diff(r_vals[1:2])
auc_sim <- apply(g_sim[r_vals <= 100, ], 2, function(x) sum(x) * diff(r_vals[1:2]))
p_val <- mean(auc_sim >= auc_obs)  # upper-tail test






# # # # # # # # # #






# # # # # # #
# adjust lamda (density of points)

# Check default kernel estimate for 'obs' pattern
obs_ppp <- split(ppp_obs)$obs
plot(density(obs_ppp), main = "Estimated Intensity of Observations (obs)")


u_ppp <- split(ppp_obs)$univ
plot(density(u_ppp), main = "Estimated Intensity of Anchor Points (u)")

lambda_obs <- density(obs_ppp, sigma = 500)
plot(lambda_obs, main = "Intensity of obs (sigma = 200)")


lambda_obs <- density(obs_ppp, sigma = 200, edge = TRUE)  # example: 200 km bandwidth
lambda_u   <- density(u_ppp, sigma = 200, edge = TRUE)

pcf_obs <- pcfcross.inhom(ppp_obs, i = "univ", j = "obs",
                          lambdaI = lambda_u, lambdaJ = lambda_obs,
                          correction = "trans")

plot(pcf_obs, xlim=c(0,100))




# # # # # # #
# adjust sigma...

# Split ppp_obs into univ and obs components
ppp_split <- split(ppp_obs)
u_ppp     <- ppp_split$univ
obs_ppp   <- ppp_split$obs

# Estimate intensity with sigma = 500 km
sigma_val <- 500  # You can also try 400, 600 etc. to test sensitivity

lambda_u   <- density(u_ppp, sigma = sigma_val, edge = TRUE)
lambda_obs <- density(obs_ppp, sigma = sigma_val, edge = TRUE)

pcf_obs <- pcfcross.inhom(
  X = ppp_obs,
  i = "univ", j = "obs",
  lambdaI = lambda_u,
  lambdaJ = lambda_obs,
  correction = "trans"
)

simulate_fn <- function(i) {
  uv_rand <- all_u[sample(1:nrow(dfr), anchor_n), ]
  prep_k_data(uv_rand, gv, sp2)
}

env <- envelope(
  Y = ppp_obs,
  fun = pcfcross.inhom,
  simulate = simulate_fn,
  i = "univ", j = "obs",
  nsim = 199,
  correction = "trans",
  lambdaI = lambda_u,
  lambdaJ = lambda_obs,
  savefuns = TRUE,
  verbose = TRUE
)

plot(env, main = paste("Envelope (σ =", sigma_val, "km): Top", anchor_n, "u vs Random"))
abline(h = 1, lty = 3, col = "black")

