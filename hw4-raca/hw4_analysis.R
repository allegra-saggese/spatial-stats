# HW-4: Bayesian Gaussian Process — Soil Carbon  (MD / NJ / DE / PA)
# ─────────────────────────────────────────────────────────────────────
# HW question mapping:
#   Q1 → Sections 1-3  : Data extraction + EDA (MD/NJ/DE/PA, measured only)
#   Q2 → Section 4     : Full Bayesian GP + MCMC; covariate significance
#   Q3 → Section 4     : Priors (Jeffreys 1/σ², reference 1/φ, δ from REML)
#   Q4 → Section 5     : Posterior predictive surface (mean + uncertainty)
#   Q5 → Section 6     : MAP estimate; comparison to MCMC
#   Q6 → Section 7     : Goodness of fit (LOO-CV, PIT, CRPS)

# ── 0. Setup ──────────────────────────────────────────────────────────────────
library(soilDB); library(geoR); library(readr); library(dplyr)
library(ggplot2); library(maps); library(sf); library(terra); library(httr2)

setwd("/Users/allegrasaggese/Documents/GitHub/spatial-stats/hw4-raca")
dir.create("data_clean", showWarnings = FALSE)
dir.create("figures",    showWarnings = FALSE)

STATES      <- c("MD", "NJ", "DE", "PA")
STATES_FULL <- c("maryland", "new jersey", "delaware", "pennsylvania")
KAPPA       <- 1.5          # Matérn smoothness (best from HW-3 analysis)
N_BURN      <- 2000
N_ITER      <- 12000
THIN        <- 2
RAW_FILE    <- "data_clean/MNJDEPA_raw.csv"
CLEAN_FILE  <- "data_clean/MNJDEPA_measured_clean.csv"
MCMC_FILE   <- "data_clean/mcmc_samples.rds"

# ── 1. Data download (skips if cached)  [HW Q1] ──────────────────────────────
if (!file.exists(RAW_FILE)) {
  fetch_state <- function(st) {
    cat(sprintf("Fetching RaCA: %s\n", st))
    soc  <- fetchRaCA(state = st)
    samp <- soc$sample |> select(any_of(c("sample_id","rcapid","soc","soc_sd",
                                           "soc_measured","sample_top","sample_bottom","texture")))
    site <- soc$pedons@site |>
      select(any_of(c("rcapid","elev_field","x","y","landuse"))) |>
      rename(elevation = elev_field, long = x, lat = y)
    samp |> left_join(site, by = "rcapid") |>
      filter(sample_top == 0) |> mutate(state = st)
  }
  raw <- bind_rows(lapply(STATES, fetch_state))

  # Fill missing elevation with SRTM via OpenTopography
  n_miss <- sum(is.na(raw$elevation))
  if (n_miss > 0) {
    cat(sprintf("Fetching SRTM elevation for %d missing rows...\n", n_miss))
    api_key <- Sys.getenv("OPENTOPO_API_KEY")
    if (!nchar(api_key)) stop("OPENTOPO_API_KEY not in ~/.Renviron")
    dem_file <- tempfile(fileext = ".tif")
    request("https://portal.opentopography.org/API/globaldem") |>
      req_url_query(demtype = "SRTMGL1",
                    west  = floor(min(raw$long, na.rm=TRUE)*10)/10 - 0.1,
                    south = floor(min(raw$lat,  na.rm=TRUE)*10)/10 - 0.1,
                    east  = ceiling(max(raw$long, na.rm=TRUE)*10)/10 + 0.1,
                    north = ceiling(max(raw$lat,  na.rm=TRUE)*10)/10 + 0.1,
                    outputFormat = "GTiff", API_Key = api_key) |>
      req_timeout(600) |> req_perform() |>
      resp_body_raw() |> writeBin(dem_file)
    srtm <- terra::extract(rast(dem_file),
                           vect(raw, geom = c("long","lat"), crs = "EPSG:4326"))[,2]
    raw  <- raw |> mutate(elevation = ifelse(is.na(elevation), srtm, elevation))
    unlink(dem_file)
  }
  write_csv(raw, RAW_FILE)
  cat(sprintf("Saved raw data: %d rows\n", nrow(raw)))
}

# ── 2. Prepare  [HW Q1] ───────────────────────────────────────────────────────
if (!file.exists(CLEAN_FILE)) {
  raw <- read_csv(RAW_FILE, show_col_types = FALSE)
  dat <- raw |>
    filter(soc_measured == "measured", !is.na(soc), soc > 0,
           !is.na(long), !is.na(lat)) |>
    select(sample_id, rcapid, soc, soc_sd, soc_measured,
           depth_top = sample_top, depth_bot = sample_bottom,
           texture, elevation_m = elevation,
           lon = long, lat, landuse, state)
  ref_lon <- mean(dat$lon); ref_lat <- mean(dat$lat)
  dat <- dat |> mutate(
    x_km    = (lon - ref_lon) * cos(ref_lat * pi / 180) * 111.32,
    y_km    = (lat - ref_lat) * 111.32,
    lsoc    = log(soc),
    elev_sc = as.numeric(scale(elevation_m))
  )
  write_csv(dat, CLEAN_FILE)
}

dat <- read_csv(CLEAN_FILE, show_col_types = FALSE) |>
  mutate(lu = factor(landuse))

cat(sprintf("n = %d measured observations across %s\n",
            nrow(dat), paste(STATES, collapse = "/")))
cat("Land use:\n"); print(table(dat$landuse))

# ── 3. EDA  [HW Q1] ───────────────────────────────────────────────────────────
border <- map_data("state", region = STATES_FULL)

# 3a. Site map
p_map <- ggplot() +
  geom_polygon(data = border, aes(long, lat, group = group),
               fill = "grey92", colour = "grey55", linewidth = 0.3) +
  geom_point(data = dat, aes(lon, lat, colour = lsoc), size = 2.2, alpha = 0.85) +
  scale_colour_viridis_c(name = "log(SOC)") +
  coord_fixed(1.3) +
  labs(title = "RaCA measured sites — MD/NJ/DE/PA", x = "Longitude", y = "Latitude") +
  theme_bw(11)
ggsave("figures/01_map.png", p_map, width = 7, height = 5.5, dpi = 150)

# 3b. Normality check
png("figures/02_normality.png", width = 900, height = 420, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(dat$soc,  breaks = 30, col = "steelblue3", border = "white", freq = FALSE,
     main = "SOC (g/kg)", xlab = "SOC")
curve(dnorm(x, mean(dat$soc), sd(dat$soc)), add = TRUE, col = "firebrick", lwd = 2)
hist(dat$lsoc, breaks = 30, col = "steelblue3", border = "white", freq = FALSE,
     main = "log(SOC)", xlab = "log(SOC)")
curve(dnorm(x, mean(dat$lsoc), sd(dat$lsoc)), add = TRUE, col = "firebrick", lwd = 2)
dev.off()

# 3c. Covariates
p_lu <- ggplot(dat, aes(lu, lsoc, fill = lu)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1) + scale_fill_brewer(palette = "Set2") +
  labs(title = "log(SOC) by land use", x = "Land use", y = "log(SOC)") +
  theme_bw(11) + theme(legend.position = "none")
ggsave("figures/03_landuse.png", p_lu, width = 5, height = 4, dpi = 150)

p_elev <- ggplot(dat, aes(elevation_m, lsoc)) +
  geom_point(alpha = 0.55, colour = "steelblue4") +
  geom_smooth(method = "lm", colour = "firebrick", se = TRUE) +
  labs(title = "log(SOC) vs elevation", x = "Elevation (m)", y = "log(SOC)") +
  theme_bw(11)
ggsave("figures/04_elevation.png", p_elev, width = 5, height = 4, dpi = 150)

# 3d. Variogram of log(SOC) residuals
set.seed(42)
# OLS residuals using the same mean function as the GP model (intercept + elev + land use)
# Using these instead of a spatial polynomial keeps the variogram consistent with the model
lu_vg    <- model.matrix(~ lu - 1, data = dat)[, -1, drop = FALSE]
F_vg     <- cbind(1, dat$elev_sc, lu_vg)
bols     <- solve(crossprod(F_vg), crossprod(F_vg, dat$lsoc))
resid_vg <- dat$lsoc - as.vector(F_vg %*% bols)

# gd carries original lsoc — needed by MCMC in Section 4
gd <- as.geodata(dat, coords.col = c("x_km","y_km"), data.col = "lsoc")
gd <- jitterDupCoords(gd, max = 0.5)

# Variogram on covariate residuals using jittered coords from gd
gd_vg      <- gd
gd_vg$data <- resid_vg
vario <- variog(gd_vg, trend = "cte", max.dist = 0.5 * max(dist(gd$coords)), messages = FALSE)

png("figures/05_variogram.png", width = 600, height = 480, res = 130)
plot(vario, main = "Empirical variogram — log(SOC) OLS residuals (intercept + elev + land use)",
     xlab = "Distance (km)", ylab = "Semivariance", pch = 16, col = "steelblue4")
dev.off()

cat("EDA figures saved.\n")

# ── 4. MCMC  [HW Q2: Full Bayesian GP + MCMC; Q3: Priors + measurement error] ─
# Model: Y ~ N(Fβ, σ²(R_φ + δI))
# Q2: Collapsed block sampler — MH on log(φ), composition draws for σ² and β.
#     Covariate significance assessed from 95% credible intervals.
# Q3: Priors: β ∝ 1 (flat), σ² ∝ 1/σ² (Jeffreys), φ ∝ 1/φ (reference).
#     δ = relative nugget (τ²/σ²) fixed via pilot REML, incorporating the
#     measurement error SD reported in the RaCA data.

Y      <- gd$data
n      <- length(Y)
coords <- gd$coords
D      <- as.matrix(dist(coords))

lu_mat <- model.matrix(~ lu - 1, data = dat)[, -1, drop = FALSE]
# Spatial polynomial dropped: it competed with the GP covariance
# and collapsed φ→0 The GP now handles all spatial structure
F_mat  <- cbind(1, dat$elev_sc, lu_mat)
p <- ncol(F_mat)
beta_names <- c("intercept", "elev", colnames(lu_mat))

# Measurement error on log scale (delta method): var(log Y_obs - log Y_true) ≈ CV²
tau2_meas <- median((dat$soc_sd / dat$soc)^2, na.rm = TRUE)
cat(sprintf("Measurement error τ²_meas = %.5f  (median delta-method)\n", tau2_meas))

# Pilot REML to estimate relative nugget δ = τ²/σ² (regularises near-singular R)
cat("Fitting pilot REML for relative nugget δ...\n")
pilot_fit <- tryCatch(
  likfit(gd, trend = "cte", ini.cov.pars = c(0.25, 50),
         nugget = 0.2, fix.nugget = FALSE,
         cov.model = "matern", kappa = KAPPA,
         lik.method = "REML", messages = FALSE),
  error = function(e) NULL
)
if (!is.null(pilot_fit)) {
  tau2 <- pilot_fit$nugget / (pilot_fit$cov.pars[1] + pilot_fit$nugget)
  cat(sprintf("δ = %.3f  (REML: τ²=%.4f σ²=%.4f; measurement error = %.5f)\n",
              tau2, pilot_fit$nugget, pilot_fit$cov.pars[1], tau2_meas))
} else {
  tau2 <- 0.40   # fallback: moderate relative nugget typical for soil carbon
  cat(sprintf("δ = %.3f  (fallback — pilot REML failed)\n", tau2))
}
cat(sprintf("Fixed relative nugget δ = %.4f  (used in V_tilde = R + δI)\n", tau2))

# Combined function: returns marginal log-likelihood + GLS cache in one Cholesky pass.
# Caching this result avoids recomputing the O(n³) Cholesky on rejected MH steps.
gls_ll <- function(phi) {
  if (phi <= 0 || phi > max(D)) return(NULL)
  R <- matrix(matern(as.vector(D), phi = phi, kappa = KAPPA), n, n); diag(R) <- 1
  V <- R + tau2 * diag(n)
  L <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  Li_F <- backsolve(L, F_mat, transpose = TRUE)
  Li_Y <- backsolve(L, Y,     transpose = TRUE)
  LF   <- tryCatch(chol(crossprod(Li_F)), error = function(e) NULL)
  if (is.null(LF)) return(NULL)
  bh   <- backsolve(LF, backsolve(LF, crossprod(Li_F, Li_Y), transpose = TRUE))
  S2   <- sum((Li_Y - Li_F %*% bh)^2)
  if (S2 <= 0) return(NULL)
  ll   <- -0.5 * 2*sum(log(diag(L))) - 0.5 * 2*sum(log(diag(LF))) - (n-p)/2 * log(S2)
  list(ll=ll, V=V, L=L, LF=LF, beta_hat=bh, S2=S2)
}

# Lognormal prior on φ: log(φ) ~ N(MU_LOG_PHI, SD_LOG_PHI²)
# Parker & Sansó (2023, UCSC-SOE-23-13) fit land-use specific semivariograms to RaCA
# and report practical ranges of ~50–60 km across categories. We centre the prior at
# 60 km with a tight SD=0.3 on the log scale (95% prior mass ≈ 33–110 km).
# Sampling log(φ) with a Normal prior requires NO Jacobian in the MH ratio.
MU_LOG_PHI <- log(60)    # prior mean on log scale (~60 km, per Parker & Sansó 2023)
SD_LOG_PHI <- 0.3        # prior SD on log scale (tight; 95% mass ≈ 33–110 km)

if (!file.exists(MCMC_FILE)) {
  set.seed(2024)
  phi_chain  <- numeric(N_ITER - N_BURN)
  s2_chain   <- numeric(N_ITER - N_BURN)
  beta_chain <- matrix(NA, N_ITER - N_BURN, p)
  sig_prop   <- 0.4
  lp_cur     <- MU_LOG_PHI   # initialise at prior mean (~100 km)
  # Cache GLS quantities at current accepted φ — only recomputed on acceptance.
  # On rejected steps, the cached result is reused for the composition draw,
  # avoiding a redundant O(n³) Cholesky decomposition.
  cache_cur  <- gls_ll(exp(lp_cur))
  phi_acc <- 0L; win_acc <- 0L

  cat(sprintf("Running MCMC: %d burn-in + %d iterations...\n", N_BURN, N_ITER - N_BURN))
  pb <- txtProgressBar(0, N_ITER, style = 3)

  for (iter in seq_len(N_ITER)) {
    # Propose new φ and evaluate marginal likelihood (one Cholesky per iteration)
    lp_prop    <- lp_cur + rnorm(1, 0, sig_prop)
    cache_prop <- gls_ll(exp(lp_prop))
    ll_prop    <- if (is.null(cache_prop)) -Inf else cache_prop$ll

    log_ar <- (ll_prop - cache_cur$ll) +
              dnorm(lp_prop, MU_LOG_PHI, SD_LOG_PHI, log = TRUE) -
              dnorm(lp_cur,  MU_LOG_PHI, SD_LOG_PHI, log = TRUE)
    if (log(runif(1)) < log_ar) {
      lp_cur    <- lp_prop
      cache_cur <- cache_prop   # accepted: swap in the already-computed cache
      win_acc   <- win_acc + 1L
      if (iter > N_BURN) phi_acc <- phi_acc + 1L
    }
    # Rejected: cache_cur unchanged — no Cholesky recomputed for composition draw

    # Adaptive tuning during burn-in (every 200 steps)
    if (iter <= N_BURN && iter %% 200 == 0) {
      sig_prop <- sig_prop * exp(ifelse(win_acc / 200 > 0.4, 0.1, -0.1))
      sig_prop <- max(0.05, min(sig_prop, 2.0))
      win_acc  <- 0L
    }

    if (iter > N_BURN) {
      # Composition draws using cached Cholesky at current accepted φ
      # σ² | φ, Y ~ Inv-Gamma((n-p)/2, S²/2)
      s2  <- 1 / rgamma(1, shape = (n - p) / 2, rate = cache_cur$S2 / 2)
      # β | φ, σ², Y ~ N(β̂, σ²(F'V⁻¹F)⁻¹)
      C_b <- s2 * chol2inv(cache_cur$LF)
      bdr <- cache_cur$beta_hat + t(chol(C_b)) %*% rnorm(p)

      k <- iter - N_BURN
      phi_chain[k]   <- exp(lp_cur)
      s2_chain[k]    <- s2
      beta_chain[k,] <- bdr
    }
    setTxtProgressBar(pb, iter)
  }
  close(pb)

  # Thin
  idx <- seq(1, N_ITER - N_BURN, by = THIN)
  out <- list(phi = phi_chain[idx], sigma2 = s2_chain[idx],
              beta = beta_chain[idx,], tau2 = tau2, kappa = KAPPA)
  colnames(out$beta) <- beta_names
  saveRDS(out, MCMC_FILE)
  cat(sprintf("φ acceptance rate: %.1f%% | Posterior draws: %d\n",
              100 * phi_acc / (N_ITER - N_BURN), length(idx)))
}

mc  <- readRDS(MCMC_FILE)
np  <- length(mc$phi)

cat(sprintf("\nφ  posterior: mean=%.2f  95%%CI=[%.2f, %.2f] km\n",
            mean(mc$phi), quantile(mc$phi,.025), quantile(mc$phi,.975)))
cat(sprintf("σ² posterior: mean=%.4f  95%%CI=[%.4f, %.4f]\n",
            mean(mc$sigma2), quantile(mc$sigma2,.025), quantile(mc$sigma2,.975)))
cat("\nβ posterior (mean, 95% CI):\n")
bs <- data.frame(param=colnames(mc$beta),
                 mean=colMeans(mc$beta),
                 lo=apply(mc$beta,2,quantile,.025),
                 hi=apply(mc$beta,2,quantile,.975)) |>
  mutate(sig = !(lo < 0 & hi > 0))
print(bs, digits=3, row.names=FALSE)

# MCMC trace plots
png("figures/06_mcmc_traces.png", width = 960, height = 560, res = 130)
par(mfrow = c(2, 2), mar = c(3.5, 4, 2.5, 1))
plot(mc$phi,    type="l", col="steelblue4", lwd=.6, main="Trace: φ (km)",  ylab="φ")
plot(mc$sigma2, type="l", col="steelblue4", lwd=.6, main="Trace: σ²",      ylab="σ²")
plot(density(mc$phi),    main="Posterior: φ",  xlab="φ (km)",  col="steelblue4", lwd=2)
abline(v=quantile(mc$phi, c(.025,.975)), lty=2, col="firebrick")
plot(density(mc$sigma2), main="Posterior: σ²", xlab="σ²",      col="steelblue4", lwd=2)
abline(v=quantile(mc$sigma2, c(.025,.975)), lty=2, col="firebrick")
dev.off()

# β coefficient plot
png("figures/07_posterior_beta.png", width = 780, height = 540, res = 130)
par(mar = c(4, 8, 3, 2))
ord <- order(bs$mean)
plot(bs$mean[ord], seq_along(ord), pch=19, col="steelblue4",
     xlim=range(c(bs$lo, bs$hi)), yaxt="n", ylab="",
     xlab="Posterior mean ± 95% CI", main="Covariate effects on log(SOC)")
segments(bs$lo[ord], seq_along(ord), bs$hi[ord], seq_along(ord),
         col="steelblue4", lwd=2)
axis(2, at=seq_along(ord), labels=bs$param[ord], las=2, cex.axis=.85)
abline(v=0, lty=2, col="grey50")
dev.off()

cat("MCMC figures saved.\n")

# ── 5. Posterior predictive surface  [HW Q4: Predictive surface + uncertainty] ─
# Build grid, mask to state outlines via sf
lon_seq  <- seq(-80.6, -73.8, by = 0.12)   # PA west to NJ coast; ~10 km E-W (Parker & Sansó 2023)
lat_seq  <- seq(37.8,  42.6,  by = 0.09)   # MD south to PA north; ~10 km N-S
grid_ll  <- expand.grid(lon = lon_seq, lat = lat_seq)
state_sf <- st_as_sf(maps::map("state", regions = STATES_FULL,
                               fill = TRUE, plot = FALSE)) |>
  sf::st_set_crs(4326)
grid_sf  <- st_as_sf(grid_ll, coords = c("lon","lat"), crs = 4326)
in_st    <- as.logical(st_within(grid_sf, st_union(state_sf), sparse = FALSE))
grid_ll  <- grid_ll[in_st, ]
M        <- nrow(grid_ll)
cat(sprintf("Prediction grid: %d points\n", M))

ref_lon <- mean(dat$lon); ref_lat <- mean(dat$lat)
gkm <- data.frame(
  x_km = (grid_ll$lon - ref_lon) * cos(ref_lat * pi/180) * 111.32,
  y_km = (grid_ll$lat - ref_lat) * 111.32
)

# Precompute cross-distance matrix D_cross (M × n) — same for all draws
D_cross <- apply(gkm, 1, function(g)
  sqrt((coords[,1]-g[1])^2 + (coords[,2]-g[2])^2)) |> t()  # M × n

# Dominant land use for grid (plug-in Forest = most common)
dom_lu  <- names(which.max(table(dat$landuse)))
glu_mat <- matrix(0, M, ncol(lu_mat)); colnames(glu_mat) <- colnames(lu_mat)
if (dom_lu %in% colnames(glu_mat)) glu_mat[, dom_lu] <- 1

# IDW-interpolate elevation from observed sites to grid (power = 2)
# Gives spatially varying predictions even when φ is short.
elev_mn <- mean(dat$elevation_m, na.rm = TRUE)
elev_sd <- sd(dat$elevation_m,   na.rm = TRUE)
g_elev_sc <- apply(D_cross, 1, function(d) {
  w <- 1 / pmax(d^2, 1e-4)
  sum(w * dat$elevation_m, na.rm = TRUE) / sum(w)
})
g_elev_sc <- (g_elev_sc - elev_mn) / elev_sd

F_g <- cbind(1, g_elev_sc, glu_mat)   # intercept, IDW elev, land use

n_draw  <- min(np, 300)
draw_idx <- round(seq(1, np, length.out = n_draw))
mu_mat  <- matrix(NA, n_draw, M)
var_mat <- matrix(NA, n_draw, M)

cat(sprintf("Posterior predictive: %d draws × %d grid points...\n", n_draw, M))
pb <- txtProgressBar(0, n_draw, style = 3)

for (k in seq_len(n_draw)) {
  ii     <- draw_idx[k]
  phi    <- mc$phi[ii]; s2 <- mc$sigma2[ii]; beta <- mc$beta[ii,]
  R <- matrix(matern(as.vector(D), phi=phi, kappa=KAPPA), n, n); diag(R) <- 1
  # Under Y ~ N(Fβ, σ²(R+δI)), V_tilde = R + δI does not depend on σ².
  # Kriging mean: F₀'β + R_cross·V_tilde⁻¹·(Y−Fβ)   [σ² cancels]
  # Kriging var:  σ²·(1 + δ − diag(R_cross·V_tilde⁻¹·R_cross'))
  V <- R + tau2 * diag(n)
  L <- tryCatch(chol(V), error=function(e) NULL)
  if (is.null(L)) next

  Yr   <- Y - F_mat %*% beta
  Vi_r <- backsolve(L, backsolve(L, Yr, transpose=TRUE))   # V_tilde⁻¹(Y−Fβ)

  R_cross <- matrix(matern(as.vector(D_cross), phi=phi, kappa=KAPPA), M, n)

  mu_mat[k,]  <- as.vector(F_g %*% beta) + R_cross %*% Vi_r   # σ² cancels in mean

  tmp          <- backsolve(L, t(R_cross), transpose = TRUE)   # n×M
  var_mat[k,] <- pmax(s2 * (1 + tau2 - colSums(tmp^2)), 0)    # σ²·(1+δ−r'V⁻¹r)

  setTxtProgressBar(pb, k)
}
close(pb)

pred_mean <- colMeans(mu_mat, na.rm=TRUE)
pred_var  <- colMeans(var_mat, na.rm=TRUE) + apply(mu_mat, 2, var, na.rm=TRUE)
pred_sd   <- sqrt(pred_var)

pred_df <- data.frame(lon=grid_ll$lon, lat=grid_ll$lat,
                      log_soc=pred_mean, pred_sd=pred_sd)

p_pmean <- ggplot() +
  geom_tile(data=pred_df, aes(lon, lat, fill=log_soc)) +
  geom_polygon(data=border, aes(long, lat, group=group),
               fill=NA, colour="grey35", linewidth=0.35) +
  geom_point(data=dat, aes(lon,lat), colour="white", size=0.8, alpha=0.5) +
  scale_fill_viridis_c(name="log(SOC)") + coord_fixed(1.3) +
  labs(title="Posterior predictive mean — log(SOC)",
       subtitle="MD / NJ / DE / PA   Bayesian GP, Matérn κ=1.5",
       x="Longitude", y="Latitude") + theme_bw(11)
ggsave("figures/08_pred_mean.png", p_pmean, width=8, height=6, dpi=150)

p_psd <- ggplot() +
  geom_tile(data=pred_df, aes(lon, lat, fill=pred_sd)) +
  geom_polygon(data=border, aes(long, lat, group=group),
               fill=NA, colour="grey35", linewidth=0.35) +
  geom_point(data=dat, aes(lon,lat), colour="white", size=0.8, alpha=0.5) +
  scale_fill_distiller(name="SD[log(SOC)]", palette="YlOrRd", direction=1) +
  coord_fixed(1.3) +
  labs(title="Posterior predictive SD — log(SOC)",
       subtitle="White dots = observation sites", x="Longitude", y="Latitude") +
  theme_bw(11)
ggsave("figures/09_pred_sd.png", p_psd, width=8, height=6, dpi=150)

cat("Prediction maps saved.\n")

# ── 6. MAP estimate  [HW Q5: Compare MAP to MCMC] ────────────────────────────
phi_grid <- seq(0.5, max(D)*0.9, length.out = 300)
ll_grid  <- sapply(phi_grid, function(p) {
  r <- gls_ll(p)
  if (is.null(r)) -Inf else r$ll + dnorm(log(p), MU_LOG_PHI, SD_LOG_PHI, log = TRUE)
})
ll_norm  <- ll_grid - max(ll_grid, na.rm=TRUE)
phi_map  <- phi_grid[which.max(ll_grid)]

# REML σ² at φ_MAP
gl_map     <- gls_ll(phi_map)
sigma2_map <- gl_map$S2 / (n - p)
beta_map   <- as.vector(gl_map$beta_hat)

cat(sprintf("\nMAP:  φ_MAP = %.2f km   σ²_MAP = %.4f\n", phi_map, sigma2_map))
ci_phi <- range(phi_grid[!is.na(ll_norm) & ll_norm > -1.92])
cat(sprintf("95%% CI (grid): [%.2f, %.2f] km\n", ci_phi[1], ci_phi[2]))

# Figure 10a: MAP — marginal log-likelihood profile for φ
png("figures/10a_map_phi.png", width = 700, height = 520, res = 130)
par(mar = c(4.5, 4.5, 3.5, 1.5))
ylim_phi <- c(max(-10, min(ll_norm, na.rm=TRUE)), 0.3)
plot(phi_grid, ll_norm, type="l", lwd=2.5, col="steelblue4",
     main=sprintf("MAP — log-posterior profile (likelihood + Lognormal prior)\nφ_MAP = %.2f km", phi_map),
     xlab="φ (km)", ylab="Normalised log-posterior", ylim=ylim_phi)
abline(h=-1.92, lty=2, col="grey60", lwd=1.2)
abline(v=phi_map, col="firebrick", lty=2, lwd=1.8)
legend("topright", bty="n", lwd=c(2.5,1.8,1.2), lty=c(1,2,2),
       col=c("steelblue4","firebrick","grey60"),
       legend=c("Log-likelihood (REML)", sprintf("φ_MAP = %.2f km", phi_map),
                "−1.92 threshold (95% CI"), cex=0.82)
dev.off()

# Figure 10b: MCMC posterior densities for φ and σ²
png("figures/10b_mcmc_posteriors.png", width = 960, height = 480, res = 130)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

d_phi <- density(mc$phi, bw="SJ")
plot(d_phi, lwd=2.5, col="steelblue4",
     main=sprintf("MCMC posterior: φ\nmean=%.2f km  95%%CI=[%.2f, %.2f]",
                  mean(mc$phi), quantile(mc$phi,0.025), quantile(mc$phi,0.975)),
     xlab="φ (km)", ylab="Density")
abline(v=phi_map, col="firebrick", lty=2, lwd=1.8)
abline(v=mean(mc$phi), col="darkorange", lty=1, lwd=1.5)
legend("topright", bty="n", lwd=c(2.5,1.8,1.5), lty=c(1,2,1),
       col=c("steelblue4","firebrick","darkorange"),
       legend=c("MCMC posterior","φ_MAP","Posterior mean"), cex=0.82)

d_s2 <- density(mc$sigma2, bw="SJ")
plot(d_s2, lwd=2.5, col="steelblue4",
     main=sprintf("MCMC posterior: σ²\nmean=%.3f  σ²_MAP=%.3f", mean(mc$sigma2), sigma2_map),
     xlab="σ²", ylab="Density")
abline(v=sigma2_map, col="firebrick", lty=2, lwd=1.8)
abline(v=mean(mc$sigma2), col="darkorange", lty=1, lwd=1.5)
legend("topright", bty="n", lwd=c(2.5,1.8,1.5), lty=c(1,2,1),
       col=c("steelblue4","firebrick","darkorange"),
       legend=c("MCMC posterior","σ²_MAP","Posterior mean"), cex=0.82)
dev.off()

cat("MAP comparison saved (10a + 10b).\n")

# ── 7. Goodness of fit  [HW Q6: Convince reader model fits reasonably] ────────
# LOO-CV using matrix inversion lemma (exact, no re-fitting)
phi_pm <- mean(mc$phi); s2_pm <- mean(mc$sigma2); beta_pm <- colMeans(mc$beta)

R_pm   <- matrix(matern(as.vector(D), phi=phi_pm, kappa=KAPPA), n, n); diag(R_pm) <- 1
V_pm   <- s2_pm * (R_pm + tau2 * diag(n))
L_pm   <- chol(V_pm)
Vi_pm  <- chol2inv(L_pm)                         # V⁻¹  (n×n)
mu_pm  <- as.vector(F_mat %*% beta_pm)
Vi_r   <- as.vector(Vi_pm %*% (Y - mu_pm))       # V⁻¹(Y − Fβ)
vi_diag <- diag(Vi_pm)                            # diagonal of V⁻¹

loo_var  <- 1 / vi_diag
loo_mean <- Y - Vi_r / vi_diag
resid_l  <- Y - loo_mean
rmse_loo <- sqrt(mean(resid_l^2))
bias_loo <- mean(resid_l)
cov95    <- mean(abs(resid_l) < 1.96 * sqrt(loo_var))

# PIT
pit <- pnorm(Y, loo_mean, sqrt(loo_var))
ks  <- ks.test(pit, "punif")

# CRPS (LEC11): for Gaussian N(μ,σ²), closed-form is
#   crps = σ[1/√π − 2φ(z) − z(2Φ(z)−1)]  where z = (y−μ)/σ
# Using plug-in posterior mean parameters for the LOO predictive distribution.
crps_gauss <- function(y, mu, sigma) {
  z <- (y - mu) / sigma
  sigma * (z * (2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi))
}
crps_vals <- crps_gauss(Y, loo_mean, sqrt(loo_var))
mean_crps <- mean(crps_vals)

# Posterior predictive replicates
set.seed(99)
Yrep <- replicate(200, as.vector(mu_pm + t(L_pm) %*% rnorm(n)))

cat(sprintf("\nGoodness of fit:\n"))
cat(sprintf("  LOO-RMSE     = %.4f\n", rmse_loo))
cat(sprintf("  LOO-Bias     = %.5f\n", bias_loo))
cat(sprintf("  95%% coverage = %.1f%%  (nominal 95%%)\n", 100*cov95))
cat(sprintf("  PIT KS p     = %.4f  (%s)\n", ks$p.value,
            ifelse(ks$p.value > 0.05, "well-calibrated", "some miscalibration")))
cat(sprintf("  Mean CRPS    = %.4f  (lower = better; in log-SOC units)\n", mean_crps))

png("figures/11_gof.png", width = 1280, height = 340, res = 130)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))

# LOO scatter
lim <- range(c(Y, loo_mean))
plot(Y, loo_mean, pch=19, cex=0.75, col="steelblue4", xlim=lim, ylim=lim,
     main="LOO cross-validation", xlab="Observed log(SOC)", ylab="LOO predicted")
abline(0, 1, col="firebrick", lwd=1.5, lty=2)
segments(Y, loo_mean - 1.96*sqrt(loo_var), Y, loo_mean + 1.96*sqrt(loo_var),
         col=adjustcolor("steelblue4", 0.3))
legend("topleft", bty="n", cex=0.75,
       legend=c(sprintf("RMSE=%.3f", rmse_loo), sprintf("Cov95=%.0f%%", 100*cov95)))

# PIT
hist(pit, breaks=15, col="steelblue3", border="white", freq=FALSE,
     main=sprintf("PIT  (KS p=%.3f)", ks$p.value), xlab="PIT value")
abline(h=1, col="firebrick", lty=2, lwd=1.5)

# Posterior predictive check
dens_obs <- density(Y, bw="SJ")
dens_rep <- lapply(seq_len(ncol(Yrep)), function(k) density(Yrep[,k], bw="SJ"))
ylim_ppc <- c(0, max(c(dens_obs$y, sapply(dens_rep, \(d) max(d$y)))))
plot(NA, xlim=range(c(Y, Yrep)), ylim=ylim_ppc,
     main="Posterior predictive check", xlab="log(SOC)", ylab="Density")
for (d in dens_rep) lines(d, col=adjustcolor("steelblue3",.12), lwd=0.8)
lines(dens_obs, col="firebrick", lwd=2.5)
legend("topright", bty="n", lwd=c(0.8,2.5),
       col=c("steelblue3","firebrick"),
       legend=c("Predictive replicates","Observed"), cex=0.8)

# CRPS spatial map
plot(dat$x_km, dat$y_km, pch=19, cex=0.7,
     col=colorRampPalette(c("steelblue4","gold","firebrick"))(100)[
       pmin(100, pmax(1, round((crps_vals - min(crps_vals)) /
                               (max(crps_vals) - min(crps_vals)) * 99) + 1))],
     main=sprintf("LOO CRPS (mean=%.3f)", mean_crps),
     xlab="x (km)", ylab="y (km)")
dev.off()

cat("GOF figures saved.\n")
cat("\nAll done. Figures in hw4-raca/figures/\n")
