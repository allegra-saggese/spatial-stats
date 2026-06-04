# 06_profile_likelihood.R
# Profile likelihood for sigma^2 (sill) and phi (range) in the Matern family
# kappa = 0.5, 1.0, 1.5, 2.5
#
# Method:
#   1. Fit likfit() to get MLE of trend parameters beta and nugget tau^2
#   2. Fix beta at MLE; fix tau^2 at MLE (or WLS estimate for degenerate cases)
#   3. Evaluate log-likelihood on a grid of (sigma^2, phi) via Cholesky
#   4. Profile out each nuisance parameter by taking the max over the other

library(readr)
library(dplyr)
library(geoR)

va   <- read_csv("VA_soil_carbon_clean.csv", show_col_types = FALSE)
va_m <- filter(va, soc_measured == "measured") |> mutate(lsoc = log(soc))

ref_lon <- mean(va_m$lon); ref_lat <- mean(va_m$lat)
va_m <- va_m |>
  mutate(x_km = (lon - ref_lon) * cos(ref_lat * pi / 180) * 111.32,
         y_km = (lat - ref_lat) * 111.32)

set.seed(42)
gd <- as.geodata(va_m, coords.col = c("x_km", "y_km"), data.col = "lsoc")
gd <- jitterDupCoords(gd, max = 0.5)

kappas <- c(0.5, 1.0, 1.5, 2.5)
cols   <- c("#2166ac", "#4dac26", "#d6604d", "#8B1A1A")
ltys   <- 1:4

# ── 1. Fit ML models to get beta and nugget MLEs ──────────────────────────────
cat("Fitting ML models for trend parameters...\n")

fits_ml <- lapply(kappas, function(k) {
  fix_nug <- k >= 1.5
  tryCatch(
    likfit(gd, trend = "2nd",
           ini.cov.pars = c(0.20, 40),
           nugget       = if (fix_nug) 0.212 else 0.03,
           fix.nugget   = fix_nug,
           cov.model    = "matern", kappa = k,
           lik.method   = "ML", messages = FALSE),
    error = function(e) NULL
  )
})
names(fits_ml) <- paste0("k", kappas)

for (i in seq_along(kappas)) {
  f <- fits_ml[[i]]
  if (!is.null(f))
    cat(sprintf("kappa=%.1f  sigma2=%.4f  phi=%.2f  nugget=%.4f  logLik=%.2f\n",
                kappas[i], f$cov.pars[1], f$cov.pars[2], f$nugget, f$loglik))
}

# ── 2. Precompute distances and design matrix (shared across all kappas) ───────
Y    <- gd$data
n    <- length(Y)
D    <- as.matrix(dist(gd$coords))   # n x n distance matrix

# Design matrix for 2nd order trend (intercept + x + y + x^2 + y^2 + xy)
coords <- gd$coords
F_mat <- cbind(1,
               coords[,1], coords[,2],
               coords[,1]^2, coords[,2]^2,
               coords[,1] * coords[,2])

# ── 3. Log-likelihood function (profiling out beta via GLS) ───────────────────
log_lik_grid <- function(sigmasq, phi, nugget, kappa, D, Y, F_mat) {
  n <- nrow(D)
  # Correlation matrix via Matern
  R     <- matrix(matern(as.vector(D), phi = phi, kappa = kappa), n, n)
  diag(R) <- 1
  Sigma <- sigmasq * R + nugget * diag(n)

  # Cholesky decomposition
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(L)) return(-Inf)

  log_det <- 2 * sum(log(diag(L)))

  # GLS estimate of beta (profiles out beta analytically)
  Li_F <- backsolve(L, F_mat, transpose = TRUE)   # L^{-T} F
  Li_Y <- backsolve(L, Y,     transpose = TRUE)   # L^{-T} Y
  beta_hat <- solve(crossprod(Li_F), crossprod(Li_F, Li_Y))

  resid  <- Li_Y - Li_F %*% beta_hat
  quad   <- sum(resid^2)

  -n/2 * log(2 * pi) - 0.5 * log_det - 0.5 * quad
}

# ── 4. Evaluate likelihood over grid ─────────────────────────────────────────
sill_grid  <- seq(0.02, 0.80, length.out = 60)
range_grid <- seq(2,    180,  length.out = 60)

cat("\nComputing likelihood grids (60x60 per kappa)...\n")

ll_grids <- lapply(seq_along(kappas), function(i) {
  k   <- kappas[i]
  f   <- fits_ml[[i]]
  nug <- if (!is.null(f)) f$nugget else if (k >= 1.5) 0.212 else 0.03

  cat(sprintf("  kappa = %.1f  (nugget plugged in = %.4f)\n", k, nug))

  mat <- outer(sill_grid, range_grid, FUN = Vectorize(function(s, p)
    log_lik_grid(s, p, nug, k, D, Y, F_mat)
  ))
  mat
})
names(ll_grids) <- paste0("k", kappas)

# ── 5. Extract 1D profile likelihoods ─────────────────────────────────────────
# pl(sigma^2) = max over phi  [row maximum of likelihood grid]
# pl(phi)     = max over sigma^2  [column maximum]

prof_sill  <- lapply(ll_grids, function(m) apply(m, 1, max))
prof_range <- lapply(ll_grids, function(m) apply(m, 2, max))

# Normalise so peak = 0
prof_sill  <- lapply(prof_sill,  function(v) v - max(v, na.rm = TRUE))
prof_range <- lapply(prof_range, function(v) v - max(v, na.rm = TRUE))

# ── 6. Plot profile likelihood for sill ───────────────────────────────────────
png("figures/14_profile_sill.png", width = 820, height = 560, res = 120)

ylim_s <- range(unlist(prof_sill), na.rm = TRUE)
ylim_s[1] <- max(ylim_s[1], -10)

plot(NA, xlim = range(sill_grid), ylim = ylim_s,
     main = expression("Profile log-likelihood — sill " * sigma^2),
     xlab = expression(sigma^2 ~ "(partial sill)"),
     ylab = "Profile log-likelihood (relative to maximum)")
abline(h = -1.92, lty = 3, col = "grey50")
text(0.78, -1.92 + 0.3, "95% CI", col = "grey50", adj = 1, cex = 0.8)

for (i in seq_along(kappas))
  lines(sill_grid, prof_sill[[i]], col = cols[i], lwd = 2, lty = ltys[i])

legend("bottomright", bty = "n", lwd = 2, lty = ltys, col = cols,
       legend = sprintf("kappa = %.1f", kappas))
dev.off()
cat("Saved: figures/14_profile_sill.png\n")

# ── 7. Plot profile likelihood for range ──────────────────────────────────────
png("figures/15_profile_range.png", width = 820, height = 560, res = 120)

ylim_r <- range(unlist(prof_range), na.rm = TRUE)
ylim_r[1] <- max(ylim_r[1], -10)

plot(NA, xlim = range(range_grid), ylim = ylim_r,
     main = "Profile log-likelihood — range φ",
     xlab = "Range φ (km)",
     ylab = "Profile log-likelihood (relative to maximum)")
abline(h = -1.92, lty = 3, col = "grey50")
text(175, -1.92 + 0.3, "95% CI", col = "grey50", adj = 1, cex = 0.8)

for (i in seq_along(kappas))
  lines(range_grid, prof_range[[i]], col = cols[i], lwd = 2, lty = ltys[i])

legend("bottomright", bty = "n", lwd = 2, lty = ltys, col = cols,
       legend = sprintf("kappa = %.1f", kappas))
dev.off()
cat("Saved: figures/15_profile_range.png\n")

# ── 8. Report MLEs and 95% CIs from profile likelihoods ──────────────────────
cat("\n── Profile likelihood summary ────────────────────────────────\n")
for (i in seq_along(kappas)) {
  ci_s <- range(sill_grid[prof_sill[[i]]   > -1.92], na.rm = TRUE)
  ci_r <- range(range_grid[prof_range[[i]] > -1.92], na.rm = TRUE)
  mle_s <- sill_grid[which.max(prof_sill[[i]])]
  mle_r <- range_grid[which.max(prof_range[[i]])]
  cat(sprintf("kappa=%.1f  sigma2: MLE=%.3f  95%%CI=[%.3f,%.3f]  |  phi: MLE=%.1f km  95%%CI=[%.1f,%.1f]\n",
              kappas[i], mle_s, ci_s[1], ci_s[2], mle_r, ci_r[1], ci_r[2]))
}
