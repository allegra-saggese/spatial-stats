# 05_matern_fit.R
# Fit Matérn family covariograms with smoothness kappa = 0.5, 1, 1.5, 2.5
# using weighted least squares (variofit). Compare LSE values and plot.
#
# Note: fitting γ(h) [variogram] is equivalent to fitting C(h) [covariogram]
# since C(h) = C(0) - γ(h); the LS estimates are identical in both spaces.
# We plot both representations.

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

vario <- variog(gd, trend = "2nd", max.dist = 291, messages = FALSE)

# ── Fit Matérn models with kappa = 0.5, 1, 1.5, 2.5 ─────────────────────────
# kappa = 0.5 → exponential (roughest, 0 mean-square derivatives)
# kappa = 1.5 → 1 derivative (moderately smooth)
# kappa = 2.5 → 2 derivatives (smoother)
# kappa → inf → Gaussian (infinitely smooth)

kappas <- c(0.5, 1.0, 1.5, 2.5)
ini    <- c(sigmasq = 0.25, phi = 80)

fits <- lapply(kappas, function(k) {
  variofit(vario,
           ini.cov.pars = ini,
           nugget       = 0.05,
           cov.model    = "matern",
           kappa        = k,
           weights      = "npairs",
           messages     = FALSE)
})
names(fits) <- paste0("kappa_", kappas)

# ── Results table ─────────────────────────────────────────────────────────────
tbl <- do.call(rbind, lapply(seq_along(kappas), function(i) {
  f <- fits[[i]]
  sigmasq <- f$cov.pars[1]   # partial sill
  phi     <- f$cov.pars[2]   # range parameter
  data.frame(
    kappa        = kappas[i],
    nugget       = round(f$nugget,          4),
    partial_sill = round(sigmasq,           4),
    range_km     = round(phi,               2),
    total_sill   = round(f$nugget + sigmasq, 4),
    LSE          = round(f$value,           6)
  )
}))

cat("\n── Matérn family — WLS fit comparison ───────────────────────\n")
print(tbl, row.names = FALSE)

best_i <- which.min(tbl$LSE)
cat(sprintf("\nBest fit: kappa = %.1f  (LSE = %.6f)\n",
            tbl$kappa[best_i], tbl$LSE[best_i]))

# ── Plot 1: variogram space ───────────────────────────────────────────────────
cols <- c("#2166ac", "#4dac26", "#d6604d", "#8B1A1A")
ltys <- 1:4

png("figures/12_matern_variogram.png", width = 820, height = 560, res = 120)
plot(vario,
     main = "Matérn family fits — variogram of log(SOC)",
     xlab = "Distance (km)", ylab = "Semivariance",
     pch = 16, col = "grey40", ylim = c(0, 0.55))

for (i in seq_along(fits)) {
  lines(fits[[i]], max.dist = 291,
        col = cols[i], lwd = 2, lty = ltys[i])
}

legend("topright", bty = "n", lwd = 2, lty = ltys, col = cols,
       legend = sprintf("κ = %.1f   LSE = %.4f", tbl$kappa, tbl$LSE))
dev.off()
cat("Saved: figures/12_matern_variogram.png\n")

# ── Plot 2: covariogram space  C(h) = C(0) − γ(h) ───────────────────────────
# C(0) estimated as the total sill of the best-fitting model
C0 <- tbl$total_sill[best_i]

cov_emp <- C0 - vario$v          # empirical covariogram
d_seq   <- seq(0, 291, length.out = 300)

matern_cov <- function(h, sigmasq, phi, kappa, nugget) {
  # Matérn covariance: C(h) = sigmasq * Matérn_corr(h) for h > 0, C(0) = sigmasq + nugget
  ifelse(h == 0, sigmasq + nugget,
         sigmasq * geoR::matern(h, phi = phi, kappa = kappa))
}

png("figures/13_matern_covariogram.png", width = 820, height = 560, res = 120)
plot(vario$u, cov_emp,
     main = "Matérn family fits — covariogram of log(SOC)",
     xlab = "Distance (km)", ylab = "Covariance C(h)",
     pch = 16, col = "grey40",
     ylim = c(min(0, min(cov_emp)) - 0.02, C0 + 0.05))
abline(h = 0, lty = 3, col = "grey60")

for (i in seq_along(fits)) {
  f   <- fits[[i]]
  C_h <- matern_cov(d_seq, f$cov.pars[1], f$cov.pars[2], kappas[i], f$nugget)
  lines(d_seq, C_h, col = cols[i], lwd = 2, lty = ltys[i])
}

legend("topright", bty = "n", lwd = 2, lty = ltys, col = cols,
       legend = sprintf("κ = %.1f   LSE = %.4f", tbl$kappa, tbl$LSE))
dev.off()
cat("Saved: figures/13_matern_covariogram.png\n")
