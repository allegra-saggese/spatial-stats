# 04_variogram_fit.R
# Fit theoretical variogram models to log(SOC) — measured sites only
# Uses geoR::likfit() (ML on full data) rather than WLS on binned empirical variogram
# Compares exponential / spherical / Matern, isotropic and anisotropic

library(readr)
library(dplyr)
library(geoR)

va    <- read_csv("VA_soil_carbon_clean.csv", show_col_types = FALSE)
va_m  <- filter(va, soc_measured == "measured") |> mutate(lsoc = log(soc))

ref_lon <- mean(va_m$lon); ref_lat <- mean(va_m$lat)
va_m <- va_m |>
  mutate(x_km = (lon - ref_lon) * cos(ref_lat * pi / 180) * 111.32,
         y_km = (lat - ref_lat) * 111.32)

set.seed(42)
gd <- as.geodata(va_m, coords.col = c("x_km", "y_km"), data.col = "lsoc")
gd <- jitterDupCoords(gd, max = 0.5)

vario_omni <- variog(gd, trend = "2nd", max.dist = 291, messages = FALSE)

# Starting values from empirical variogram:
#   nugget  ≈ 0.05 (first bin)
#   partial sill (sigmasq) ≈ 0.25
#   range (phi) ≈ 80 km
ini <- c(sigmasq = 0.25, phi = 80)
nug <- 0.05

# ── Fit isotropic models ──────────────────────────────────────────────────────
cat("Fitting isotropic models...\n")

fit_exp <- likfit(gd, trend = "2nd", ini.cov.pars = ini, nugget = nug,
                  cov.model = "exponential", lik.method = "ML", messages = FALSE)

fit_sph <- likfit(gd, trend = "2nd", ini.cov.pars = ini, nugget = nug,
                  cov.model = "spherical",   lik.method = "ML", messages = FALSE)

fit_mat <- likfit(gd, trend = "2nd", ini.cov.pars = ini, nugget = nug,
                  cov.model = "matern", kappa = 1.5,
                  lik.method = "ML", messages = FALSE)

# ── Model comparison table ────────────────────────────────────────────────────
summarise_fit <- function(fit, label) {
  data.frame(
    model        = label,
    nugget       = round(fit$nugget, 4),
    partial_sill = round(fit$sigmasq, 4),
    range_km     = round(fit$phi, 2),
    total_sill   = round(fit$nugget + fit$sigmasq, 4),
    logLik       = round(fit$loglik, 2),
    AIC          = round(-2 * fit$loglik + 2 * fit$npars, 2)
  )
}

tbl <- rbind(
  summarise_fit(fit_exp, "Exponential"),
  summarise_fit(fit_sph, "Spherical"),
  summarise_fit(fit_mat, "Matern k=1.5")
)

cat("\n── Model comparison (ML) ─────────────────────────────────────\n")
print(tbl, row.names = FALSE)

# nugget ratio = nugget / total sill
best <- tbl[which.min(tbl$AIC), ]
cat(sprintf("\nBest model: %s  (AIC = %.1f)\n", best$model, best$AIC))
cat(sprintf("  Nugget ratio (tau2/sill): %.3f\n",
            best$nugget / best$total_sill))

# ── Plot: empirical variogram + all three model fits ─────────────────────────
png("figures/11_variogram_fitted.png", width = 820, height = 560, res = 120)

plot(vario_omni,
     main = "Empirical variogram + fitted models — log(SOC)",
     xlab = "Distance (km)", ylab = "Semivariance",
     pch = 16, col = "grey40")

lines(fit_exp, max.dist = 291, col = "steelblue", lwd = 2, lty = 1)
lines(fit_sph, max.dist = 291, col = "firebrick", lwd = 2, lty = 2)
lines(fit_mat, max.dist = 291, col = "darkgreen", lwd = 2, lty = 3)

legend("topright", bty = "n", lwd = 2, lty = 1:3,
       col = c("steelblue", "firebrick", "darkgreen"),
       legend = c(
         sprintf("Exponential   AIC = %.1f", tbl$AIC[1]),
         sprintf("Spherical     AIC = %.1f", tbl$AIC[2]),
         sprintf("Matern k=1.5  AIC = %.1f", tbl$AIC[3])
       ))

dev.off()
cat("Saved: figures/11_variogram_fitted.png\n")
