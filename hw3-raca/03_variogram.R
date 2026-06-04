# 03_variogram.R
# Residual plots and empirical variogram (omnidirectional + directional)
# Data: measured SOC only, log-transformed, 2nd order spatial trend via geoR

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(maps)
library(geoR)

va      <- read_csv("VA_soil_carbon_clean.csv", show_col_types = FALSE)
va_m    <- filter(va, soc_measured == "measured") |> mutate(lsoc = log(soc))
va_state <- map_data("state") |> filter(region == "virginia")

dir.create("figures", showWarnings = FALSE)

# ── Convert lon/lat to km ──────────────────────────────────────────────────────
# geoR uses Euclidean distance on raw coords — degrees are not isotropic
# (1° lon ≠ 1° lat in km at mid-latitudes). Flat-earth approximation centred on VA.
ref_lon <- mean(va_m$lon)
ref_lat <- mean(va_m$lat)

va_m <- va_m |>
  mutate(
    x_km = (lon - ref_lon) * cos(ref_lat * pi / 180) * 111.32,
    y_km = (lat - ref_lat) * 111.32
  )

# ── Build geoR geodata object ─────────────────────────────────────────────────
# RaCA coordinates are truncated to 2 decimal places (~1 km) by USDA,
# producing many co-located points. Jitter by up to 0.5 km so variog4()
# can compute directional angles without NA errors.
set.seed(42)
gd <- as.geodata(va_m,
                 coords.col = c("x_km", "y_km"),
                 data.col   = "lsoc")
gd <- jitterDupCoords(gd, max = 0.5)

# ── 1. Residuals from 2nd order trend ────────────────────────────────────────
# Fit the same 2nd order surface to extract residuals for plotting.
# geoR uses trend="2nd" internally; we replicate with lm() here for ggplot.
fit2 <- lm(lsoc ~ lon + lat + I(lon^2) + I(lat^2) + lon:lat, data = va_m)
va_m$resid <- residuals(fit2)

cat(sprintf("Residual summary — mean: %.4f  sd: %.3f  skew: %.3f\n",
            mean(va_m$resid), sd(va_m$resid),
            (mean(va_m$resid) - median(va_m$resid)) / sd(va_m$resid)))

pal_div <- scale_color_distiller(palette = "RdBu", name = "Residual",
                                 limits = \(x) c(-max(abs(x)), max(abs(x))))

p_rmap <- ggplot() +
  geom_polygon(data = va_state, aes(long, lat, group = group),
               fill = "grey92", color = "grey60", linewidth = 0.3) +
  geom_point(data = va_m,
             aes(lon, lat, color = resid, size = abs(resid)),
             alpha = 0.85) +
  pal_div +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  coord_fixed(1.3) +
  labs(title = "Residuals from 2nd order trend (measured SOC)",
       x = "Longitude", y = "Latitude")

p_rhist <- ggplot(va_m, aes(resid)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Residual distribution", x = "Residual (log SOC)", y = "Count")

p_rqq <- ggplot(va_m, aes(sample = resid)) +
  stat_qq(alpha = 0.6) + stat_qq_line(color = "red") +
  labs(title = "Residual Q-Q plot", x = "Theoretical", y = "Sample")

ggsave("figures/08_residuals.png",
       p_rmap / (p_rhist | p_rqq),
       width = 9, height = 9, dpi = 150)
cat("Saved: figures/08_residuals.png\n")

# ── 2. Omnidirectional empirical variogram ────────────────────────────────────
# max.dist: rule of thumb = ~half the max pairwise distance
max_d <- max(dist(cbind(va_m$x_km, va_m$y_km))) / 2
cat(sprintf("Max pairwise distance: %.0f km  →  cutoff: %.0f km\n",
            max_d * 2, max_d))

vario_omni <- variog(gd,
                     trend    = "2nd",
                     max.dist = max_d,
                     messages = FALSE)

png("figures/09_variogram_omni.png", width = 700, height = 500, res = 120)
plot(vario_omni,
     main = "Omnidirectional variogram — log(SOC), 2nd order trend removed",
     xlab = "Distance (km)", ylab = "Semivariance",
     pch = 16, col = "steelblue")
dev.off()
cat("Saved: figures/09_variogram_omni.png\n")

# ── 3. Directional variogram (anisotropy check) ───────────────────────────────
# variog4() computes variograms in 4 directions: 0°, 45°, 90°, 135°
# If the curves differ in sill/range → geometric anisotropy
vario4 <- variog4(gd,
                  trend    = "2nd",
                  max.dist = max_d,
                  messages = FALSE)

png("figures/10_variogram_directional.png", width = 800, height = 600, res = 120)
plot(vario4,
     main = "Directional variograms — log(SOC), 2nd order trend removed",
     xlab = "Distance (km)", ylab = "Semivariance",
     legend = TRUE)
dev.off()
cat("Saved: figures/10_variogram_directional.png\n")
