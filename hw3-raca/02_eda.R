# 02_eda.R
# Graphical exploration of VA soil carbon data
# Addresses: spatial trend, normality, measured vs. modeled, elevation covariate

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(maps)

va <- read_csv("VA_soil_carbon_clean.csv", show_col_types = FALSE)
va_state <- map_data("state") |> filter(region == "virginia")

dir.create("figures", showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))
pal <- c(measured = "#2166ac", modeled = "#d6604d")

# ── 1. Spatial distribution map ───────────────────────────────────────────────
# First look: where are the samples, colored by SOC value
p_map <- ggplot() +
  geom_polygon(data = va_state, aes(long, lat, group = group),
               fill = "grey92", color = "grey60", linewidth = 0.3) +
  geom_point(data = va, aes(lon, lat, color = soc, shape = soc_measured),
             size = 1.8, alpha = 0.8) +
  scale_color_viridis_c(name = "SOC (%)", option = "plasma", trans = "log") +
  scale_shape_manual(values = c(measured = 16, modeled = 17),
                     name = "Type") +
  coord_fixed(1.3) +
  labs(title = "Soil organic carbon across Virginia",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right")

ggsave("figures/01_map_soc.png", p_map, width = 8, height = 4, dpi = 150)

# ── 2. First / second order spatial trend ────────────────────────────────────
# SOC vs. longitude and latitude with linear + quadratic smooths
# Evidence of a trend: if the loess departs systematically from flat

p_lon <- ggplot(va, aes(lon, soc, color = soc_measured)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x,        se = TRUE,
              aes(group = 1), color = "black",        linewidth = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE,
              aes(group = 1), color = "black", linetype = "dashed",
              linewidth = 0.8) +
  scale_color_manual(values = pal, name = "Type") +
  labs(title = "SOC vs. longitude", x = "Longitude", y = "SOC (%)") +
  theme(legend.position = "none")

p_lat <- ggplot(va, aes(lat, soc, color = soc_measured)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x,        se = TRUE,
              aes(group = 1), color = "black",        linewidth = 0.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE,
              aes(group = 1), color = "black", linetype = "dashed",
              linewidth = 0.8) +
  scale_color_manual(values = pal, name = "Type") +
  labs(title = "SOC vs. latitude", x = "Latitude", y = "SOC (%)") +
  theme(legend.position = "right")

ggsave("figures/02_trend_location.png",
       p_lon + p_lat + plot_layout(guides = "collect"),
       width = 10, height = 4, dpi = 150)

# ── 3. Normality: raw vs. log-transformed ────────────────────────────────────
p_hist_raw <- ggplot(va, aes(soc)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "SOC — raw", x = "SOC (%)", y = "Count")

p_hist_log <- ggplot(va, aes(log(soc))) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "SOC — log", x = "log(SOC)", y = "Count")

p_qq_raw <- ggplot(va, aes(sample = soc)) +
  stat_qq(alpha = 0.5) + stat_qq_line(color = "red") +
  labs(title = "Q-Q: raw SOC", x = "Theoretical", y = "Sample")

p_qq_log <- ggplot(va, aes(sample = log(soc))) +
  stat_qq(alpha = 0.5) + stat_qq_line(color = "red") +
  labs(title = "Q-Q: log(SOC)", x = "Theoretical", y = "Sample")

ggsave("figures/03_normality.png",
       (p_hist_raw | p_hist_log) / (p_qq_raw | p_qq_log),
       width = 9, height = 7, dpi = 150)

# ── 4. Measured vs. modeled comparison ───────────────────────────────────────
p_box <- ggplot(va, aes(soc_measured, soc, fill = soc_measured)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7, width = 0.5) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "SOC by type", x = NULL, y = "SOC (%)")

p_dens <- ggplot(va, aes(log(soc), fill = soc_measured, color = soc_measured)) +
  geom_density(alpha = 0.35) +
  scale_fill_manual(values = pal, name = "Type") +
  scale_color_manual(values = pal, name = "Type") +
  labs(title = "log(SOC) density by type", x = "log(SOC)", y = "Density")

# Map: spatial clustering by type?
p_map2 <- ggplot() +
  geom_polygon(data = va_state, aes(long, lat, group = group),
               fill = "grey92", color = "grey60", linewidth = 0.3) +
  geom_point(data = va, aes(lon, lat, color = soc_measured),
             size = 1.5, alpha = 0.7) +
  scale_color_manual(values = pal, name = "Type") +
  coord_fixed(1.3) +
  labs(title = "Sample locations by type", x = "Longitude", y = "Latitude")

ggsave("figures/04_measured_vs_modeled.png",
       (p_box | p_dens) / p_map2 + plot_layout(heights = c(1, 1.2)),
       width = 10, height = 8, dpi = 150)

# ── 5. Elevation as covariate ─────────────────────────────────────────────────

# 5a. Standalone elevation map
p_elev_map <- ggplot() +
  geom_polygon(data = va_state, aes(long, lat, group = group),
               fill = "grey92", color = "grey60", linewidth = 0.3) +
  geom_point(data = va, aes(lon, lat, color = elevation_m),
             size = 2, alpha = 0.85) +
  scale_color_viridis_c(name = "Elevation (m)", option = "cividis") +
  coord_fixed(1.3) +
  labs(title = "Elevation at Virginia RaCA sample sites",
       x = "Longitude", y = "Latitude")

ggsave("figures/05a_elevation_map.png", p_elev_map,
       width = 9, height = 4.5, dpi = 150)

# 5b. SOC vs elevation — 2-panel facet (raw | log)
va_long <- bind_rows(
  va |> mutate(response = soc,       scale = "SOC (%)"),
  va |> mutate(response = log(soc),  scale = "log(SOC)")
) |>
  mutate(scale = factor(scale, levels = c("SOC (%)", "log(SOC)")))

p_elev_facet <- ggplot(va_long, aes(elevation_m, response, color = soc_measured)) +
  geom_point(alpha = 0.45, size = 1.4) +
  geom_smooth(method = "loess", span = 0.6, se = TRUE,
              aes(group = 1), color = "black", linewidth = 0.9) +
  scale_color_manual(values = pal, name = "Type") +
  facet_wrap(~ scale, scales = "free_y", nrow = 1) +
  labs(title = "SOC vs. elevation",
       x = "Elevation (m)", y = NULL) +
  theme(strip.text = element_text(size = 11))

ggsave("figures/05b_elevation_soc.png", p_elev_facet,
       width = 10, height = 4.5, dpi = 150)

# ── 6. Proper 2D spatial trend surface assessment ────────────────────────────
# Analyze measured and modeled separately — they are different processes.
# For each group fit null / 1st order / 2nd order trend on log(soc) as a
# joint function of both coordinates, then compare via AIC and F-test.

va_m   <- filter(va, soc_measured == "measured") |> mutate(lsoc = log(soc))
va_mod <- filter(va, soc_measured == "modeled")  |> mutate(lsoc = log(soc))

fit_models <- function(df) {
  list(
    null = lm(lsoc ~ 1,                                                   data = df),
    ord1 = lm(lsoc ~ lon + lat,                                            data = df),
    ord2 = lm(lsoc ~ lon + lat + I(lon^2) + I(lat^2) + lon:lat,           data = df)
  )
}

fits_m   <- fit_models(va_m)
fits_mod <- fit_models(va_mod)

model_table <- function(fits, label) {
  nms <- c("null", "ord1", "ord2")
  data.frame(
    group  = label,
    model  = c("Null (constant)", "1st order", "2nd order"),
    df     = sapply(fits, \(f) df.residual(f)),
    AIC    = round(sapply(fits, AIC), 2),
    R2     = round(sapply(fits, \(f) summary(f)$r.squared), 3),
    pval_vs_prev = c(
      NA,
      round(anova(fits$null, fits$ord1)$`Pr(>F)`[2], 4),
      round(anova(fits$ord1, fits$ord2)$`Pr(>F)`[2], 4)
    )
  )
}

tbl <- rbind(model_table(fits_m, "measured"), model_table(fits_mod, "modeled"))
cat("\n── Trend surface model comparison (response: log SOC) ──\n")
print(tbl, row.names = FALSE)

# Prediction grid — clipped to Virginia state boundary via sf
library(sf)
sf_use_s2(FALSE)
va_poly <- st_as_sf(map("state", "virginia", plot = FALSE, fill = TRUE)) |>
  st_set_crs(4326)

lon_seq <- seq(min(va$lon) - 0.05, max(va$lon) + 0.05, length.out = 120)
lat_seq <- seq(min(va$lat) - 0.05, max(va$lat) + 0.05, length.out = 80)
grid    <- expand.grid(lon = lon_seq, lat = lat_seq)

# Keep only grid cells that fall inside Virginia
grid_sf  <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
in_va    <- st_within(grid_sf, va_poly, sparse = FALSE)[, 1]
grid     <- grid[in_va, ]

grid$trend1_m   <- predict(fits_m$ord1,   grid)
grid$trend2_m   <- predict(fits_m$ord2,   grid)
grid$trend1_mod <- predict(fits_mod$ord1, grid)
grid$trend2_mod <- predict(fits_mod$ord2, grid)

# Note: response is log(SOC) — transformation chosen in normality analysis
# (fig 03) because raw SOC is heavily right-skewed; all model-based
# analysis uses log scale for valid inference.
map_base <- list(
  geom_polygon(data = va_state, aes(long, lat, group = group),
               fill = NA, color = "grey50", linewidth = 0.4),
  coord_fixed(1.3),
  scale_fill_distiller(palette = "RdYlBu", name = "log(SOC)"),
  theme(axis.title = element_blank())
)

make_surface_plot <- function(grid_col, pts, title) {
  ggplot() +
    geom_tile(data = grid, aes(lon, lat, fill = .data[[grid_col]]),
              alpha = 0.9) +
    map_base +
    geom_point(data = pts, aes(lon, lat), size = 1.2,
               color = "black", alpha = 0.6) +
    labs(title = title)
}

p_s1m  <- make_surface_plot("trend1_m",   va_m,   "Measured — 1st order")
p_s2m  <- make_surface_plot("trend2_m",   va_m,   "Measured — 2nd order")
p_s1md <- make_surface_plot("trend1_mod", va_mod, "Modeled — 1st order")
p_s2md <- make_surface_plot("trend2_mod", va_mod, "Modeled — 2nd order")

ggsave("figures/06_trend_surfaces.png",
       (p_s1m | p_s2m) / (p_s1md | p_s2md) +
         plot_annotation(title = "Fitted spatial trend surfaces (log SOC)"),
       width = 12, height = 9, dpi = 150)

# ── 7. Residuals from best trend model mapped spatially ───────────────────────
# Spatial pattern in residuals → remaining structure for the variogram to capture

va_m$resid_null <- residuals(fits_m$null)
va_m$resid_ord1 <- residuals(fits_m$ord1)
va_m$resid_ord2 <- residuals(fits_m$ord2)

make_resid_map <- function(pts, resid_col, title) {
  ggplot() +
    geom_polygon(data = va_state, aes(long, lat, group = group),
                 fill = "grey92", color = "grey60", linewidth = 0.3) +
    geom_point(data = pts,
               aes(lon, lat, color = .data[[resid_col]], size = abs(.data[[resid_col]])),
               alpha = 0.8) +
    scale_color_distiller(palette = "RdBu", name = "Residual",
                          limits = \(x) c(-max(abs(x)), max(abs(x)))) +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    coord_fixed(1.3) +
    labs(title = title, x = "Longitude", y = "Latitude")
}

p_r0 <- make_resid_map(va_m, "resid_null", "Residuals: null (demeaned)")
p_r1 <- make_resid_map(va_m, "resid_ord1", "Residuals: 1st order trend")
p_r2 <- make_resid_map(va_m, "resid_ord2", "Residuals: 2nd order trend")

ggsave("figures/07_residual_maps.png",
       p_r0 / p_r1 / p_r2 +
         plot_annotation(title = "Measured SOC — spatial residuals by trend order"),
       width = 9, height = 13, dpi = 150)

cat("Figures saved to figures/\n")
cat("  01_map_soc.png              — spatial distribution\n")
cat("  02_trend_location.png       — marginal trend (EDA)\n")
cat("  03_normality.png            — raw vs log normality\n")
cat("  04_measured_vs_modeled.png  — observed vs simulated\n")
cat("  05_elevation_covariate.png  — elevation as covariate\n")
cat("  06_trend_surfaces.png       — fitted 1st/2nd order surfaces\n")
cat("  07_residual_maps.png        — spatial residuals by trend order\n")
