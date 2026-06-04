# 01_merge.R
# Loads VA_soil_carbon.csv (with elevation filled by 00_download.R)
# and prepares a clean analysis-ready dataset.

library(readr)
library(dplyr)

va <- read_csv("VA_soil_carbon.csv", show_col_types = FALSE)

# Split into observed (lab-measured SOC) and simulated (modeled SOC)
va_measured <- va |> filter(soc_measured == "measured")
va_modeled  <- va |> filter(soc_measured == "modeled")

cat(sprintf("Measured: %d rows | Modeled: %d rows\n",
            nrow(va_measured), nrow(va_modeled)))
cat(sprintf("Elevation complete: %d / %d\n",
            sum(!is.na(va$elevation)), nrow(va)))

# Analysis-ready dataset: all rows with key variables named clearly
va_clean <- va |>
  select(
    sample_id,
    rcapid,
    soc,
    soc_sd,
    soc_measured,
    depth_top    = sample_top,
    depth_bot    = sample_bottom,
    texture,
    elevation_m  = elevation,
    lon          = long,
    lat,
    landuse
  )

write_csv(va_clean, "VA_soil_carbon_clean.csv")
cat("Saved: VA_soil_carbon_clean.csv\n")
