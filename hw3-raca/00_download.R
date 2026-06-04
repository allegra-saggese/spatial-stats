# 00_download.R
# Downloads USGS NED elevation for Virginia and merges with VA soil carbon data.
# Elevation is extracted at each sample point (both measured and modeled SOC).
# API key stored in ~/.Renviron as OPENTOPO_API_KEY

library(readr)
library(httr2)
library(terra)
library(dplyr)

# ── Load VA soil carbon data ──────────────────────────────────────────────────
va <- read_csv("VA_soil_carbon.csv", show_col_types = FALSE)

cat(sprintf("Loaded %d rows (%d measured, %d modeled); %d missing elevation\n",
            nrow(va),
            sum(va$soc_measured == "measured"),
            sum(va$soc_measured == "modeled"),
            sum(is.na(va$elevation))))

# ── Download USGSNED raster for Virginia ──────────────────────────────────────
# Bounding box derived from data extent + 0.1° buffer
api_key <- Sys.getenv("OPENTOPO_API_KEY")
if (nchar(api_key) == 0) stop("OPENTOPO_API_KEY not found in ~/.Renviron")

west  <- floor(min(va$long) * 10) / 10 - 0.1
east  <- ceiling(max(va$long) * 10) / 10 + 0.1
south <- floor(min(va$lat)  * 10) / 10 - 0.1
north <- ceiling(max(va$lat)  * 10) / 10 + 0.1

cat(sprintf("Downloading USGSNED for bbox [%.2f, %.2f, %.2f, %.2f]...\n",
            west, south, east, north))

dem_file <- tempfile(fileext = ".tif")

resp <- request("https://portal.opentopography.org/API/globaldem") |>
  req_url_query(
    demtype      = "SRTMGL1",   # NASA/USGS SRTM ~30m, same resolution as USGSNED
    west         = west,
    south        = south,
    east         = east,
    north        = north,
    outputFormat = "GTiff",
    API_Key      = api_key
  ) |>
  req_timeout(600) |>
  req_perform()

writeBin(resp_body_raw(resp), dem_file)
cat(sprintf("Downloaded: %.1f MB\n", file.size(dem_file) / 1e6))

# ── Extract elevation at sample points ────────────────────────────────────────
dem  <- rast(dem_file)
pts  <- vect(va, geom = c("long", "lat"), crs = "EPSG:4326")
vals <- terra::extract(dem, pts)[, 2]

va <- va |>
  mutate(elevation = as.numeric(vals))

cat(sprintf("Elevation filled: %d / %d rows (%.1f%%)\n",
            sum(!is.na(va$elevation)), nrow(va),
            100 * mean(!is.na(va$elevation))))

# ── Save ──────────────────────────────────────────────────────────────────────
write_csv(va, "VA_soil_carbon.csv")
cat("Saved: VA_soil_carbon.csv\n")

unlink(dem_file)
