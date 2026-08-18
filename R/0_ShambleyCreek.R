# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Shambley Creek HGF Delineation
# Coder: Nate Jones (cnjones7@ua.edu)
# Date: 8/18/2026
# Purpose: Expore longitudinal profile of Shambley creek
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Section 1: Setup worksapce
# Section 2: Delineate watershed
# Section 3: Delineate stream network

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.0 Setup workspace ----------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clear workspace
remove(list = ls())

# Load packages
library(tidyverse)   # join the cult (pipe, dplyr, etc.)
library(sf)
library(terra)
library(elevatr)
library(whitebox)
library(mapview)

# Initialize whitebox (only needs to run once per install)
# wbt_init()

# Define working directory for whitebox I/O
wbt_wd <- tempdir()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.0 Delineate watershed ------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1 Create spatial point from outlet coordinates -----------------------------
# Shambley Creek watershed outlet
outlet <- data.frame(
  lon = -88.013343,
  lat =  32.984109
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# 2.2 Download DEM with elevatr ------------------------------------------------
# z = 12-14 for watershed-scale work; z = 14 ~ 1 arc-sec / ~10 m at this latitude.
# Pull a buffered area so the full contributing watershed is captured.
dem <- get_elev_raster(
  locations = st_buffer(outlet, dist = 5000),  # 5 km buffer around outlet
  z = 14,
  clip = "bbox"
)

# Coerce to terra SpatRaster and project to a metric CRS (UTM 16N for W AL/E MS)
dem <- rast(dem) %>%
  project("EPSG:26916")

# Write raw DEM to disk for whitebox
dem_raw <- file.path(wbt_wd, "dem_raw.tif")
writeRaster(dem, dem_raw, overwrite = TRUE)

# 2.3 Smooth the DEM -----------------------------------------------------------
# Feature-preserving smoothing to reduce noise before hydro-conditioning
dem_smooth <- file.path(wbt_wd, "dem_smooth.tif")
wbt_feature_preserving_smoothing(
  dem    = dem_raw,
  output = dem_smooth,
  filter = 11,
  wd     = wbt_wd
)

# 2.4 Fill single-cell pits ----------------------------------------------------
dem_fill <- file.path(wbt_wd, "dem_fill.tif")
wbt_fill_single_cell_pits(
  dem    = dem_smooth,
  output = dem_fill,
  wd     = wbt_wd
)

# 2.5 Breach depressions -------------------------------------------------------
# Breaching preserves flow paths better than filling for larger depressions
dem_breach <- file.path(wbt_wd, "dem_breach.tif")
wbt_breach_depressions_least_cost(
  dem    = dem_fill,
  output = dem_breach,
  dist   = 100,
  fill   = TRUE,
  wd     = wbt_wd
)

# 2.6 Delineate watershed upstream of outlet -----------------------------------
# Flow direction (D8)
d8_pntr <- file.path(wbt_wd, "d8_pntr.tif")
wbt_d8_pointer(
  dem    = dem_breach,
  output = d8_pntr,
  wd     = wbt_wd
)

# Flow accumulation (to build a stream raster for snapping the pour point)
d8_accum <- file.path(wbt_wd, "d8_accum.tif")
wbt_d8_flow_accumulation(
  input  = dem_breach,
  output = d8_accum,
  wd     = wbt_wd
)

# Extract streams (threshold in n cells — tune to Shambley; placeholder 5000)
streams <- file.path(wbt_wd, "streams.tif")
wbt_extract_streams(
  flow_accum = d8_accum,
  output     = streams,
  threshold  = 2000,
  wd         = wbt_wd
)

# Write outlet point to disk (project to match DEM first)
outlet_utm <- st_transform(outlet, 26916)
outlet_shp <- file.path(wbt_wd, "outlet.shp")
st_write(outlet_utm, outlet_shp, delete_dsn = TRUE, quiet = TRUE)

# Snap pour point to the stream network so it sits on a high-accumulation cell
outlet_snap <- file.path(wbt_wd, "outlet_snap.shp")
wbt_jenson_snap_pour_points(
  pour_pts = outlet_shp,
  streams  = streams,
  output   = outlet_snap,
  snap_dist = 100,
  wd       = wbt_wd
)

# Delineate everything upstream of the (snapped) outlet
watershed <- file.path(wbt_wd, "watershed.tif")
wbt_watershed(
  d8_pntr = d8_pntr,
  pour_pts = outlet_snap,
  output   = watershed,
  wd       = wbt_wd
)

# Read back in and vectorize for inspection
ws_rast <- rast(watershed)
ws_poly <- as.polygons(ws_rast) %>% st_as_sf()

# Quick check
plot(dem)
plot(st_geometry(ws_poly), add = TRUE, border = "red", lwd = 2)
plot(st_geometry(outlet_utm), add = TRUE, pch = 16, col = "blue")

# 2.7 Interactive visualization ------------------------------------------------
# Read snapped outlet back in for comparison against the original
outlet_snap_sf <- st_read(outlet_snap, quiet = TRUE)

mapview(ws_poly,        col.regions = "steelblue", alpha.regions = 0.3,
        layer.name = "Watershed") +
  mapview(outlet_utm,     col.regions = "red",   cex = 6, layer.name = "Outlet (raw)") +
  mapview(outlet_snap_sf, col.regions = "yellow", cex = 6, layer.name = "Outlet (snapped)")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.0 Delineate stream network -------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The D8 pointer, flow accumulation, and stream raster already exist from
# Section 1 (built for pour-point snapping). Here we clip the stream network to
# the delineated watershed, vectorize it, and tag Strahler stream order.

# 2.1 Clip stream raster to watershed ------------------------------------------
# Mask the watershed-wide stream raster to just our basin so we don't carry
# neighboring drainages into the network.
streams_clip <- file.path(wbt_wd, "streams_clip.tif")
wbt_multiply(
  input1 = streams,
  input2 = watershed,   # watershed raster = 1 inside basin, NA outside
  output = streams_clip,
  wd     = wbt_wd
)

# 2.2 Strahler stream order ----------------------------------------------------
strahler <- file.path(wbt_wd, "strahler.tif")
wbt_strahler_stream_order(
  d8_pntr = d8_pntr,
  streams = streams_clip,
  output  = strahler,
  wd      = wbt_wd
)

# 2.3 Convert stream raster to vector ------------------------------------------
# Vectorize the ordered network. wbt writes the Strahler value into the
# STRM_VAL / FID attributes of the resulting lines.
streams_vec <- file.path(wbt_wd, "streams_vec.shp")
wbt_raster_streams_to_vector(
  streams = strahler,
  d8_pntr = d8_pntr,
  output  = streams_vec,
  wd      = wbt_wd
)

# 2.4 Read network back in -----------------------------------------------------
# wbt_raster_streams_to_vector doesn't write a .prj, so set the CRS to match
# the DEM (UTM 16N) explicitly.
streams_sf <- st_read(streams_vec, quiet = TRUE) %>%
  st_set_crs(26916)

# 2.5 Interactive visualization ------------------------------------------------
mapview(ws_poly,     col.regions = "steelblue", alpha.regions = 0.3,
        layer.name = "Watershed") +
  mapview(streams_sf, zcol = "STRM_VAL", legend = TRUE,
          layer.name = "Stream order") +
  mapview(outlet_snap_sf, col.regions = "yellow", cex = 6,
          layer.name = "Outlet (snapped)")









