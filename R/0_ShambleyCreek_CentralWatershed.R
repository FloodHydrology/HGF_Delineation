# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Shambley Creek HGF Delineation
# Coder: Nate Jones (cnjones7@ua.edu)
# Date: 8/18/2026
# Purpose: Longitudinal profile of Shambley Creek main stem, compared against
#          Ashleigh's field-mapped hydrogeomorphic features (HGFs).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Section 1: Setup workspace
# Section 2: Delineate watershed
# Section 3: Extract main stem
# Section 4: Longitudinal profile vs. field HGFs

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
library(zoo)         # rolling smooth
library(patchwork)   # stack elevation / slope panels

# Initialize whitebox (only needs to run once per install)
# wbt_init()

# Working directory for whitebox I/O
wbt_wd <- tempdir()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.0 Delineate watershed ------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shambley Creek watershed outlet
outlet <- data.frame(lon = -88.013343, lat = 32.984109) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Download DEM (elevatr z = 14 ~ 10 m at this latitude), project to UTM 16N
dem <- get_elev_raster(st_buffer(outlet, dist = 5000), z = 14, clip = "bbox") %>%
  rast() %>%
  project("EPSG:26916")

dem_raw <- file.path(wbt_wd, "dem_raw.tif")
writeRaster(dem, dem_raw, overwrite = TRUE)

# Hydro-condition: feature-preserving smooth -> fill single-cell pits -> breach
dem_smooth <- file.path(wbt_wd, "dem_smooth.tif")
wbt_feature_preserving_smoothing(dem_raw, dem_smooth, filter = 11, wd = wbt_wd)

dem_fill <- file.path(wbt_wd, "dem_fill.tif")
wbt_fill_single_cell_pits(dem_smooth, dem_fill, wd = wbt_wd)

dem_breach <- file.path(wbt_wd, "dem_breach.tif")
wbt_breach_depressions_least_cost(dem_fill, dem_breach, dist = 100,
                                  fill = TRUE, wd = wbt_wd)

# Flow direction, accumulation, and stream raster (for pour-point snapping)
d8_pntr <- file.path(wbt_wd, "d8_pntr.tif")
wbt_d8_pointer(dem_breach, d8_pntr, wd = wbt_wd)

d8_accum <- file.path(wbt_wd, "d8_accum.tif")
wbt_d8_flow_accumulation(dem_breach, d8_accum, wd = wbt_wd)

streams <- file.path(wbt_wd, "streams.tif")
wbt_extract_streams(d8_accum, streams, threshold = 2000, wd = wbt_wd)

# Snap outlet to the stream network, then delineate everything upstream
outlet_utm <- st_transform(outlet, 26916)
outlet_shp <- file.path(wbt_wd, "outlet.shp")
st_write(outlet_utm, outlet_shp, delete_dsn = TRUE, quiet = TRUE)

outlet_snap <- file.path(wbt_wd, "outlet_snap.shp")
wbt_jenson_snap_pour_points(outlet_shp, streams, outlet_snap,
                            snap_dist = 100, wd = wbt_wd)

watershed <- file.path(wbt_wd, "watershed.tif")
wbt_watershed(d8_pntr, outlet_snap, watershed, wd = wbt_wd)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.0 Extract main stem --------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wbt_longest_flowpath returns the single longest channel in the basin -- no
# forking. This is the spine we build the longitudinal profile along.
longest <- file.path(wbt_wd, "longest.shp")
wbt_longest_flowpath(dem_breach, watershed, longest, wd = wbt_wd)

mainstem <- st_read(longest, quiet = TRUE) %>%
  st_set_crs(26916) %>%
  arrange(desc(LENGTH)) %>%   # LENGTH attribute written by the tool
  slice(1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.0 Longitudinal profile vs. field HGFs --------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HGFs are segments of the river corridor -- the whole valley cross-section
# (channel, riparian, floodplain, contributing hillslopes), not just the
# channel thread. Here we profile elevation and slope down the valley and
# compare against Ashleigh's field-mapped HGF reaches.

# 4.1 Sample elevation along the main stem -------------------------------------
step_m <- 10   # point spacing (m); ~ DEM resolution

mainstem_pts <- mainstem %>%
  st_line_sample(density = 1 / step_m) %>%
  st_cast("POINT") %>%
  st_as_sf() %>%
  mutate(elev = terra::extract(rast(dem_breach), vect(.))[, 2])

# 4.2 Build the distance-ordered profile ---------------------------------------
# Cumulative distance along the stem, oriented so distance increases from the
# outlet (lowest elevation) upstream.
coords <- st_coordinates(mainstem_pts)

profile <- mainstem_pts %>%
  st_drop_geometry() %>%
  mutate(x = coords[, 1], y = coords[, 2],
         dist = c(0, cumsum(sqrt(diff(x)^2 + diff(y)^2))))

if (profile$elev[1] > profile$elev[nrow(profile)]) {
  profile <- profile %>%
    arrange(desc(dist)) %>%
    mutate(dist = max(dist) - dist)
}

# 4.3 Compute slope ------------------------------------------------------------
# Point slope: raw two-point gradient between adjacent cells (the noisy signal,
# matching the point slopes shown in DMP WRR Fig 1C).
# Smoothed slope: local elevation~distance regression over a rolling window
# (the line we fit through that noise). smooth_win is the smoothing knob.
# Both taken as magnitude (|dz/dx|) so they share one sign convention -- the
# profile rises up-valley, so slope is reported positive.
smooth_win <- 10   # cells (~100 m at 10 m spacing); odd centers cleanly

profile <- profile %>%
  mutate(slope_pt = abs(c(NA, diff(elev) / diff(dist))))   # raw point-to-point

roll_slope <- zoo::rollapply(
  profile[, c("dist", "elev")], width = smooth_win,
  by.column = FALSE, fill = NA, align = "center",
  FUN = function(w) coef(lm(elev ~ dist, data = as.data.frame(w)))[["dist"]]
)

profile <- profile %>%
  mutate(slope_sm = abs(roll_slope))   # smoothed magnitude

# 4.4 Read + project Ashleigh's field HGFs onto the profile --------------------
# Central-site field polygons; filenames encode class (Inc/Wet/Typ).
ak_files <- list.files("data/StreamPoly", pattern = "^central.*\\.shp$",
                       full.names = TRUE)

ak_hgf <- ak_files %>%
  map(function(f) {
    nm <- tools::file_path_sans_ext(basename(f))
    st_read(f, quiet = TRUE) %>%
      mutate(
        source = nm,
        hgf = case_when(
          str_detect(nm, regex("Inc", ignore_case = TRUE)) ~ "Incised",
          str_detect(nm, regex("Wet", ignore_case = TRUE)) ~ "Wetland-stream",
          str_detect(nm, regex("Typ", ignore_case = TRUE)) ~ "Intact riparian",
          TRUE ~ "Unknown"
        )
      ) %>%
      select(source, hgf)
  }) %>%
  bind_rows() %>%
  st_transform(26916)

# Assign each profile point to the nearest field HGF (her polygons were drawn on
# a slightly offset DEM, so we snap to nearest rather than require containment).
# Points beyond max_snap are off her surveyed reaches -> Unclassified.
max_snap <- 50   # m; median snap ~14 m, cleanly separates on- vs off-reach

segment_pts <- profile %>% st_as_sf(coords = c("x", "y"), crs = 26916)
nn      <- st_nearest_feature(segment_pts, ak_hgf)
nn_dist <- as.numeric(st_distance(segment_pts, ak_hgf[nn, ], by_element = TRUE))

profile <- profile %>%
  mutate(
    hgf    = if_else(nn_dist > max_snap, "Unclassified", ak_hgf$hgf[nn]),
    source = if_else(nn_dist > max_snap, NA_character_, ak_hgf$source[nn])
  )

# Field HGF reach extents (bands) + the surveyed distance window
hgf_bands <- profile %>%
  filter(hgf != "Unclassified") %>%
  group_by(source, hgf) %>%
  summarise(d_start = min(dist), d_end = max(dist), .groups = "drop")

ak_window <- profile %>%
  filter(hgf != "Unclassified") %>%
  summarise(d_min = min(dist), d_max = max(dist))

# 4.5 Longitudinal profile figure ----------------------------------------------
hgf_cols <- c("Incised"         = "#d7191c",
              "Wetland-stream"  = "#2c7bb6",
              "Intact riparian" = "#fdae61")

# Elevation panel
p_elev <- profile %>%
  ggplot(aes(x = dist, y = elev)) +
  geom_rect(data = hgf_bands, inherit.aes = FALSE,
            aes(xmin = d_start, xmax = d_end,
                ymin = -Inf, ymax = Inf, fill = hgf),
            alpha = 0.35, color = NA) +
  scale_fill_manual(values = hgf_cols, name = "Field HGF") +
  geom_line(lwd = 0.75, col = "black") +
  coord_cartesian(xlim = c(ak_window$d_min, ak_window$d_max)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10)) +
  xlab(NULL) +
  ylab("Elevation [masl]")

# Slope panel -- raw point slopes (grey scatter) under the smoothed slope
# (black line), matching DMP WRR Fig 1C. Linear y-axis so the near-zero and
# occasionally reversed point slopes at the flats stay visible.
p_slope <- profile %>%
  ggplot(aes(x = dist, y = slope_sm)) +
  geom_rect(data = hgf_bands, inherit.aes = FALSE,
            aes(xmin = d_start, xmax = d_end,
                ymin = -Inf, ymax = Inf, fill = hgf),
            alpha = 0.35, color = NA) +
  scale_fill_manual(values = hgf_cols, name = "Field HGF") +
  geom_point(aes(y = slope_pt), size = 1.2, alpha = 0.40, col = "grey30") +
  geom_line(lwd = 0.75, col = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10)) +
  coord_cartesian(xlim = c(ak_window$d_min, ak_window$d_max),
                  ylim = c(0, 0.045)) +
  xlab("Distance from outlet [m]") +
  ylab("Slope [m/m]")

# Stack, shared legend
fig <- (p_elev / p_slope) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Central Watershed",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )
fig

# 4.6 Export -------------------------------------------------------------------
ggsave("output/shambley_longitudinal_profile.png", fig,
       width = 8, height = 5, dpi = 300)
