.libPaths(c(.libPaths(), "/usr/lib/R/site-library"))
library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(png)
library(grid)

tokens <- c(Sys.getenv("BIRDWEATHER_TOKEN"), Sys.getenv("BIRDWEATHER_TOKEN_2"))
tokens <- tokens[nchar(tokens) > 0]
stopifnot("Missing BIRDWEATHER_TOKEN secret" = length(tokens) > 0)

# 1. Fetch 48-hour species summary for a single station token
fetch_station_species <- function(tok) {
  species_url <- paste0(
    "https://", paste(c("app", "birdweather", "com"), collapse = "."), "/",
    paste(c("api", "v1", "stations", tok, "species"), collapse = "/")
  )
  
  cutoff_48h <- format(as.POSIXct(Sys.time() - as.difftime(48, units = "hours"), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  
  res <- GET(
    url = species_url,
    query = list(from = cutoff_48h, limit = 200),
    add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0")
  )
  stop_for_status(res)
  parsed <- fromJSON(content(res, as = "text", encoding = "UTF-8"), flatten = TRUE)
  sp_df <- if (!is.null(parsed$species)) parsed$species else as.data.frame(parsed)
  
  if (is.null(sp_df) || nrow(sp_df) == 0) return(tibble())
  
  common_col <- grep("commonName|common_name", names(sp_df), value = TRUE)
  sci_col    <- grep("scientificName|scientific_name", names(sp_df), value = TRUE)
  count_col  <- grep("detections_count|detections.total|count|total", names(sp_df), value = TRUE)
  
  sp_df %>%
    transmute(
      common_name = as.character(.data[[common_col]]),
      sci_name    = as.character(.data[[sci_col]]),
      detections  = as.numeric(.data[[count_col]])
    ) %>%
    filter(!is.na(common_name), !is.na(sci_name), detections > 0)
}

# Fetch and combine detections across both devices (keeping BOTH birds and bats)
all_stations_df <- bind_rows(lapply(tokens, fetch_station_species))
stopifnot("No detections found across stations" = nrow(all_stations_df) > 0)

flock_df <- all_stations_df %>%
  mutate(clean_sci = sub("^(\\S+\\s+\\S+).*", "\\1", sci_name)) %>%
  group_by(common_name, clean_sci) %>%
  summarise(detections = sum(detections), .groups = "drop") %>%
  mutate(
    is_bat = grepl("\\bbat\\b|myotis|flying fox", common_name, ignore.case = TRUE) |
             grepl("^(Myotis|Lasionycteris|Eptesicus|Lasiurus|Corynorhinus|Antrozous|Tadarida)\\b", clean_sci)
  ) %>%
  arrange(desc(detections)) %>%
  mutate(
    pose = ifelse(row_number() %% 2 == 1, 1, 2),
    scale_r = (detections^0.32)
  )

flock_df <- flock_df %>%
  mutate(
    radius = 0.90 + 1.40 * ((scale_r - min(scale_r)) / max(1e-6, max(scale_r) - min(scale_r))),
    label_txt = paste0(common_name, " (", detections, ")")
  )

# 2. Fetch official 6-letter eBird species codes for birds
ebird_tax_url <- paste0("https://", paste(c("api", "ebird", "org"), collapse = "."), "/v2/ref/taxonomy/ebird")
tax_res <- try(GET(
  url = ebird_tax_url,
  query = list(fmt = "json", cat = "species"),
  add_headers("User-Agent" = "Mozilla/5.0"),
  timeout(15)
), silent = TRUE)

ebird_lookup <- tibble(sciName = character(), comName = character(), speciesCode = character())
if (!inherits(tax_res, "try-error") && status_code(tax_res) == 200) {
  ebird_lookup <- fromJSON(content(tax_res, as = "text", encoding = "UTF-8")) %>%
    select(sciName, comName, speciesCode)
}

flock_df$ebird_code <- vapply(seq_len(nrow(flock_df)), function(idx) {
  if (flock_df$is_bat[idx]) return("")
  sci <- flock_df$clean_sci[idx]
  com <- flock_df$common_name[idx]
  m1 <- ebird_lookup$speciesCode[ebird_lookup$sciName == sci]
  if (length(m1) > 0 && !is.na(m1[1])) return(m1[1])
  m2 <- ebird_lookup$speciesCode[tolower(ebird_lookup$comName) == tolower(com)]
  if (length(m2) > 0 && !is.na(m2[1])) return(m2[1])
  return("")
}, character(1))

# 3. Fetch full-color woodblock PNGs for birds with gentle pacing & retries
get_bird_art <- function(sci_name, common_name, pose = 1, is_bat = FALSE) {
  if (is_bat) return(NULL)
  cutout_url <- paste0(
    "https://", paste(c("bird", "onethreenine", "net"), collapse = "."),
    "/avian/api/cutout.php"
  )
  
  poses_to_try <- unique(c(pose, 1, 2))
  for (p_try in poses_to_try) {
    for (attempt in 1:2) {
      Sys.sleep(0.15)
      res <- try(GET(
        url = cutout_url,
        query = list(sci = sci_name, com = common_name, pose = p_try),
        add_headers("User-Agent" = "Mozilla/5.0"),
        timeout(15)
      ), silent = TRUE)
      
      if (!inherits(res, "try-error") && status_code(res) == 200) {
        img <- try(readPNG(content(res, as = "raw")), silent = TRUE)
        if (!inherits(img, "try-error")) return(img)
      }
      Sys.sleep(0.35)
    }
  }
  return(NULL)
}

rasters <- mapply(get_bird_art, flock_df$clean_sci, flock_df$common_name, flock_df$pose, flock_df$is_bat, SIMPLIFY = FALSE)
has_art <- !sapply(rasters, is.null)

# 4. Compute exact rectangular bounding box
n <- nrow(flock_df)
bird_w <- numeric(n)
bird_h <- numeric(n)
box_hw <- numeric(n)
box_hh <- numeric(n)
label_band <- 0.45

for (i in seq_len(n)) {
  if (has_art[i]) {
    h_px <- nrow(rasters[[i]])
    w_px <- ncol(rasters[[i]])
    sf <- (flock_df$radius[i] * 1.95) / sqrt(h_px^2 + w_px^2)
    bird_w[i] <- w_px * sf
    bird_h[i] <- h_px * sf
  } else if (flock_df$is_bat[i]) {
    bird_w[i] <- flock_df$radius[i] * 1.55
    bird_h[i] <- flock_df$radius[i] * 0.85
  } else {
    bird_w[i] <- flock_df$radius[i] * 1.25
    bird_h[i] <- flock_df$radius[i] * 1.05
  }
  est_text_w <- nchar(flock_df$label_txt[i]) * 0.165
  box_hw[i] <- max(bird_w[i], est_text_w) / 2 + 0.12
  box_hh[i] <- (bird_h[i] + label_band) / 2 + 0.10
}

# 5. Concentric Center-First Rectangle Packing
cx <- numeric(n)
cy <- numeric(n)

if (n > 1) {
  radii_steps <- seq(0.1, 30, by = 0.08)
  angles <- seq(0, 2 * pi, length.out = 180)[-180]
  
  for (i in 2:n) {
    placed <- FALSE
    rot_angles <- angles[((seq_along(angles) + (i * 37)) %% length(angles)) + 1]
    
    for (rad in radii_steps) {
      for (ang in rot_angles) {
        tx <- rad * cos(ang) * 1.55
        ty <- rad * sin(ang)
        overlap_x <- abs(tx - cx[1:(i - 1)]) < (box_hw[i] + box_hw[1:(i - 1)])
        overlap_y <- abs(ty - cy[1:(i - 1)]) < (box_hh[i] + box_hh[1:(i - 1)])
        if (!any(overlap_x & overlap_y)) {
          cx[i] <- tx
          cy[i] <- ty
          placed <- TRUE
          break
        }
      }
      if (placed) break
    }
  }
}

flock_df <- flock_df %>%
  mutate(
    x_wide = cx,
    y = cy,
    bird_w = bird_w,
    bird_h = bird_h,
    box_hw = box_hw,
    box_hh = box_hh,
    has_art = has_art
  )

x_limits <- c(min(flock_df$x_wide - flock_df$box_hw) - 0.4, max(flock_df$x_wide + flock_df$box_hw) + 0.4)
y_limits <- c(min(flock_df$y - flock_df$box_hh) - 0.3, max(flock_df$y + flock_df$box_hh) + 0.3)

# Generic perched bird silhouette coordinates
bird_poly_template <- tibble(
  px = c(
    -0.48, -0.25, -0.08,  0.08,  0.22,  0.34,  0.48,  0.35,  0.32,
     0.24,  0.10,  0.06,  0.12,  0.02, -0.04, -0.02, -0.16, -0.45, -0.36
  ),
  py = c(
    -0.28, -0.08,  0.14,  0.22,  0.38,  0.34,  0.28,  0.22,  0.08,
    -0.10, -0.22, -0.44, -0.46, -0.46, -0.24, -0.45, -0.22, -0.35, -0.26
  )
)

# Flying bat silhouette coordinates
bat_poly_template <- tibble(
  px = c(
     0.00, -0.04, -0.07, -0.14, -0.32, -0.50, -0.38, -0.26, -0.14,
    -0.06,  0.00,  0.06,  0.14,  0.26,  0.38,  0.50,  0.32,  0.14,
     0.07,  0.04
  ),
  py = c(
     0.18,  0.34,  0.22,  0.32,  0.42,  0.12, -0.12, -0.05, -0.22,
    -0.15, -0.35, -0.15, -0.22, -0.05, -0.12,  0.12,  0.42,  0.32,
     0.22,  0.34
  )
)

# 6. Render warm cream woodblock collage
bc_time <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%b %d, %I:%M %p %Z")
total_calls <- sum(flock_df$detections)

p <- ggplot(flock_df, aes(x = x_wide, y = y)) +
  geom_blank() +
  coord_equal(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "off") +
  labs(
    title = "Myrtle Rocks Bird & Bat Detections",
    subtitle = paste0("Last 48 Hours (", length(tokens), " Stations): ", nrow(flock_df), " Species (", total_calls, " Detections) • Updated ", bc_time)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#f6f0e2", color = NA),
    panel.background = element_rect(fill = "#f6f0e2", color = NA),
    plot.title = element_text(
      hjust = 0.5, family = "serif", size = 20, face = "bold",
      color = "#2b241d", margin = margin(t = 12, b = 2)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, family = "serif", size = 11.5, color = "#6e6051", margin = margin(b = 8)
    ),
    plot.margin = margin(12, 16, 14, 16)
  )

for (i in seq_len(n)) {
  y_bird_center <- flock_df$y[i] + (label_band / 2)
  if (flock_df$has_art[i]) {
    img <- rasters[[i]]
    half_w <- flock_df$bird_w[i] / 2
    half_h <- flock_df$bird_h[i] / 2
    
    p <- p + annotation_raster(
      img,
      xmin = flock_df$x_wide[i] - half_w,
      xmax = flock_df$x_wide[i] + half_w,
      ymin = y_bird_center - half_h,
      ymax = y_bird_center + half_h
    )
  } else if (flock_df$is_bat[i]) {
    poly_i <- bat_poly_template %>%
      mutate(
        x = flock_df$x_wide[i] + px * flock_df$bird_w[i],
        y = y_bird_center + py * flock_df$bird_h[i]
      )
    p <- p + geom_polygon(
      data = poly_i,
      aes(x = x, y = y),
      fill = "#4a4e5c",
      color = "#262933",
      linewidth = 0.45
    )
  } else {
    poly_i <- bird_poly_template %>%
      mutate(
        x = flock_df$x_wide[i] + px * flock_df$bird_w[i],
        y = y_bird_center + py * flock_df$bird_h[i]
      )
    p <- p + geom_polygon(
      data = poly_i,
      aes(x = x, y = y),
      fill = "#7d6b57",
      color = "#3b3228",
      linewidth = 0.45
    )
  }
}

p <- p + geom_text(
  data = flock_df,
  aes(x = x_wide, y = y - (bird_h / 2) - 0.06, label = label_txt),
  size = 2.9, family = "serif", color = "#3b3228", vjust = 1
)

ggsave("HeckfordBirds.png", plot = p, width = 13, height = 8, dpi = 220, bg = "#f6f0e2")

