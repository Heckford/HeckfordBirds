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

fetch_station_species <- function(tok, from_iso) {
  species_url <- paste0(
    "https", "://", paste(c("app", "birdweather", "com"), collapse = "."), "/",
    paste(c("api", "v1", "stations", tok, "species"), collapse = "/")
  )
  res <- GET(url = species_url, query = list(from = from_iso, limit = 500), add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0"))
  stop_for_status(res)
  parsed <- fromJSON(content(res, as = "text", encoding = "UTF-8"), flatten = TRUE)
  sp_df <- if (!is.null(parsed$species)) parsed$species else as.data.frame(parsed)
  if (is.null(sp_df) || nrow(sp_df) == 0) return(tibble())
  
  common_col <- grep("commonName|common_name", names(sp_df), value = TRUE)
  sci_col    <- grep("scientificName|scientific_name", names(sp_df), value = TRUE)
  count_col  <- grep("detections_count|detections.total|count|total", names(sp_df), value = TRUE)
  
  df <- tibble(common_name = as.character(sp_df[[common_col]]), sci_name = as.character(sp_df[[sci_col]]), detections = as.numeric(sp_df[[count_col]]))
  df[!is.na(df$common_name) & !is.na(df$sci_name) & df$detections > 0, ]
}

ebird_tax_url <- paste0("https", "://", paste(c("api", "ebird", "org"), collapse = "."), "/v2/ref/taxonomy/ebird")
tax_res <- try(GET(url = ebird_tax_url, query = list(fmt = "json", cat = "species"), add_headers("User-Agent" = "Mozilla/5.0"), timeout(15)), silent = TRUE)
ebird_lookup <- data.frame(sciName = character(), comName = character(), speciesCode = character(), stringsAsFactors = FALSE)
if (!inherits(tax_res, "try-error") && status_code(tax_res) == 200) ebird_lookup <- fromJSON(content(tax_res, as = "text", encoding = "UTF-8"))

art_cache <- new.env(parent = emptyenv())
get_bird_art <- function(sci_name, common_name, pose = 1, is_bat = FALSE) {
  if (is_bat) return(NULL)
  cache_key <- paste0(sci_name, "_", pose)
  if (exists(cache_key, envir = art_cache)) return(get(cache_key, envir = art_cache))
  
  cutout_url <- paste0("https", "://", paste(c("bird", "onethreenine", "net"), collapse = "."), "/avian/api/cutout.php")
  for (p_try in unique(c(pose, 1, 2))) {
    for (attempt in 1:2) {
      Sys.sleep(0.12)
      res <- try(GET(url = cutout_url, query = list(sci = sci_name, com = common_name, pose = p_try), add_headers("User-Agent" = "Mozilla/5.0"), timeout(15)), silent = TRUE)
      if (!inherits(res, "try-error") && status_code(res) == 200) {
        img <- try(readPNG(content(res, as = "raw")), silent = TRUE)
        if (!inherits(img, "try-error")) {
          assign(cache_key, img, envir = art_cache)
          return(img)
        }
      }
      Sys.sleep(0.3)
    }
  }
  assign(cache_key, NULL, envir = art_cache)
  return(NULL)
}

bird_poly <- data.frame(px = c(-0.48,-0.25,-0.08,0.08,0.22,0.34,0.48,0.35,0.32,0.24,0.10,0.06,0.12,0.02,-0.04,-0.02,-0.16,-0.45,-0.36), py = c(-0.28,-0.08,0.14,0.22,0.38,0.34,0.28,0.22,0.08,-0.10,-0.22,-0.44,-0.46,-0.46,-0.24,-0.45,-0.22,-0.35,-0.26))
bat_poly  <- data.frame(px = c(0.00,-0.04,-0.07,-0.14,-0.32,-0.50,-0.38,-0.26,-0.14,-0.06,0.00,0.06,0.14,0.26,0.38,0.50,0.32,0.14,0.07,0.04), py = c(0.18,0.34,0.22,0.32,0.42,0.12,-0.12,-0.05,-0.22,-0.15,-0.35,-0.15,-0.22,-0.05,-0.12,0.12,0.42,0.32,0.22,0.34))
build_collage <- function(from_iso, period_label, png_file, json_file, compact_mode = FALSE) {
  raw_df <- bind_rows(lapply(tokens, fetch_station_species, from_iso = from_iso))
  stopifnot("No detections found" = nrow(raw_df) > 0)
  
  raw_df$clean_sci <- sub("^(\\S+\\s+\\S+).*", "\\1", raw_df$sci_name)
  flock_df <- aggregate(detections ~ common_name + clean_sci, data = raw_df, FUN = sum)
  flock_df <- flock_df[order(-flock_df$detections), ]
  rownames(flock_df) <- NULL
  
  flock_df$is_bat <- grepl("\\bbat\\b|myotis|flying fox", flock_df$common_name, ignore.case = TRUE) |
                     grepl("^(Myotis|Lasionycteris|Eptesicus|Lasiurus|Corynorhinus|Antrozous|Tadarida)\\b", flock_df$clean_sci)
  flock_df$pose <- ifelse(seq_len(nrow(flock_df)) %% 2 == 1, 1, 2)
  flock_df$scale_r <- flock_df$detections^ifelse(compact_mode, 0.26, 0.32)
  
  r_min <- ifelse(compact_mode, 0.68, 0.90)
  r_span <- ifelse(compact_mode, 0.95, 1.40)
  flock_df$radius <- r_min + r_span * ((flock_df$scale_r - min(flock_df$scale_r)) / max(1e-6, max(flock_df$scale_r) - min(flock_df$scale_r)))
  flock_df$label_txt <- paste0(flock_df$common_name, " (", flock_df$detections, ")")
  
  flock_df$ebird_code <- vapply(seq_len(nrow(flock_df)), function(idx) {
    if (flock_df$is_bat[idx]) return("")
    m1 <- ebird_lookup$speciesCode[ebird_lookup$sciName == flock_df$clean_sci[idx]]
    if (length(m1) > 0 && !is.na(m1)) return(m1)
    m2 <- ebird_lookup$speciesCode[tolower(ebird_lookup$comName) == tolower(flock_df$common_name[idx])]
    if (length(m2) > 0 && !is.na(m2)) return(m2)
    return("")
  }, character(1))
  
  rasters <- mapply(get_bird_art, flock_df$clean_sci, flock_df$common_name, flock_df$pose, flock_df$is_bat, SIMPLIFY = FALSE)
  has_art <- !sapply(rasters, is.null)
  
  n <- nrow(flock_df)
  bird_w <- numeric(n); bird_h <- numeric(n); box_hw <- numeric(n); box_hh <- numeric(n)
  label_band <- ifelse(compact_mode, 0.36, 0.45)
  char_w <- ifelse(compact_mode, 0.13, 0.165)
  
  for (i in seq_len(n)) {
    if (has_art[i]) {
      sf <- (flock_df$radius[i] * 1.95) / sqrt(nrow(rasters[[i]])^2 + ncol(rasters[[i]])^2)
      bird_w[i] <- ncol(rasters[[i]]) * sf; bird_h[i] <- nrow(rasters[[i]]) * sf
    } else {
      bird_w[i] <- flock_df$radius[i] * ifelse(flock_df$is_bat[i], 1.55, 1.25)
      bird_h[i] <- flock_df$radius[i] * ifelse(flock_df$is_bat[i], 0.85, 1.05)
    }
    box_hw[i] <- max(bird_w[i], nchar(flock_df$label_txt[i]) * char_w) / 2 + 0.10
    box_hh[i] <- (bird_h[i] + label_band) / 2 + 0.08
  }
  
  cx <- numeric(n); cy <- numeric(n)
  if (n > 1) {
    radii_steps <- seq(0.1, 45, by = 0.08)
    angles <- seq(0, 2 * pi, length.out = 180)[-180]
    for (i in 2:n) {
      placed <- FALSE
      rot_angles <- angles[((seq_along(angles) + (i * 37)) %% length(angles)) + 1]
      for (rad in radii_steps) {
        for (ang in rot_angles) {
          tx <- rad * cos(ang) * 1.62; ty <- rad * sin(ang)
          if (!any((abs(tx - cx[1:(i - 1)]) < (box_hw[i] + box_hw[1:(i - 1)])) & (abs(ty - cy[1:(i - 1)]) < (box_hh[i] + box_hh[1:(i - 1)])))) {
            cx[i] <- tx; cy[i] <- ty; placed <- TRUE; break
          }
        }
        if (placed) break
      }
    }
  }
  
  flock_df$x_wide <- cx; flock_df$y <- cy
  x_min <- min(cx - box_hw) - 0.4; x_max <- max(cx + box_hw) + 0.4
  y_min <- min(cy - box_hh) - 0.3; y_max <- max(cy + box_hh) + 0.3
  
  bc_time <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%b %d, %I:%M %p %Z")
  txt_size <- ifelse(compact_mode, 2.35, 2.9)
  
  p <- ggplot(flock_df, aes(x = x_wide, y = y)) + geom_blank() +
    coord_equal(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE, clip = "off") +
    labs(title = "Myrtle Rocks Bird & Bat Detections", subtitle = paste0(period_label, " (", length(tokens), " Stations): ", n, " Species (", sum(flock_df$detections), " Detections) • Updated ", bc_time)) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "#f6f0e2", color = NA), panel.background = element_rect(fill = "#f6f0e2", color = NA),
      plot.title = element_text(hjust = 0.5, family = "serif", size = 20, face = "bold", color = "#2b241d", margin = margin(t = 12, b = 2)),
      plot.subtitle = element_text(hjust = 0.5, family = "serif", size = 11.5, color = "#6e6051", margin = margin(b = 8)),
      plot.margin = margin(12, 16, 14, 16)
    )
  
  for (i in seq_len(n)) {
    yc <- cy[i] + (label_band / 2)
    if (has_art[i]) {
      p <- p + annotation_raster(rasters[[i]], xmin = cx[i] - bird_w[i]/2, xmax = cx[i] + bird_w[i]/2, ymin = yc - bird_h[i]/2, ymax = yc + bird_h[i]/2)
    } else {
      tmpl <- if (flock_df$is_bat[i]) bat_poly else bird_poly
      p <- p + geom_polygon(data = data.frame(x = cx[i] + tmpl$px * bird_w[i], y = yc + tmpl$py * bird_h[i]), aes(x = x, y = y), fill = ifelse(flock_df$is_bat[i], "#4a4e5c", "#7d6b57"), color = "#2b241d", linewidth = 0.45)
    }
  }
  
  p <- p + geom_text(data = flock_df, aes(x = x_wide, y = y - (bird_h / 2) - 0.05, label = label_txt), size = txt_size, family = "serif", color = "#3b3228", vjust = 1)
  ggsave(png_file, plot = p, width = 13, height = 8, dpi = 220, bg = "#f6f0e2")
  
  W <- 13; H <- 8; p_left <- 0.22; p_right <- W - 0.22; p_top <- 0.92; p_bot <- H - 0.19
  sc <- min((p_right - p_left) / (x_max - x_min), (p_bot - p_top) / (y_max - y_min))
  ox <- p_left + ((p_right - p_left) - (x_max - x_min) * sc) / 2
  oy <- p_top + ((p_bot - p_top) - (y_max - y_min) * sc) / 2
  
  inat_base <- paste0("https", "://", paste(c("www", "inaturalist", "org"), collapse = "."), "/taxa/search?q=")
  ebird_sp  <- paste0("https", "://", paste(c("ebird", "org"), collapse = "."), "/species/")
  ebird_ex  <- paste0("https", "://", paste(c("ebird", "org"), collapse = "."), "/explore?q=")
  
  hotspots <- data.frame(
    common_name = flock_df$common_name, detections = flock_df$detections,
    left_pct   = round(((ox + ((cx - box_hw) - x_min) * sc) / W) * 100, 3),
    top_pct    = round(((oy + (y_max - (cy + box_hh)) * sc) / H) * 100, 3),
    width_pct  = round(((2 * box_hw * sc) / W) * 100, 3),
    height_pct = round(((2 * box_hh * sc) / H) * 100, 3),
    site_name  = ifelse(flock_df$is_bat, "iNaturalist", "eBird"),
    target_url = ifelse(flock_df$is_bat, paste0(inat_base, URLencode(flock_df$clean_sci)), ifelse(nchar(flock_df$ebird_code) > 0, paste0(ebird_sp, flock_df$ebird_code), paste0(ebird_ex, URLencode(flock_df$common_name)))),
    stringsAsFactors = FALSE
  )
  write_json(hotspots, json_file, pretty = TRUE, auto_unbox = TRUE)
}

cutoff_48h <- format(as.POSIXct(Sys.time() - as.difftime(48, units = "hours"), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
build_collage(cutoff_48h, "Last 48 Hours", "HeckfordBirds.png", "species_map.json", compact_mode = FALSE)

current_year <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%Y")
cutoff_ytd <- paste0(current_year, "-01-01T00:00:00Z")
build_collage(cutoff_ytd, paste0("Year to Date (", current_year, ")"), "HeckfordBirds_YTD.png", "species_map_ytd.json", compact_mode = TRUE)
