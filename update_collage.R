# ==============================================================================
# SECTION 1: LIBRARIES, STATION TOKENS & BIRDWEATHER API FETCHERS
# ==============================================================================
.libPaths(c(.libPaths(), "/usr/lib/R/site-library"))
library(httr); library(jsonlite); library(dplyr); library(ggplot2); library(png); library(grid)

# Load BirdWeather API tokens from GitHub Secrets (supports 1 or 2 stations)
tokens <- c(Sys.getenv("BIRDWEATHER_TOKEN"), Sys.getenv("BIRDWEATHER_TOKEN_2"))
tokens <- tokens[nchar(tokens) > 0]
stopifnot("Missing BIRDWEATHER_TOKEN secret" = length(tokens) > 0)

# Helper: parse HTTP response only if it is valid JSON (guards against HTML 200 pages)
safe_from_json <- function(res) {
  if (is.null(res) || inherits(res, "try-error") || status_code(res) != 200) return(NULL)
  txt <- trimws(content(res, as = "text", encoding = "UTF-8"))
  if (nchar(txt) == 0 || grepl("^<", txt)) return(NULL)
  parsed <- try(fromJSON(txt, flatten = TRUE), silent = TRUE)
  if (inherits(parsed, "try-error")) return(NULL)
  parsed
}

# Fetch full species summary table (used for "All" confidence level so 0 species drop out)
fetch_species_summary <- function(tok, from_iso) {
  empty_res <- tibble(common_name = character(), sci_name = character(), detections = numeric(), tax_grp = character())
  sp_url <- paste0("https", "://", paste(c("app", "birdweather", "com"), collapse = "."), "/", paste(c("api", "v1", "stations", tok, "species"), collapse = "/"))
  all_pages <- list(); pg <- 0
  repeat {
    pg <- pg + 1
    parsed <- NULL
    for (att in 1:3) {
      r_try <- try(GET(url = sp_url, query = list(from = from_iso, limit = 100, page = pg), add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0"), timeout(20)), silent = TRUE)
      parsed <- safe_from_json(r_try)
      if (!is.null(parsed)) break
      Sys.sleep(0.5 * att)
    }
    if (is.null(parsed)) break
    sp_df <- if (!is.null(parsed$species)) parsed$species else parsed
    if (!is.data.frame(sp_df) || nrow(sp_df) == 0) break
    c_col <- grep("commonName|common_name", names(sp_df), value = TRUE)[1]
    s_col <- grep("scientificName|scientific_name", names(sp_df), value = TRUE)[1]
    n_col <- grep("detections_count|detections\\.total|count|total", names(sp_df), value = TRUE)[1]
    g_col <- grep("group|taxon|class|category", names(sp_df), value = TRUE)[1]
    if (is.na(c_col) || is.na(s_col) || is.na(n_col)) break
    grp_vals <- if (!is.na(g_col)) as.character(sp_df[[g_col]]) else rep("", nrow(sp_df))
    df <- tibble(common_name = as.character(sp_df[[c_col]]), sci_name = as.character(sp_df[[s_col]]), detections = as.numeric(sp_df[[n_col]]), tax_grp = grp_vals)
    df <- df[!is.na(df$common_name) & !is.na(df$sci_name) & !is.na(df$detections) & df$detections > 0, ]
    all_pages[[pg]] <- df
    if (nrow(sp_df) < 100) break
  }
  if (length(all_pages) == 0) return(empty_res)
  bind_rows(all_pages)
}

# Fetch individual timestamped detections >= 60% confidence (used for 60/70/80/90% toggles)
fetch_all_detections <- function(tok, from_iso, min_conf = 0.60) {
  empty_res <- tibble(common_name = character(), sci_name = character(), conf = numeric(), ts = as.POSIXct(character()))
  det_url <- paste0("https", "://", paste(c("app", "birdweather", "com"), collapse = "."), "/", paste(c("api", "v1", "stations", tok, "detections"), collapse = "/"))
  all_rows <- list(); cursor_id <- NULL; pg <- 0
  repeat {
    pg <- pg + 1
    q_list <- list(from = from_iso, confidence = min_conf, limit = 100, order = "desc")
    if (!is.null(cursor_id)) q_list$cursor <- cursor_id
    parsed <- NULL
    for (att in 1:3) {
      r_try <- try(GET(url = det_url, query = q_list, add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0"), timeout(25)), silent = TRUE)
      parsed <- safe_from_json(r_try)
      if (!is.null(parsed)) break
      Sys.sleep(0.5 * att)
    }
    if (is.null(parsed)) break
    d_df <- if (!is.null(parsed$detections)) parsed$detections else parsed
    if (!is.data.frame(d_df) || nrow(d_df) == 0) break
    c_col <- grep("species\\.commonName|commonName|common_name", names(d_df), value = TRUE)[1]
    s_col <- grep("species\\.scientificName|scientificName|scientific_name", names(d_df), value = TRUE)[1]
    t_col <- grep("^timestamp$|createdAt|created_at", names(d_df), value = TRUE)[1]
    conf_col <- grep("^confidence$", names(d_df), value = TRUE)[1]
    if (is.na(c_col) || is.na(s_col)) break
    ts_vals <- if (!is.na(t_col)) as.POSIXct(sub("Z$", "", d_df[[t_col]]), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") else rep(Sys.time(), nrow(d_df))
    conf_vals <- if (!is.na(conf_col)) as.numeric(d_df[[conf_col]]) else rep(1, nrow(d_df))
    sub_df <- tibble(common_name = as.character(d_df[[c_col]]), sci_name = as.character(d_df[[s_col]]), conf = conf_vals, ts = ts_vals)
    sub_df <- sub_df[!is.na(sub_df$common_name) & !is.na(sub_df$sci_name) & !is.na(sub_df$conf) & sub_df$conf >= min_conf, ]
    all_rows[[pg]] <- sub_df
    if (nrow(d_df) < 100 || !"id" %in% names(d_df)) break
    new_cursor <- tail(d_df$id, 1)
    if (!is.null(cursor_id) && identical(new_cursor, cursor_id)) break
    cursor_id <- new_cursor
    if (pg %% 25 == 0) Sys.sleep(0.2)
  }
  if (length(all_rows) == 0) return(empty_res)
  bind_rows(all_rows)
}

# ==============================================================================
# SECTION 2: EBIRD TAXONOMY, WOODBLOCK ART CACHE & FALLBACK SILHOUETTES
# ==============================================================================
# Download official eBird taxonomy once per run (provides speciesCode links & taxonomic order)
ebird_tax_url <- paste0("https", "://", paste(c("api", "ebird", "org"), collapse = "."), "/v2/ref/taxonomy/ebird")
tax_res <- try(GET(url = ebird_tax_url, query = list(fmt = "json", cat = "species"), add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0"), timeout(15)), silent = TRUE)
ebird_lookup <- data.frame(sciName = character(), comName = character(), speciesCode = character(), order = character(), stringsAsFactors = FALSE)
tax_parsed <- safe_from_json(tax_res)
if (is.data.frame(tax_parsed) && nrow(tax_parsed) > 0) ebird_lookup <- tax_parsed
)

# Cache woodblock bird illustrations in memory so each species is downloaded only once
art_cache <- new.env(parent = emptyenv())
get_bird_art <- function(sci_name, common_name, pose = 1, is_bat = FALSE) {
  if (is_bat) return(NULL)
  ck <- paste0(sci_name, "_", pose)
  if (exists(ck, envir = art_cache)) return(get(ck, envir = art_cache))
  cutout_url <- paste0("https", "://", paste(c("bird", "onethreenine", "net"), collapse = "."), "/avian/api/cutout.php")
  for (p_try in unique(c(pose, 1, 2))) {
    for (attempt in 1:2) {
      Sys.sleep(0.10)
      res <- try(GET(url = cutout_url, query = list(sci = sci_name, com = common_name, pose = p_try), add_headers("User-Agent" = "Mozilla/5.0"), timeout(15)), silent = TRUE)
      if (!inherits(res, "try-error") && status_code(res) == 200) {
        img <- try(readPNG(content(res, as = "raw")), silent = TRUE)
        if (!inherits(img, "try-error")) { assign(ck, img, envir = art_cache); return(img) }
      }
      Sys.sleep(0.25)
    }
  }
  assign(ck, NULL, envir = art_cache)
  return(NULL)
}

# Vector silhouettes used for bats and any bird lacking a woodblock PNG
bird_poly <- data.frame(px = c(-0.48,-0.25,-0.08,0.08,0.22,0.34,0.48,0.35,0.32,0.24,0.10,0.06,0.12,0.02,-0.04,-0.02,-0.16,-0.45,-0.36), py = c(-0.28,-0.08,0.14,0.22,0.38,0.34,0.28,0.22,0.08,-0.10,-0.22,-0.44,-0.46,-0.46,-0.24,-0.45,-0.22,-0.35,-0.26))
bat_poly  <- data.frame(px = c(0.00,-0.04,-0.07,-0.14,-0.32,-0.50,-0.38,-0.26,-0.14,-0.06,0.00,0.06,0.14,0.26,0.38,0.50,0.32,0.14,0.07,0.04), py = c(0.18,0.34,0.22,0.32,0.42,0.12,-0.12,-0.05,-0.22,-0.15,-0.35,-0.15,-0.22,-0.05,-0.12,0.12,0.42,0.32,0.22,0.34))
# ==============================================================================
# SECTION 3: TAXONOMIC GROUPING & HORIZONTAL FIELD-GUIDE ROW PACKER
# ==============================================================================
# Assign every species into a deterministic group (catch-all ensures 0 species drop out)
classify_species_group <- function(common_name, clean_sci, tax_grp = "", ord_str = "") {
  is_bat <- grepl("bat|myotis|flying fox|pipistrelle", common_name, ignore.case = TRUE) |
    grepl("^(Myotis|Aeorestes|Lasionycteris|Eptesicus|Lasiurus|Corynorhinus|Antrozous|Tadarida|Perimyotis|Parastrellus|Nyctinomops|Eumops|Chiroptera)\\b", clean_sci, ignore.case = TRUE) |
    grepl("bat|chiroptera|mammalia", tax_grp, ignore.case = TRUE)
  if (is_bat) return("Bats")
  if (ord_str == "Passeriformes") return("Songbirds")
  if (ord_str %in% c("Accipitriformes", "Falconiformes", "Strigiformes")) return("Raptors")
  if (ord_str %in% c("Anseriformes", "Podicipediformes", "Gaviiformes", "Procellariiformes", "Suliformes", "Pelecaniformes")) return("Waterfowl & Seabirds")
  if (ord_str == "Charadriiformes") return("Shorebirds & Gulls")
  return("Woodpeckers & Other Birds")
}

# Arrange species left-to-right in ranked rows inside each group box (+1.25 vertical padding)
pack_group_cluster <- function(sub_idx, box_hw, box_hh, aspect_bias = 1.15) {
  k_n <- length(sub_idx)
  lx <- numeric(k_n); ly <- numeric(k_n)
  widths <- box_hw[sub_idx] * 2; heights <- box_hh[sub_idx] * 2
  target_w <- max(max(widths) + 0.4, sqrt(sum((widths + 0.26) * (heights + 0.24)) * aspect_bias) * 1.15)
  cur_x <- 0; cur_y <- 0; row_h <- 0; row_ids <- integer(0)
  for (k in seq_len(k_n)) {
    if (length(row_ids) > 0 && (cur_x + widths[k] > target_w)) {
      row_w <- cur_x - 0.26
      lx[row_ids] <- lx[row_ids] - (row_w / 2); ly[row_ids] <- cur_y - (row_h / 2)
      cur_y <- cur_y - row_h - 0.24; cur_x <- 0; row_h <- 0; row_ids <- integer(0)
    }
    lx[k] <- cur_x + widths[k] / 2; row_h <- max(row_h, heights[k]); row_ids <- c(row_ids, k); cur_x <- cur_x + widths[k] + 0.26
  }
  if (length(row_ids) > 0) {
    row_w <- cur_x - 0.26
    lx[row_ids] <- lx[row_ids] - (row_w / 2); ly[row_ids] <- cur_y - (row_h / 2)
  }
  hw_span <- max(lx + box_hw[sub_idx]) - min(lx - box_hw[sub_idx])
  hh_span <- max(ly + box_hh[sub_idx]) - min(ly - box_hh[sub_idx])
  x_mid <- (max(lx + box_hw[sub_idx]) + min(lx - box_hw[sub_idx])) / 2
  y_mid <- (max(ly + box_hh[sub_idx]) + min(ly - box_hh[sub_idx])) / 2
  list(lx = lx - x_mid, ly = ly - y_mid, w = hw_span + 0.52, h = hh_span + 1.25)
}
# ==============================================================================
# SECTION 4A: TAXONOMY MATCHING, SIZING & GROUP PANEL POSITIONING
# ==============================================================================
render_flock_df <- function(flock_df, period_label, conf_label, png_file, json_file, compact_mode = FALSE) {
  bc_time <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%b %d, %I:%M %p %Z")
  if (nrow(flock_df) == 0) {
    p_empty <- ggplot() + geom_text(aes(x = 0, y = 0, label = paste0("No detections (", conf_label, ") in ", period_label)), family = "serif", size = 6, color = "#5c4e3f") + theme_void() + theme(plot.background = element_rect(fill = "#f6f0e2", color = NA))
    ggsave(png_file, plot = p_empty, width = 13, height = 8, dpi = 220, bg = "#f6f0e2")
    write_json(list(), json_file, pretty = TRUE, auto_unbox = TRUE)
    return(invisible(NULL))
  }
  if (!"tax_grp" %in% names(flock_df)) flock_df$tax_grp <- ""
  flock_df <- flock_df[order(-flock_df$detections), ]; rownames(flock_df) <- NULL
  flock_df$ebird_code <- character(nrow(flock_df)); flock_df$ebird_order <- character(nrow(flock_df))
  for (idx in seq_len(nrow(flock_df))) {
    m1 <- which(ebird_lookup$sciName == flock_df$clean_sci[idx])
    if (length(m1) == 0) m1 <- which(tolower(ebird_lookup$comName) == tolower(flock_df$common_name[idx]))
    if (length(m1) > 0) { flock_df$ebird_code[idx] <- ebird_lookup$speciesCode[m1[1]]; flock_df$ebird_order[idx] <- ebird_lookup$order[m1[1]] }
  }
  flock_df$group_name <- mapply(classify_species_group, flock_df$common_name, flock_df$clean_sci, flock_df$tax_grp, flock_df$ebird_order)
  flock_df$is_bat <- flock_df$group_name == "Bats"
  flock_df$pose <- ifelse(seq_len(nrow(flock_df)) %% 2 == 1, 1, 2)
  flock_df$scale_r <- flock_df$detections^ifelse(compact_mode, 0.25, 0.32)
  r_min <- ifelse(compact_mode, 0.60, 0.88); r_span <- ifelse(compact_mode, 0.85, 1.35)
  flock_df$radius <- r_min + r_span * ((flock_df$scale_r - min(flock_df$scale_r)) / max(1e-6, max(flock_df$scale_r) - min(flock_df$scale_r)))
  
  # Two-line stacked labels: Species Name on line 1, formatted (count) on line 2
  count_txt <- paste0("(", format(flock_df$detections, big.mark = ",", trim = TRUE), ")")
  flock_df$label_txt <- paste0(flock_df$common_name, "\n", count_txt)
  max_chars <- pmax(nchar(flock_df$common_name), nchar(count_txt))
  rasters <- mapply(get_bird_art, flock_df$clean_sci, flock_df$common_name, flock_df$pose, flock_df$is_bat, SIMPLIFY = FALSE)
  has_art <- !sapply(rasters, is.null)
  
  # Compute bounding boxes for each species icon + 2-line label
  n <- nrow(flock_df); bird_w <- numeric(n); bird_h <- numeric(n); box_hw <- numeric(n); box_hh <- numeric(n)
  label_band <- ifelse(compact_mode, 0.56, 0.70); char_w <- ifelse(compact_mode, 0.16, 0.20)
  for (i in seq_len(n)) {
    if (has_art[i]) {
      sf <- (flock_df$radius[i] * 1.95) / sqrt(nrow(rasters[[i]])^2 + ncol(rasters[[i]])^2)
      bird_w[i] <- ncol(rasters[[i]]) * sf; bird_h[i] <- nrow(rasters[[i]]) * sf
    } else {
      bird_w[i] <- flock_df$radius[i] * ifelse(flock_df$is_bat[i], 1.55, 1.25)
      bird_h[i] <- flock_df$radius[i] * ifelse(flock_df$is_bat[i], 0.85, 1.05)
    }
    box_hw[i] <- max(bird_w[i], max_chars[i] * char_w) / 2 + 0.14
    box_hh[i] <- (bird_h[i] + label_band) / 2 + 0.10
  }

  # ============================================================================
  # SECTION 4B: PANEL PLACEMENT, GGPLOT DRAWING & HOTSPOT JSON EXPORT
  # ============================================================================
  cx <- numeric(n); cy <- numeric(n); panel_boxes <- list()
  song_idx <- which(flock_df$group_name == "Songbirds")
  right_groups <- c("Raptors", "Waterfowl & Seabirds", "Shorebirds & Gulls", "Woodpeckers & Other Birds", "Bats")
  active_right <- right_groups[sapply(right_groups, function(g) any(flock_df$group_name == g))]
  song_pack <- if (length(song_idx) > 0) pack_group_cluster(song_idx, box_hw, box_hh, aspect_bias = 1.05) else NULL
  right_packs <- lapply(active_right, function(g) { idx <- which(flock_df$group_name == g); pk <- pack_group_cluster(idx, box_hw, box_hh, aspect_bias = 1.25); pk$idx <- idx; pk$name <- g; pk })
  r_cols <- if (length(right_packs) > 2) 2 else 1
  col_widths <- numeric(r_cols); col_heights <- numeric(r_cols)
  for (m in seq_along(right_packs)) {
    c_id <- ((m - 1) %% r_cols) + 1
    col_widths[c_id] <- max(col_widths[c_id], right_packs[[m]]$w)
    col_heights[c_id] <- col_heights[c_id] + right_packs[[m]]$h + ifelse(col_heights[c_id] > 0, 0.25, 0)
  }
  total_h <- max(ifelse(!is.null(song_pack), song_pack$h, 0), max(c(0, col_heights)))
  cur_x <- 0
  if (!is.null(song_pack)) {
    s_top <- total_h / 2; s_bot <- s_top - max(song_pack$h, total_h)
    s_xmid <- cur_x + song_pack$w / 2; s_ymid <- (s_top + s_bot) / 2 - 0.16
    cx[song_idx] <- s_xmid + song_pack$lx; cy[song_idx] <- s_ymid + song_pack$ly
    panel_boxes[[length(panel_boxes) + 1]] <- list(xmin = cur_x, xmax = cur_x + song_pack$w, ymin = s_bot, ymax = s_top, title = "Songbirds")
    cur_x <- cur_x + song_pack$w + 0.28
  }
  if (length(right_packs) > 0) {
    col_x_starts <- cumsum(c(cur_x, head(col_widths, -1) + 0.25))
    col_y_curs <- rep(total_h / 2, r_cols)
    for (m in seq_along(right_packs)) {
      c_id <- ((m - 1) %% r_cols) + 1; pk <- right_packs[[m]]
      p_xmin <- col_x_starts[c_id]; p_xmax <- p_xmin + col_widths[c_id]
      p_ymax <- col_y_curs[c_id]; p_ymin <- p_ymax - pk$h
      cx[pk$idx] <- (p_xmin + p_xmax) / 2 + pk$lx; cy[pk$idx] <- (p_ymin + p_ymax) / 2 - 0.16 + pk$ly
      panel_boxes[[length(panel_boxes) + 1]] <- list(xmin = p_xmin, xmax = p_xmax, ymin = p_ymin, ymax = p_ymax, title = pk$name)
      col_y_curs[c_id] <- p_ymin - 0.25
    }
  }
  
  # Render ggplot canvas, panel boxes, illustrations, labels, and JSON clickable map
  flock_df$x_wide <- cx; flock_df$y <- cy
  x_min <- min(sapply(panel_boxes, `[[`, "xmin")) - 0.25; x_max <- max(sapply(panel_boxes, `[[`, "xmax")) + 0.25
  y_min <- min(sapply(panel_boxes, `[[`, "ymin")) - 0.25; y_max <- max(sapply(panel_boxes, `[[`, "ymax")) + 0.25
  bird_cnt <- sum(!flock_df$is_bat); bat_cnt <- sum(flock_df$is_bat)
  p <- ggplot(flock_df, aes(x = x_wide, y = y)) + geom_blank() + coord_equal(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE, clip = "off") +
    labs(title = "Qathet: Myrtle Rock Bird and Bat Detections", subtitle = paste0(period_label, " • ", conf_label, " (", length(tokens), " Stations): ", bird_cnt, " Bird & ", bat_cnt, " Bat Species (", sum(flock_df$detections), " Detections) • Updated ", bc_time)) +
    theme_void() + theme(plot.background = element_rect(fill = "#f6f0e2", color = NA), panel.background = element_rect(fill = "#f6f0e2", color = NA), plot.title = element_text(hjust = 0.5, family = "serif", size = 20, face = "bold", color = "#2b241d", margin = margin(t = 12, b = 2)), plot.subtitle = element_text(hjust = 0.5, family = "serif", size = 11.5, color = "#6e6051", margin = margin(b = 8)), plot.margin = margin(12, 16, 14, 16))
  for (pb in panel_boxes) {
    box_fill <- ifelse(pb$title == "Bats", "#ebe1ce", "#f1eadb")
    p <- p + annotate("rect", xmin = pb$xmin, xmax = pb$xmax, ymin = pb$ymin, ymax = pb$ymax, fill = box_fill, color = "#cbbda6", linewidth = 0.55) +
      annotate("text", x = (pb$xmin + pb$xmax) / 2, y = pb$ymax - 0.34, label = pb$title, family = "serif", fontface = "bold", size = 3.3, color = "#3b3228")
  }
  for (i in seq_len(n)) {
    yc <- cy[i] + (label_band / 2)
    if (has_art[i]) { p <- p + annotation_raster(rasters[[i]], xmin = cx[i] - bird_w[i]/2, xmax = cx[i] + bird_w[i]/2, ymin = yc - bird_h[i]/2, ymax = yc + bird_h[i]/2) }
    else { tmpl <- if (flock_df$is_bat[i]) bat_poly else bird_poly; p <- p + geom_polygon(data = data.frame(x = cx[i] + tmpl$px * bird_w[i], y = yc + tmpl$py * bird_h[i]), aes(x = x, y = y), fill = ifelse(flock_df$is_bat[i], "#4a4e5c", "#7d6b57"), color = "#2b241d", linewidth = 0.45) }
  }
  p <- p + geom_text(data = flock_df, aes(x = x_wide, y = y - (bird_h / 2) - 0.05, label = label_txt), size = ifelse(compact_mode, 2.25, 2.9), lineheight = 0.92, family = "serif", color = "#3b3228", vjust = 1)
  ggsave(png_file, plot = p, width = 13, height = 8, dpi = 220, bg = "#f6f0e2")
  W <- 13; H <- 8; p_left <- 0.22; p_right <- W - 0.22; p_top <- 0.92; p_bot <- H - 0.19
  sc <- min((p_right - p_left) / (x_max - x_min), (p_bot - p_top) / (y_max - y_min))
  ox <- p_left + ((p_right - p_left) - (x_max - x_min) * sc) / 2; oy <- p_top + ((p_bot - p_top) - (y_max - y_min) * sc) / 2
  inat_base <- paste0("https", "://", paste(c("www", "inaturalist", "org"), collapse = "."), "/taxa/search?q=")
  ebird_sp  <- paste0("https", "://", paste(c("ebird", "org"), collapse = "."), "/species/")
  ebird_ex  <- paste0("https", "://", paste(c("ebird", "org"), collapse = "."), "/explore?q=")
  hotspots <- data.frame(common_name = flock_df$common_name, detections = flock_df$detections, left_pct = round(((ox + ((cx - box_hw) - x_min) * sc) / W) * 100, 3), top_pct = round(((oy + (y_max - (cy + box_hh)) * sc) / H) * 100, 3), width_pct = round(((2 * box_hw * sc) / W) * 100, 3), height_pct = round(((2 * box_hh * sc) / H) * 100, 3), site_name = ifelse(flock_df$is_bat, "iNaturalist", "eBird"), target_url = ifelse(flock_df$is_bat, paste0(inat_base, URLencode(flock_df$clean_sci, reserved = TRUE)), ifelse(nchar(flock_df$ebird_code) > 0, paste0(ebird_sp, flock_df$ebird_code), paste0(ebird_ex, URLencode(flock_df$common_name, reserved = TRUE)))), stringsAsFactors = FALSE)
  write_json(hotspots, json_file, pretty = TRUE, auto_unbox = TRUE)
}
# ==============================================================================
# SECTION 5: TIMEFRAME CUTOFFS & EXECUTION LOOP (ALL, 60%, 70%, 80%, 90%)
# ==============================================================================
current_year <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%Y")
cutoff_ytd <- paste0(current_year, "-01-01T00:00:00Z")
cutoff_48h_iso <- format(as.POSIXct(Sys.time() - as.difftime(48, units = "hours"), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
cutoff_48h_utc <- as.POSIXct(Sys.time() - as.difftime(48, units = "hours"), tz = "UTC")

# 1. Generate "All Detections" (c = 0) from /species endpoint so 100% of species are included
sp_48h_raw <- bind_rows(lapply(tokens, fetch_species_summary, from_iso = cutoff_48h_iso))
sp_48h_all <- if (nrow(sp_48h_raw) > 0) {
  sp_48h_raw$clean_sci <- sub("^(\\S+\\s+\\S+).*", "\\1", sp_48h_raw$sci_name)
  aggregate(detections ~ common_name + clean_sci + tax_grp, data = sp_48h_raw, FUN = sum)
} else {
  data.frame(common_name = character(), clean_sci = character(), tax_grp = character(), detections = numeric())
}
render_flock_df(sp_48h_all, "Last 48 Hours", "All Detections", "HeckfordBirds_48h_0.png", "species_map_48h_0.json", compact_mode = FALSE)

sp_ytd_raw <- bind_rows(lapply(tokens, fetch_species_summary, from_iso = cutoff_ytd))
sp_ytd_all <- if (nrow(sp_ytd_raw) > 0) {
  sp_ytd_raw$clean_sci <- sub("^(\\S+\\s+\\S+).*", "\\1", sp_ytd_raw$sci_name)
  aggregate(detections ~ common_name + clean_sci + tax_grp, data = sp_ytd_raw, FUN = sum)
} else {
  data.frame(common_name = character(), clean_sci = character(), tax_grp = character(), detections = numeric())
}
render_flock_df(sp_ytd_all, paste0("Year to Date (", current_year, ")"), "All Detections", "HeckfordBirds_ytd_0.png", "species_map_ytd_0.json", compact_mode = TRUE)

# 2. Fetch individual detections >= 60% once, then filter 60%, 70%, 80%, 90% locally
ytd_conf <- bind_rows(lapply(tokens, fetch_all_detections, from_iso = cutoff_ytd, min_conf = 0.60))
det_48h_conf <- ytd_conf[is.na(ytd_conf$ts) | ytd_conf$ts >= cutoff_48h_utc, ]

for (c_pct in c(60, 70, 80, 90)) {
  sub_48 <- det_48h_conf[det_48h_conf$conf >= (c_pct / 100), ]
  df_48 <- if (nrow(sub_48) > 0) { sub_48$clean_sci <- sub("^(\\S+\\s+\\S+).*", "\\1", sub_48$sci_name); a <- aggregate(conf ~ common_name + clean_sci, data = sub_48, FUN = length); colnames(a) <- c("common_name", "clean_sci", "detections"); a } else data.frame(common_name = character(), clean_sci = character(), detections = numeric())
  render_flock_df(df_48, "Last 48 Hours", paste0("≥", c_pct, "% Confidence"), paste0("HeckfordBirds_48h_", c_pct, ".png"), paste0("species_map_48h_", c_pct, ".json"), compact_mode = FALSE)
  
  sub_ytd <- ytd_conf[ytd_conf$conf >= (c_pct / 100), ]
  df_ytd <- if (nrow(sub_ytd) > 0) { sub_ytd$clean_sci <- sub("^(\\S+\\s+\\S+).*", "\\1", sub_ytd$sci_name); a <- aggregate(conf ~ common_name + clean_sci, data = sub_ytd, FUN = length); colnames(a) <- c("common_name", "clean_sci", "detections"); a } else data.frame(common_name = character(), clean_sci = character(), detections = numeric())
  render_flock_df(df_ytd, paste0("Year to Date (", current_year, ")"), paste0("≥", c_pct, "% Confidence"), paste0("HeckfordBirds_ytd_", c_pct, ".png"), paste0("species_map_ytd_", c_pct, ".json"), compact_mode = TRUE)
}
