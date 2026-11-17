if (!requireNamespace("png", quietly = TRUE)) {
  install.packages("png", repos = "https://p3m.dev")
}

library(httr)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(png)
library(grid)

station_token <- Sys.getenv("BIRDWEATHER_TOKEN")
stopifnot("Missing BIRDWEATHER_TOKEN secret" = nchar(station_token) > 0)

# 1. Fetch ALL detections from the last 24 hours (with pagination)
api_url <- paste0(
  "https://", paste(c("app", "birdweather", "com"), collapse = "."), "/",
  paste(c("api", "v1", "stations", station_token, "detections"), collapse = "/")
)

cutoff_utc <- Sys.time() - as.difftime(24, units = "hours")
all_pages <- list()
cursor_id <- NULL

repeat {
  q <- list(limit = 500, period = "day")
  if (!is.null(cursor_id)) q$before <- cursor_id
  
  res <- GET(
    url = api_url,
    query = q,
    add_headers("Accept" = "application/json", "User-Agent" = "Mozilla/5.0")
  )
  stop_for_status(res)
  parsed <- fromJSON(content(res, as = "text", encoding = "UTF-8"), flatten = TRUE)
  batch <- parsed$detections
  
  if (is.null(batch) || nrow(batch) == 0) break
  all_pages[[length(all_pages) + 1]] <- batch
  
  ts_col <- grep("^timestamp$", names(batch), value = TRUE)
  if (length(ts_col) > 0) {
    oldest_ts <- suppressWarnings(as.POSIXct(min(batch[[ts_col]], na.rm = TRUE), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"))
    if (!is.na(oldest_ts) && oldest_ts < cutoff_utc) break
  }
  if (nrow(batch) < 500 || !"id" %in% names(batch)) break
  cursor_id <- tail(batch$id, 1)
}

detections_df <- bind_rows(all_pages)

if ("timestamp" %in% names(detections_df)) {
  ts_parsed <- suppressWarnings(as.POSIXct(detections_df$timestamp, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"))
  if (any(!is.na(ts_parsed))) {
    detections_df <- detections_df[!is.na(ts_parsed) & ts_parsed >= cutoff_utc, ]
  }
}

# 2. Summarize ALL detected species & apply Avian Visitors scaling (count^0.65)
common_col <- grep("commonName|common_name", names(detections_df), value = TRUE)
sci_col    <- grep("scientificName|scientific_name", names(detections_df), value = TRUE)

flock_df <- detections_df %>%
  rename(common_name = all_of(common_col), sci_name = all_of(sci_col)) %>%
  filter(!is.na(common_name), !is.na(sci_name)) %>%
  mutate(clean_sci = sub("^(\\S+\\s+\\S+).*", "\\1", sci_name)) %>%
  count(common_name, clean_sci, name = "detections") %>%
  arrange(desc(detections)) %>%
  mutate(
    pose = ifelse(row_number() %% 2 == 1, 1, 2),
    area_score = detections^0.65
  )

# 3. Fetch full-color woodblock bird illustrations from Avian Visitors API
get_bird_art <- function(sci_name, common_name, pose = 1) {
  cutout_url <- paste0(
    "https://", paste(c("bird", "onethreenine", "net"), collapse = "."),
    "/avian/api/cutout.php"
  )
  res <- try(GET(
    url = cutout_url,
    query = list(sci = sci_name, com = common_name, pose = pose),
    add_headers("User-Agent" = "Mozilla/5.0")
  ), silent = TRUE)
  if (!inherits(res, "try-error") && status_code(res) == 200) {
    img <- try(readPNG(content(res, as = "raw")), silent = TRUE)
    if (!inherits(img, "try-error")) return(img)
  }
  return(NULL)
}

rasters <- mapply(get_bird_art, flock_df$clean_sci, flock_df$common_name, flock_df$pose, SIMPLIFY = FALSE)
valid <- !sapply(rasters, is.null)
flock_df <- flock_df[valid, ]
rasters <- rasters[valid]

# 4. Built-in center-out spiral circle packing (replaces external packcircles package)
pack_spiral <- function(areas) {
  r <- sqrt(areas / pi)
  n <- length(r)
  x <- numeric(n)
  y <- numeric(n)
  if (n > 1) {
    angles <- seq(0, 25 * pi, length.out = 2500)
    for (i in 2:n) {
      step <- r[1] * 0.04
      rad <- 0
      for (ang in angles) {
        cx <- rad * cos(ang) * 1.5 # 1.5x horizontal landscape bias
        cy <- rad * sin(ang)
        dists <- sqrt(((cx - x[1:(i - 1)]) / 1.5)^2 + (cy - y[1:(i - 1)])^2)
        if (all(dists >= (r[i] + r[1:(i - 1)]))) {
          x[i] <- cx
          y[i] <- cy
          break
        }
        rad <- rad + step
      }
    }
  }
  tibble(x_wide = x, y = y, radius = r)
}

packing <- pack_spiral(flock_df$area_score * 1.35)
flock_df <- bind_cols(flock_df, packing)

# 5. Render warm cream woodblock collage
bc_time <- format(as.POSIXct(Sys.time(), tz = "America/Vancouver"), "%b %d, %I:%M %p %Z")
total_calls <- sum(flock_df$detections)

p <- ggplot(flock_df, aes(x = x_wide, y = y)) +
  geom_blank() +
  coord_equal(expand = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.22)) +
  scale_y_continuous(expand = expansion(mult = 0.22)) +
  labs(
    title = "AVIAN VISITORS — ARBUTUS (#29794)",
    subtitle = paste0("Last 24 Hours: ", nrow(flock_df), " Species (", total_calls, " Detections) • Updated ", bc_time)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#f6f0e2", color = NA),
    panel.background = element_rect(fill = "#f6f0e2", color = NA),
    plot.title = element_text(
      hjust = 0.5, family = "serif", size = 18, face = "bold",
      color = "#2b241d", margin = margin(t = 15, b = 3)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, family = "serif", size = 11, color = "#6e6051", margin = margin(b = 8)
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

for (i in seq_len(nrow(flock_df))) {
  img <- rasters[[i]]
  h_px <- nrow(img)
  w_px <- ncol(img)
  scale_factor <- (flock_df$radius[i] * 1.55) / sqrt(h_px^2 + w_px^2)
  half_w <- (w_px * scale_factor) / 2
  half_h <- (h_px * scale_factor) / 2
  
  p <- p + annotation_raster(
    img,
    xmin = flock_df$x_wide[i] - half_w,
    xmax = flock_df$x_wide[i] + half_w,
    ymin = flock_df$y[i] - half_h,
    ymax = flock_df$y[i] + half_h
  )
}

p <- p + geom_text(
  data = flock_df,
  aes(x = x_wide, y = y - radius * 0.88, label = paste0(common_name, " (", detections, ")")),
  size = 2.9, family = "serif", color = "#3b3228", check_overlap = TRUE
)

ggsave("avian_visitors.png", plot = p, width = 11, height = 7.5, dpi = 200, bg = "#f6f0e2")
